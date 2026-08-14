"""
Client for the Copernicus Climate Data Store (CDS) **Retrieve API** — the
job-submission service used to order data that is not published as a
cloud-optimized store.

This is a different protocol from the rest of this package's CDS access. ERA5-Land
is read directly from ECMWF's ARCO Zarr host by HTTP byte range
(`authenticated_http_store.jl`), which needs no job. Datasets that exist only as
archived files (e.g. `satellite-albedo`) must instead be *ordered*:

    POST /processes/{dataset}/execution   → job id
    GET  /jobs/{job_id}                   → poll until status == "successful"
    GET  /jobs/{job_id}/results           → asset href
    download the asset                    → a ZIP of NetCDF files

The service implements the OGC API - Processes standard; endpoint shapes here were
verified against the live spec at
https://cds.climate.copernicus.eu/api/retrieve/v1/openapi.json

Because a job is queued server-side, calls can take minutes. Requests are also
size-limited (see [`cds_estimate_cost`](@ref)), so callers are expected to split
large requests into several jobs.

!!! note "Authentication differs from the Zarr path"
    The Retrieve API authenticates with a `PRIVATE-TOKEN` header, whereas
    `AuthenticatedHTTPStore` uses `Authorization: Bearer`. Both take the *same*
    key from [`get_cds_api_key`](@ref); only the header differs.

Julia alternative: `CDSAPI.jl` implements the same exchange. It is deliberately not
used here — it discovers credentials from `CDSAPI_KEY`/`CDSAPI_URL` rather than this
package's `CDS_API_KEY`, its latest release requires HTTP.jl v2 while this package
pins v1, and it collapses all 4xx responses into one message (losing the
licence-acceptance and constraint-violation detail handled below). If CDS changes
this protocol, adopting CDSAPI.jl is the natural fallback.
"""

using HTTP
import JSON3
import Downloads

const _CDS_RETRIEVE_BASE = "https://cds.climate.copernicus.eu/api/retrieve/v1"

# Where a user accepts a dataset's licence terms (the most common first-run failure).
const _CDS_DATASET_PAGE = "https://cds.climate.copernicus.eu/datasets"

# Where queued/failed jobs can be inspected after login.
const _CDS_REQUESTS_PAGE = "https://cds.climate.copernicus.eu/requests"

# Job statuses that mean "not finished yet"; anything else terminates the poll loop.
const _CDS_PENDING_STATUSES = ("accepted", "running", "queued", "started")

"""
    _cds_headers(token) -> Vector{Pair{String,String}}

Request headers for the CDS Retrieve API.

Note the `PRIVATE-TOKEN` header: the Retrieve API does **not** accept the
`Authorization: Bearer` scheme used by [`AuthenticatedHTTPStore`] for the ARCO Zarr
stores, even though the key itself is the same.
"""
_cds_headers(token::AbstractString) = [
    "PRIVATE-TOKEN" => String(token),
    "Content-Type" => "application/json",
    "Accept" => "application/json",
]

"""
    _cds_post(url, inputs; token) -> HTTP.Response

POST `{"inputs": inputs}` to `url`, without throwing on HTTP error status (callers
inspect the status themselves so they can produce actionable messages).
"""
_cds_post(url::AbstractString, inputs::AbstractDict; token::AbstractString) =
    HTTP.request("POST", url, _cds_headers(token),
        JSON3.write(Dict("inputs" => inputs)); status_exception=false)

"""
    _cds_error_detail(body) -> String

Extract the human-readable explanation from a CDS JSON error body, falling back to
the raw body when it is not JSON or carries no recognised field.
"""
function _cds_error_detail(body::AbstractString)
    isempty(strip(body)) && return "(empty response body)"
    try
        parsed = JSON3.read(body)
        parts = String[]
        for k in (:title, :detail, :traceback, :messages)
            if haskey(parsed, k)
                v = parsed[k]
                v === nothing && continue
                s = v isa AbstractString ? String(v) : string(v)
                isempty(strip(s)) || push!(parts, s)
            end
        end
        isempty(parts) || return join(unique(parts), " — ")
    catch
        # Not JSON; fall through to the raw body.
    end
    return String(first(body, 2000))
end

"""
    _cds_check_response(response, dataset, context; inputs=nothing)

Translate a non-2xx CDS response into an informative error. `context` names the
operation ("submit", "costing", …) and `inputs`, when given, is echoed so a rejected
request can be inspected.

Handles the three failures that dominate in practice:
- **403/401** — the dataset's licence has not been accepted for this account.
- **400** — the request violates the dataset's constraints or size limit; the CDS
  explanation is echoed verbatim because it names the offending field.
- **404** — unknown dataset id.
"""
function _cds_check_response(response, dataset::AbstractString, context::AbstractString;
                             inputs=nothing)
    200 <= response.status < 300 && return nothing

    body = String(response.body)
    detail = _cds_error_detail(body)
    request_note = isnothing(inputs) ? "" : "\nRequest sent:\n$(JSON3.write(inputs))"

    if response.status in (401, 403)
        error("""
        CDS returned HTTP $(response.status) ($(context)) for dataset "$(dataset)".

        $(detail)

        This usually means the dataset's licence has not yet been accepted for your
        account, or the API key is invalid. Accept the terms at:
            $(_CDS_DATASET_PAGE)/$(dataset)?tab=download
        then retry. Verify your key with `get_cds_api_key()`.$(request_note)
        """)
    elseif response.status == 404
        throw(ArgumentError("""
        CDS does not recognise the dataset "$(dataset)" (HTTP 404, $(context)).

        $(detail)
        """))
    elseif 400 <= response.status < 500
        throw(ArgumentError("""
        CDS rejected the request for "$(dataset)" (HTTP $(response.status), $(context)).

        $(detail)$(request_note)
        """))
    else
        error("""
        CDS server error for "$(dataset)" (HTTP $(response.status), $(context)).

        $(detail)$(request_note)
        """)
    end
end

"""
    cds_estimate_cost(dataset, inputs; token) -> (cost, limit)

Ask CDS how large a request is, **without** submitting a job. Returns the request's
`cost` and the server-enforced `limit`; a request with `cost > limit` is rejected at
submission.

For `satellite-albedo` the cost is the number of (variable × timestep) combinations
and the limit is 20 — notably, restricting `area` reduces the *volume* of the result
but not its cost. Callers therefore have to split long time series into several jobs.

This endpoint is cheap and unqueued, so it is worth calling before every submission.

# Example
```julia
cost, limit = cds_estimate_cost("satellite-albedo", request; token)
cost > limit && error("split this request")
```
"""
function cds_estimate_cost(dataset::AbstractString, inputs::AbstractDict;
                           token::AbstractString)
    url = "$(_CDS_RETRIEVE_BASE)/processes/$(dataset)/costing?request_origin=api"
    response = _cds_post(url, inputs; token=token)
    _cds_check_response(response, dataset, "costing"; inputs=inputs)
    body = JSON3.read(String(response.body))
    cost = haskey(body, :cost) ? Float64(body[:cost]) : NaN
    limit = haskey(body, :limit) ? Float64(body[:limit]) : Inf
    return (cost, limit)
end

"""
    cds_valid_options(dataset, inputs; token) -> Dict{String,Any}

Ask CDS which option values remain valid given a partial selection — the live version
of the dataset's constraint table. Useful for validating against the *current* catalogue
rather than a hardcoded copy that can go stale as new years are published.

# Example
```julia
opts = cds_valid_options("satellite-albedo", Dict("sensor" => "olci_and_slstr"); token)
opts["year"]         # ["2018", …, "2024"]
opts["nominal_day"]  # ["10", "20", "28", "29", "30", "31"]
```
"""
function cds_valid_options(dataset::AbstractString, inputs::AbstractDict;
                           token::AbstractString)
    url = "$(_CDS_RETRIEVE_BASE)/processes/$(dataset)/constraints"
    response = _cds_post(url, inputs; token=token)
    _cds_check_response(response, dataset, "constraints"; inputs=inputs)
    return copy(JSON3.read(String(response.body), Dict{String,Any}))
end

"""
    _cds_submit(dataset, inputs; token) -> String

Submit a retrieval job and return its job id. The id is read from the response body's
`jobID`, falling back to the trailing path segment of the `Location` header.
"""
function _cds_submit(dataset::AbstractString, inputs::AbstractDict;
                     token::AbstractString)
    url = "$(_CDS_RETRIEVE_BASE)/processes/$(dataset)/execution"
    response = _cds_post(url, inputs; token=token)
    _cds_check_response(response, dataset, "submit"; inputs=inputs)

    body = JSON3.read(String(response.body))
    if haskey(body, :jobID)
        return String(body[:jobID])
    end
    # Fall back to the Location header, whose last segment is the job id.
    for (name, value) in response.headers
        if lowercase(name) == "location"
            return String(last(split(rstrip(value, '/'), '/')))
        end
    end
    error("""
    CDS accepted the request for "$(dataset)" but returned no job id.
    Response body:
    $(String(first(String(response.body), 1000)))
    """)
end

"""
    _cds_wait(dataset, job_id; token, poll_interval, timeout, verbose) -> Nothing

Poll a job until it succeeds, throwing if it fails or `timeout` seconds elapse.

The poll interval grows geometrically (capped at 60 s) so a long job does not
generate needless traffic. Progress is reported when the status changes and roughly
once a minute otherwise, since jobs routinely take minutes.
"""
function _cds_wait(dataset::AbstractString, job_id::AbstractString;
                   token::AbstractString, poll_interval::Real, timeout::Real,
                   verbose::Bool)
    url = "$(_CDS_RETRIEVE_BASE)/jobs/$(job_id)"
    headers = _cds_headers(token)
    started = time()
    last_status = ""
    last_report = started
    attempt = 0

    while true
        response = HTTP.request("GET", url, headers; status_exception=false)
        _cds_check_response(response, dataset, "job status")
        body = JSON3.read(String(response.body))
        status = haskey(body, :status) ? String(body[:status]) : "unknown"

        if status == "successful"
            verbose && @info "CDS job complete" dataset job_id elapsed_s=round(time() - started; digits=1)
            return nothing
        elseif status in ("failed", "dismissed")
            error("""
            CDS job $(job_id) for "$(dataset)" $(status).

            $(_cds_error_detail(String(response.body)))

            Inspect the request (after login) at $(_CDS_REQUESTS_PAGE)
            """)
        elseif !(status in _CDS_PENDING_STATUSES)
            # Unknown status: keep waiting rather than failing, but say so once.
            status != last_status && verbose &&
                @warn "Unrecognised CDS job status; continuing to poll" dataset job_id status
        end

        elapsed = time() - started
        if elapsed > timeout
            error("""
            Timed out after $(round(elapsed; digits=1))s waiting for CDS job $(job_id)
            ("$(dataset)", last status "$(status)").

            The job may still be running server-side — it is not lost. Check
                $(_CDS_REQUESTS_PAGE)
            or retry with a larger `timeout`.
            """)
        end

        if status != last_status
            verbose && @info "CDS job $(status)" dataset job_id
            last_status = status
            last_report = time()
        elseif verbose && time() - last_report >= 60
            @info "Still waiting on CDS job" dataset job_id status elapsed_s=round(elapsed; digits=1)
            last_report = time()
        end

        attempt += 1
        sleep(min(poll_interval * 1.5^(attempt - 1), 60.0))
    end
end

"""
    _cds_result_href(dataset, job_id; token) -> String

URL of a successful job's downloadable asset.
"""
function _cds_result_href(dataset::AbstractString, job_id::AbstractString;
                          token::AbstractString)
    url = "$(_CDS_RETRIEVE_BASE)/jobs/$(job_id)/results"
    response = HTTP.request("GET", url, _cds_headers(token); status_exception=false)
    _cds_check_response(response, dataset, "job results")
    body = JSON3.read(String(response.body))
    try
        return String(body[:asset][:value][:href])
    catch
        error("""
        CDS job $(job_id) succeeded but its results carry no asset href.
        Response body:
        $(String(first(String(response.body), 1000)))
        """)
    end
end

"""
    cds_retrieve(dataset, inputs, dest; token, poll_interval=10, timeout=1800,
                 verbose=true, check_cost=true) -> String

Order `dataset` from the CDS Retrieve API with request `inputs` and save the result to
`dest`, returning `dest`.

Blocks while the job is queued and processed server-side — typically tens of seconds
to several minutes. The downloaded file is usually a ZIP of NetCDFs, but CDS may
return a bare NetCDF for a single-field request, so callers should sniff the content
rather than assume.

# Keywords
- `token`: CDS API key (see [`get_cds_api_key`](@ref)).
- `poll_interval`: initial seconds between status polls; grows geometrically to 60 s.
- `timeout`: seconds to wait before giving up. The job is not cancelled — the error
  reports its id so it can be inspected or collected later.
- `check_cost`: when `true`, call [`cds_estimate_cost`](@ref) first and fail locally
  if the request exceeds the server limit, which gives a clearer message than the
  server's rejection.
"""
function cds_retrieve(dataset::AbstractString, inputs::AbstractDict, dest::AbstractString;
                      token::AbstractString, poll_interval::Real=10, timeout::Real=1800,
                      verbose::Bool=true, check_cost::Bool=true)
    if check_cost
        cost, limit = cds_estimate_cost(dataset, inputs; token=token)
        if isfinite(cost) && isfinite(limit) && cost > limit
            throw(ArgumentError("""
            This CDS request is too large: cost $(cost) exceeds the limit of $(limit)
            for "$(dataset)". Split it into smaller requests (for `satellite-albedo`
            the cost is the number of variable × timestep combinations, and narrowing
            `area` does not reduce it).

            Request:
            $(JSON3.write(inputs))
            """))
        end
    end

    job_id = _cds_submit(dataset, inputs; token=token)
    verbose && @info "Submitted CDS request" dataset job_id
    _cds_wait(dataset, job_id; token=token, poll_interval=poll_interval,
              timeout=timeout, verbose=verbose)
    href = _cds_result_href(dataset, job_id; token=token)

    mkpath(dirname(dest))
    tmp = dest * ".part"
    try
        Downloads.download(href, tmp)
        mv(tmp, dest; force=true)
    catch e
        isfile(tmp) && rm(tmp; force=true)
        error("Failed to download CDS job $(job_id) result for \"$(dataset)\" from $(href):\n$e")
    end
    return dest
end
