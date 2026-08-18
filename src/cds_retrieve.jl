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
    CDSQueueLimitError(dataset, job_id, detail)

A job CDS answered with `status: "rejected"` because too many requests for that dataset
were already queued for this account.

This is a distinct failure from a bad request: the request itself is legal and succeeds
verbatim once the queue drains, so [`cds_retrieve`](@ref) catches this type and resubmits
rather than propagating it. It surfaces when several jobs are ordered concurrently — the
limit is per dataset, not per job — and the rejection is only visible in the *job's*
status, not in the (HTTP 200) submission response.
"""
struct CDSQueueLimitError <: Exception
    dataset::String
    job_id::String
    detail::String
end

function Base.showerror(io::IO, e::CDSQueueLimitError)
    print(io, """
    CDS rejected job $(e.job_id) for "$(e.dataset)": too many requests already queued
    for this dataset.

    $(e.detail)

    Reduce the number of concurrent jobs, or wait for the queued ones to finish.
    """)
end

# Fragment of the CDS rejection message that identifies the per-dataset queue limit, as
# opposed to any other cause of a "rejected" status.
const _CDS_QUEUE_LIMIT_HINT = "queued requests for this dataset is temporarily limited"

"""
    CDSTransientError(dataset, context, status, detail)

A CDS response that says nothing about the request: a 5xx from the service or its
fronting proxy.

Distinct from the 4xx errors, which are the caller's fault and must not be retried. A
502/503/504 typically means nginx could not reach the application for a moment; the
identical call succeeds seconds later. [`_cds_retry`](@ref) catches this and retries.

This exists because the failure it names is unavoidable over long runs and was fatal
before: a multi-hour global albedo order died on a single HTTP 502 from `/costing` — the
*cheap, unqueued* cost pre-check — discarding two batches of folded work and ~11 h of
downloads.
"""
struct CDSTransientError <: Exception
    dataset::String
    context::String
    status::Int
    detail::String
    retry_after::Float64   # server-requested wait, seconds; 0.0 when unspecified
end

CDSTransientError(dataset, context, status, detail) =
    CDSTransientError(dataset, context, status, detail, 0.0)

function Base.showerror(io::IO, e::CDSTransientError)
    print(io, """
    CDS server error for "$(e.dataset)" (HTTP $(e.status), $(e.context)).

    $(e.detail)

    This is a service-side fault, not a problem with the request; it is normally
    transient and worth retrying.
    """)
end

# Backoff bounds and attempt cap for retrying a transient (5xx / network) failure. Short
# compared with the queue-limit backoff because the blocker is a service hiccup, not other
# jobs of ours draining.
const _CDS_TRANSIENT_RETRY_BASE = 5.0
const _CDS_TRANSIENT_RETRY_MAX = 120.0
# Generous on purpose. Eight attempts is only ~6.6 min of backoff, which a CDS maintenance
# window outlasts easily — and the cost of giving up is hours of queue time plus, for the
# albedo folder, the batch of downloads behind it. ECMWF's `cdsapi` sets `retry_max=500` at
# `sleep_max=120` (~16 h), i.e. effectively "never give up"; `deadline` is what actually
# bounds us, so a high count here just stops the *attempt cap* from being the binding limit.
const _CDS_TRANSIENT_ATTEMPTS = 60

# Non-5xx statuses that are nonetheless worth retrying. 408 and 429 are 4xx by number but
# say nothing bad about the request: the service is asking us to come back. Matches the set
# ECMWF's own `cdsapi` retries (`Client.robust`: 500, 502, 503, 504, 408, 429).
const _CDS_TRANSIENT_STATUSES = (408, 429)

# Cap on an honoured `Retry-After`, so a pathological header cannot park a run for hours.
const _CDS_RETRY_AFTER_MAX = 300.0

"""
    _cds_retry_after(response) -> Float64

Seconds the server asked us to wait, from the `Retry-After` header, or `0.0` when absent
or unparseable. Only the delta-seconds form is honoured (the HTTP-date form is legal but
CDS does not use it); the value is clamped to `_CDS_RETRY_AFTER_MAX`.
"""
function _cds_retry_after(response)
    raw = HTTP.header(response, "Retry-After", "")
    isempty(strip(raw)) && return 0.0
    seconds = tryparse(Float64, strip(raw))
    isnothing(seconds) && return 0.0
    return clamp(seconds, 0.0, _CDS_RETRY_AFTER_MAX)
end

"""
    _cds_is_transient(err) -> Bool

Whether `err` is worth retrying: a [`CDSTransientError`](@ref) (a 5xx, or a 408/429 —
see `_CDS_TRANSIENT_STATUSES`) or a connection-level failure (DNS, connect, timeout,
reset, truncated body).

Deliberately narrow. `ArgumentError` (4xx: bad request, unknown dataset) and plain
`ErrorException` (a job CDS reports as `failed`/`dismissed`, or a request that timed out)
are *not* transient — retrying those is what produced the endless resubmit loop this
package documents at `satellite_albedo`.
"""
_cds_is_transient(err) =
    err isa CDSTransientError ||
    err isa HTTP.Exceptions.ConnectError ||
    err isa HTTP.Exceptions.TimeoutError ||
    err isa HTTP.Exceptions.RequestError ||
    err isa Downloads.RequestError ||
    err isa Base.IOError ||
    err isa EOFError

"""
    _cds_retry(f; dataset, context, attempts, verbose, deadline) -> f()

Call `f`, retrying while it fails transiently (see [`_cds_is_transient`](@ref)) with
exponential backoff from `_CDS_TRANSIENT_RETRY_BASE` to `_CDS_TRANSIENT_RETRY_MAX`.

Every CDS interaction goes through this — costing, submission, status polls, results and
the asset download — because any of them can catch a proxy hiccup, and a single one
killing a multi-hour run is not acceptable. Non-transient errors propagate on the first
attempt, unchanged.

`deadline` (an absolute `time()` value) bounds the retrying so a caller's own timeout is
still respected; the final failure is rethrown rather than wrapped, so the message the
user sees is the service's own.

The mechanism itself lives in [`_retry_transient`](@ref) (`utils.jl`) and is shared with
the NASA Earthdata client; this wrapper supplies only what is CDS-specific — the error
taxonomy, the `_CDS_*` constants, the `Retry-After` accessor, and the log record.

!!! note "Resubmission may duplicate a job"
    A 5xx on *submission* can in principle mean the job was created and only the
    response was lost, in which case the retry queues a second identical job. That is
    tolerable — a duplicate is dismissible and costs one queue slot — and far cheaper
    than failing a run that has hours of folded state behind it.
"""
_cds_retry(f; dataset::AbstractString, context::AbstractString,
           attempts::Integer=_CDS_TRANSIENT_ATTEMPTS, verbose::Bool=true,
           deadline::Float64=Inf) =
    _retry_transient(f; is_transient=_cds_is_transient, context=context,
        attempts=attempts, base=_CDS_TRANSIENT_RETRY_BASE,
        max_backoff=_CDS_TRANSIENT_RETRY_MAX, deadline=deadline, verbose=verbose,
        # A server-supplied `Retry-After` (429) knows better than our curve does.
        retry_after = err -> err isa CDSTransientError ? err.retry_after : 0.0,
        on_retry = (attempt, backoff, err) ->
            @warn "Transient CDS failure; retrying" dataset context attempt backoff_s = round(backoff; digits=1) error = sprint(showerror, err))

# Backoff bounds, in seconds, for resubmitting a queue-limited job. Generous because the
# blocker is other jobs of ours finishing, which takes minutes, not milliseconds.
const _CDS_QUEUE_RETRY_BASE = 30.0
const _CDS_QUEUE_RETRY_MAX = 300.0

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
    elseif response.status in _CDS_TRANSIENT_STATUSES
        # 408 Request Timeout and 429 Too Many Requests are 4xx by number but transient by
        # nature: the request is legal and the service is asking us to come back. ECMWF's
        # own `cdsapi` retries exactly these alongside the 5xx set (`Client.robust`), and a
        # 429 reaching the user as "CDS rejected the request" would be actively misleading.
        # 429 usually carries `Retry-After`; honour it instead of our own backoff.
        throw(CDSTransientError(String(dataset), String(context), Int(response.status),
                                detail * request_note, _cds_retry_after(response)))
    elseif 400 <= response.status < 500
        throw(ArgumentError("""
        CDS rejected the request for "$(dataset)" (HTTP $(response.status), $(context)).

        $(detail)$(request_note)
        """))
    else
        # 5xx says nothing about the request — it is the service or its proxy faltering, and
        # the identical call usually succeeds moments later. Typed so `_cds_retry` can pick
        # it out; a bare `error()` here once killed an 11 h run on one nginx 502.
        throw(CDSTransientError(String(dataset), String(context), Int(response.status),
                                detail * request_note, _cds_retry_after(response)))
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
                           token::AbstractString, verbose::Bool=true,
                           deadline::Float64=Inf)
    url = "$(_CDS_RETRIEVE_BASE)/processes/$(dataset)/costing?request_origin=api"
    response = _cds_retry(; dataset=dataset, context="costing", verbose=verbose,
                          deadline=deadline) do
        r = _cds_post(url, inputs; token=token)
        _cds_check_response(r, dataset, "costing"; inputs=inputs)
        r
    end
    body = JSON3.read(String(response.body))
    cost = haskey(body, :cost) ? Float64(body[:cost]) : NaN
    limit = haskey(body, :limit) ? Float64(body[:limit]) : Inf
    return (cost, limit)
end

"""
    _cds_advisory_cost(dataset, inputs; token, verbose, deadline) -> (cost, limit)

[`cds_estimate_cost`](@ref), downgraded from a gate to a courtesy: returns `(NaN, Inf)`
instead of throwing when the costing endpoint cannot be reached, so the caller submits the
request and lets the *server* judge it.

This distinction is load-bearing. The pre-check exists only to turn an eventual server
rejection into a clearer local error, and ECMWF's own `cdsapi` does not call the endpoint
at all — so a `/costing` outage must never be fatal. It was: one nginx **HTTP 502 from
this cheap, unqueued call** killed an 11 h global albedo run that had folded 8 timesteps
and downloaded 86 GB. Making the call retry ([`_cds_retry`](@ref)) narrows the window but
does not close it — exhausting the attempts would still abort a run over a check whose
answer is not needed.

A **non**-transient error still propagates: a 4xx means the request is malformed or the
licence unaccepted, and submitting it would fail the same way with a worse message.

`estimate` exists so the tests can exercise both branches of that policy without a network
round trip; callers should leave it at its default.
"""
function _cds_advisory_cost(dataset::AbstractString, inputs::AbstractDict;
                            token::AbstractString, verbose::Bool=true,
                            deadline::Float64=Inf, estimate=cds_estimate_cost)
    try
        return estimate(dataset, inputs; token=token, verbose=verbose,
                        deadline=deadline)
    catch err
        _cds_is_transient(err) || rethrow()
        verbose && @warn "CDS cost pre-check unavailable; submitting the request anyway" dataset error = sprint(showerror, err)
        return (NaN, Inf)
    end
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
                     token::AbstractString, verbose::Bool=true, deadline::Float64=Inf)
    url = "$(_CDS_RETRIEVE_BASE)/processes/$(dataset)/execution"
    response = _cds_retry(; dataset=dataset, context="submit", verbose=verbose,
                          deadline=deadline) do
        r = _cds_post(url, inputs; token=token)
        _cds_check_response(r, dataset, "submit"; inputs=inputs)
        r
    end

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
        # A poll that hits a 5xx or a dropped connection must not lose the job: it is
        # running server-side regardless, so retry the *poll* and keep waiting. Bounded by
        # this wait's own deadline so a permanently broken service still gives up.
        response = _cds_retry(; dataset=dataset, context="job status", verbose=verbose,
                              deadline=started + timeout) do
            r = HTTP.request("GET", url, headers; status_exception=false)
            _cds_check_response(r, dataset, "job status")
            r
        end
        body = JSON3.read(String(response.body))
        status = haskey(body, :status) ? String(body[:status]) : "unknown"

        if status == "successful"
            verbose && @info "CDS job complete" dataset job_id elapsed_s=round(time() - started; digits=1)
            return nothing
        elseif status == "rejected"
            # A rejection carries its reason in the *results* endpoint, not in the job
            # body, so fetch it before deciding whether this is retryable.
            detail = _cds_rejection_detail(dataset, job_id, String(response.body);
                                           token=token)
            if occursin(_CDS_QUEUE_LIMIT_HINT, detail)
                throw(CDSQueueLimitError(String(dataset), String(job_id), detail))
            end
            error("""
            CDS job $(job_id) for "$(dataset)" was rejected.

            $(detail)

            Inspect the request (after login) at $(_CDS_REQUESTS_PAGE)
            """)
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
    _cds_rejection_detail(dataset, job_id, job_body; token) -> String

Why a job was rejected. The job document itself reports only `status: "rejected"`, with
no reason; the explanation (a traceback naming e.g. the per-dataset queue limit) lives on
the `/results` endpoint, which answers 4xx for a rejected job. So query it, tolerating any
failure — the caller is already on an error path and a missing reason must not mask it.
"""
function _cds_rejection_detail(dataset::AbstractString, job_id::AbstractString,
                               job_body::AbstractString; token::AbstractString)
    try
        url = "$(_CDS_RETRIEVE_BASE)/jobs/$(job_id)/results"
        response = HTTP.request("GET", url, _cds_headers(token); status_exception=false)
        detail = _cds_error_detail(String(response.body))
        isempty(strip(detail)) || return detail
    catch err
        @debug "Could not fetch CDS rejection reason" dataset job_id exception = err
    end
    return _cds_error_detail(job_body)
end

"""
    _cds_result_href(dataset, job_id; token) -> String

URL of a successful job's downloadable asset.
"""
function _cds_result_href(dataset::AbstractString, job_id::AbstractString;
                          token::AbstractString, verbose::Bool=true,
                          deadline::Float64=Inf)
    url = "$(_CDS_RETRIEVE_BASE)/jobs/$(job_id)/results"
    response = _cds_retry(; dataset=dataset, context="job results", verbose=verbose,
                          deadline=deadline) do
        r = HTTP.request("GET", url, _cds_headers(token); status_exception=false)
        _cds_check_response(r, dataset, "job results")
        r
    end
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

Safe to call from several concurrent tasks: a job CDS rejects for exceeding the
per-dataset queued-request limit is resubmitted with backoff (see
[`CDSQueueLimitError`](@ref)) instead of failing.

Every step — costing, submission, status polls, results, download — is also retried on a
**transient** failure: a CDS 5xx (see [`CDSTransientError`](@ref)) or a connection-level
error. Bad requests (4xx) and jobs the server reports as `failed` still fail immediately.
This matters for long runs, where a momentary proxy fault is otherwise indistinguishable
from a real error and takes hours of accumulated work with it.

# Keywords
- `token`: CDS API key (see [`get_cds_api_key`](@ref)).
- `poll_interval`: initial seconds between status polls; grows geometrically to 60 s.
- `timeout`: seconds to wait before giving up, covering submission retries as well as the
  wait itself. The job is not cancelled — the error reports its id so it can be inspected
  or collected later.
- `check_cost`: when `true`, call [`cds_estimate_cost`](@ref) first and fail locally
  if the request exceeds the server limit, which gives a clearer message than the
  server's rejection. Advisory only — if the costing endpoint is unreachable the
  request is submitted anyway (see the comment at the call site).
"""
function cds_retrieve(dataset::AbstractString, inputs::AbstractDict, dest::AbstractString;
                      token::AbstractString, poll_interval::Real=10, timeout::Real=1800,
                      verbose::Bool=true, check_cost::Bool=true)
    started = time()
    deadline = timeout == Inf ? Inf : started + Float64(timeout)

    if check_cost
        cost, limit = _cds_advisory_cost(dataset, inputs; token=token, verbose=verbose,
                                         deadline=deadline)
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

    # Submit, and resubmit if CDS rejects the job for having too many of this dataset's
    # requests already queued. That rejection is about account state, not the request, so
    # the identical submission succeeds once the queue drains — retrying here is what lets
    # callers order several jobs concurrently without hand-tuning their fan-out. `timeout`
    # bounds the retry loop as well as each wait, so a saturated queue cannot spin forever.
    job_id = ""
    attempt = 0
    while true
        attempt += 1
        job_id = _cds_submit(dataset, inputs; token=token, verbose=verbose,
                             deadline=deadline)
        verbose && @info "Submitted CDS request" dataset job_id
        try
            _cds_wait(dataset, job_id; token=token, poll_interval=poll_interval,
                      timeout=max(timeout - (time() - started), 0), verbose=verbose)
            break
        catch err
            err isa CDSQueueLimitError || rethrow()
            elapsed = time() - started
            backoff = min(_CDS_QUEUE_RETRY_BASE * 2.0^(attempt - 1), _CDS_QUEUE_RETRY_MAX)
            if elapsed + backoff >= timeout
                error("""
                CDS kept rejecting this "$(dataset)" request for exceeding the per-dataset
                queued-request limit, and $(round(timeout; digits=1))s elapsed without a
                free slot (attempt $(attempt)).

                Lower the number of concurrent jobs, or raise `timeout`.
                """)
            end
            verbose && @info "CDS queue full; resubmitting after backoff" dataset attempt backoff_s = round(backoff; digits=1)
            sleep(backoff)
        end
    end
    href = _cds_result_href(dataset, job_id; token=token, verbose=verbose,
                            deadline=deadline)

    mkpath(dirname(dest))
    tmp = dest * ".part"
    # The download is the longest single transfer in the exchange (tens of GB for a global
    # albedo timestep), so it is the most likely to be cut off mid-stream. Retry it: the
    # partial file is discarded first, so each attempt starts clean.
    try
        _cds_retry(; dataset=dataset, context="download", verbose=verbose,
                   deadline=deadline) do
            isfile(tmp) && rm(tmp; force=true)
            Downloads.download(href, tmp)
        end
        mv(tmp, dest; force=true)
    catch e
        isfile(tmp) && rm(tmp; force=true)
        error("Failed to download CDS job $(job_id) result for \"$(dataset)\" from $(href):\n$e")
    end
    return dest
end
