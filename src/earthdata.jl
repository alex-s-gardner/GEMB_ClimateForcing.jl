"""
Client for NASA **Earthdata** — granule discovery via CMR and authenticated download
from a DAAC.

This is a third access protocol, distinct from the two CDS ones already in this package:

| service | discovery | auth header | unit of transfer |
|---|---|---|---|
| ECMWF ARCO Zarr (`authenticated_http_store.jl`) | chunk keys | `Authorization: Bearer` (CDS key) | one chunk, byte range |
| CDS Retrieve (`cds_retrieve.jl`) | job submission | `PRIVATE-TOKEN` (CDS key) | an ordered ZIP |
| NASA Earthdata (here) | CMR, via `EarthData.jl` | `Authorization: Bearer` (**EDL** token) | a whole granule file |

There is no job queue: a CMR query answers immediately and a granule is fetched straight
off HTTPS. The cost model is therefore bandwidth, not latency — the opposite of the CDS
albedo path, where waiting is everything and the transfer is incidental.

!!! warning "Two different bearer tokens"
    Both this and the ARCO Zarr store send `Authorization: Bearer`, but the *tokens are
    unrelated*: that one is a CDS key from [`get_cds_api_key`](@ref), this one an Earthdata
    Login (EDL) token from [`get_earthdata_token`](@ref). Passing one where the other is
    expected fails with an unhelpful 401.

## Relationship to `EarthData.jl`

**Granule discovery is delegated to [`EarthData.jl`](https://github.com/evetion/EarthData.jl)**
(JuliaGeo), not reimplemented — the HTTP query, `CMR-Search-After` pagination, the UMM-G
schema and CMR error parsing are all upstream's. `_cmr_granules` below is a shim that flattens
the typed record into the `(; id, url, bytes)` triple this package consumes.

Upstream's typed UMM-G schema is *better* than parsing the legacy `granules.json` by hand,
which is why the shim reads it rather than the older endpoint:

  - `RelatedUrls[].Type == "GET DATA"` selects the download URL exactly. MCD43A3 entries
    carry **eight** related URLs — including an `s3credentials` endpoint and browse
    imagery — so any looser filter picks the wrong one.
  - `DataGranule.ArchiveAndDistributionInformation[].{Size, SizeUnit}` reports the size
    **with its unit** (`(0.537697, "MB")`), instead of the legacy endpoint's
    unit-less-but-actually-megabytes `granule_size`.

What is still local, and why — each is an upstream gap being contributed back rather than a
deliberate fork:

| local | upstream status |
|---|---|
| `_cmr_granules` calls `EarthData.request`, not `granules` | PR 1: spatial query params (`bounding_box` &c.) |
| `_granule_data_url`, `_granule_bytes`, `_size_unit_factor` | PR 2: `Type`-based selection + unit-aware size |
| `get_earthdata_token`, `earthdata_token_from_netrc` | PR 3: EDL bearer auth (upstream is netrc-only) |
| `EarthdataTransientError`, `_earthdata_is_transient`, `_classifying_requester`, `_earthdata_download` | PR 4: retry policy + size-verified download |
| `_cmr_granules`' `page_num` loop | PR 5: **bug** — `all=true` sends `page_num` *with* `CMR-Search-After`, which CMR rejects |

Delete each block when its PR is released; none of them should grow new features here.

!!! warning "Pinned to `main`, and that is not optional"
    The registered `EarthData` v0.1.0 is the 2023 tree and **cannot parse a present-day CMR
    response at all** (a `MethodError` from the UMM-G schema on any real query). `Project.toml`
    therefore pins `main` via `[sources]`; see the comment there before touching it.
"""

import Base64      # Basic auth for the (optional, explicit) EDL token endpoints
import EarthData   # CMR granule discovery — see "Relationship to EarthData.jl" above

# ---------------------------------------------------------------------------- endpoints

# Where a user mints an EDL token, and where the app-approval / EULA state lives. Both
# appear in error messages, because a 401 and a 403 from a DAAC are otherwise silent about
# which of the two is wrong.
const _EDL_TOKEN_PAGE = "https://urs.earthdata.nasa.gov/users/tokens"
const _EDL_TOKEN_API = "https://urs.earthdata.nasa.gov/api/users"

# CMR's page size ceiling. Requesting more is an error, not a silent clamp.
const _CMR_PAGE_SIZE = 2000

# ------------------------------------------------------------------------ credentials

"""
    get_earthdata_token() -> String

Get the NASA Earthdata Login (EDL) bearer token, from the environment or a token file.

Searches, in order:
1. `ENV["EARTHDATA_TOKEN"]`
2. `~/.edl_token` (first non-empty, non-comment line)

# Throws
- `ErrorException` with setup instructions if no token is found.

!!! note "Deliberately offline"
    This never contacts the network, matching [`get_cds_api_key`](@ref). Minting a token
    is a *separate*, explicit call — see [`earthdata_token_from_netrc`](@ref) — because an
    account may hold only **two** tokens at once, so creation must never be reachable from
    a credential lookup inside a retry loop.

`EARTHDATA_TOKEN` is a widely used community convention rather than a NASA-defined
variable; nothing in the Earthdata tooling reads it for you.

# Example
```julia
ENV["EARTHDATA_TOKEN"] = "eyJ0eXAi..."
alb = compute_glacier_ice_albedo_modis(lat, lon, 2019)
```
"""
function get_earthdata_token()
    tok = get(ENV, "EARTHDATA_TOKEN", "")
    isempty(strip(tok)) || return String(strip(tok))

    path = joinpath(homedir(), ".edl_token")
    if isfile(path)
        for line in readlines(path)
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            return String(s)
        end
    end

    error("""
    No NASA Earthdata Login token found. Get one from
        $(_EDL_TOKEN_PAGE)
    and either

        export EARTHDATA_TOKEN="your-token-here"

    or write it to ~/.edl_token. A token lasts 60 days, and an account may hold at most
    two at a time — reuse an existing one rather than minting a third, which fails.

    Alternatively, with ~/.netrc credentials for urs.earthdata.nasa.gov in place:
        token = earthdata_token_from_netrc()
    """)
end

"""
    earthdata_token_from_netrc(; create=false, machine="urs.earthdata.nasa.gov") -> String

Fetch an existing EDL bearer token using the Earthdata Login username/password stored in
`~/.netrc`, minting one only if `create=true`.

This is a **network call**, deliberately kept out of [`get_earthdata_token`](@ref): a
credential lookup that can hit the network is a lookup that can hang or rate-limit inside
a retry loop.

!!! warning "Two concurrent tokens, maximum"
    EDL allows an account **two** live tokens; requesting a third returns HTTP 403. With
    `create=false` (the default) this only *lists* tokens and returns the first, so it
    cannot exhaust the quota. `create=true` uses `find_or_create_token`, which is
    idempotent when a token already exists but will consume the second slot when none
    does. Never call this from a retry loop.

Tokens are printed nowhere and returned as a plain `String`; treat it as a secret with the
same care as a password (it authenticates as the account for 60 days).
"""
function earthdata_token_from_netrc(; create::Bool=false,
                                    machine::AbstractString="urs.earthdata.nasa.gov")
    user, pass = _netrc_credentials(machine)
    auth = "Basic " * Base64.base64encode(string(user, ":", pass))
    headers = ["Authorization" => auth, "Accept" => "application/json"]

    endpoint = create ? "$(_EDL_TOKEN_API)/token" : "$(_EDL_TOKEN_API)/tokens"
    method = create ? "POST" : "GET"
    r = HTTP.request(method, endpoint, headers; status_exception=false)
    if r.status >= 400
        error("""
        Earthdata Login rejected the token request ($(method) $(endpoint), HTTP $(r.status)):

        $(String(r.body))

        HTTP 401 means the ~/.netrc credentials for "$(machine)" are wrong. HTTP 403 on a
        creation request usually means the two-token limit is already reached — list them
        at $(_EDL_TOKEN_PAGE) and reuse one.
        """)
    end

    body = JSON3.read(String(r.body))
    # `GET /tokens` answers with an array; `POST /token` with a single object.
    entry = body isa JSON3.Array ? (isempty(body) ? nothing : first(body)) : body
    if isnothing(entry) || !haskey(entry, :access_token)
        error("""
        Earthdata Login returned no usable token. With `create=false` this means the
        account currently holds none; call `earthdata_token_from_netrc(create=true)` or
        mint one at $(_EDL_TOKEN_PAGE).
        """)
    end
    return String(entry.access_token)
end

"""
    _netrc_credentials(machine) -> (login, password)

Read `login`/`password` for `machine` out of `~/.netrc`.

A minimal parser: `.netrc` is a whitespace-separated token stream, so the tokens are
scanned in order and the values following the target `machine` are taken until the next
`machine`/`default` entry. The `macdef` form is not supported (Earthdata does not use it).
"""
function _netrc_credentials(machine::AbstractString)
    path = get(ENV, "NETRC", joinpath(homedir(), ".netrc"))
    isfile(path) || error("""
        No .netrc file at $(path), so Earthdata Login credentials cannot be read.
        Add a stanza:

            machine $(machine) login YOUR_USERNAME password YOUR_PASSWORD

        and `chmod 600` it. Or skip .netrc entirely and set EARTHDATA_TOKEN.
        """)

    tokens = split(read(path, String))
    login = password = nothing
    inside = false
    i = 1
    while i <= length(tokens)
        t = tokens[i]
        if t == "machine"
            inside = i + 1 <= length(tokens) && tokens[i+1] == machine
            i += 2
            continue
        elseif t == "default"
            inside = true
            i += 1
            continue
        elseif inside && t in ("login", "password") && i + 1 <= length(tokens)
            t == "login" ? (login = tokens[i+1]) : (password = tokens[i+1])
            i += 2
            !isnothing(login) && !isnothing(password) && break
            continue
        end
        i += 1
    end

    (isnothing(login) || isnothing(password)) && error("""
        $(path) has no login/password for machine "$(machine)".
        """)
    return String(login), String(password)
end

# --------------------------------------------------------------- transient-error policy

# Backoff bounds and attempt cap, mirroring the CDS constants' reasoning: the cap must not
# be the binding limit — the caller's `deadline` is. A granule download is minutes, so a
# maintenance window is worth waiting out rather than discarding the run.
const _EARTHDATA_RETRY_BASE = 2.0
const _EARTHDATA_RETRY_MAX = 60.0
const _EARTHDATA_RETRY_ATTEMPTS = 60

# Non-5xx statuses worth retrying: legal request, service asking us to come back. Same set
# as `_CDS_TRANSIENT_STATUSES`, for the same reason.
const _EARTHDATA_TRANSIENT_STATUSES = (408, 429, 500, 502, 503, 504)

# Cap on an honoured `Retry-After`, so a pathological header cannot park a run for hours.
const _EARTHDATA_RETRY_AFTER_MAX = 300.0

"""
    EarthdataTransientError(context, status, detail, retry_after)

An Earthdata/CMR response that says nothing about the request: a 5xx, or a 408/429.

Distinct from 401/403/404, which are **permanent** here and must fail fast:

- **401** — the bearer token is missing, malformed or expired (they last 60 days).
- **403** — the token is valid but the account has not accepted the collection's EULA or
  approved the DAAC application. Retrying can never fix it; this is the same lesson as the
  CDS licence-403.
- **404** — the granule genuinely is not there.

Retrying any of those is what produces an endless request loop, so the split is asserted
in both directions by the tests.
"""
struct EarthdataTransientError <: Exception
    context::String
    status::Int
    detail::String
    retry_after::Float64   # server-requested wait, seconds; 0.0 when unspecified
end

EarthdataTransientError(context, status, detail) =
    EarthdataTransientError(context, status, detail, 0.0)

function Base.showerror(io::IO, e::EarthdataTransientError)
    print(io, """
    NASA Earthdata server error ($(e.context), HTTP $(e.status)).

    $(e.detail)

    This is a service-side fault, not a problem with the request; it is normally
    transient and worth retrying.
    """)
end

"""
    _earthdata_retry_after(response) -> Float64

Seconds the server asked us to wait, from the `Retry-After` header, or `0.0` when absent
or unparseable. Only the delta-seconds form is honoured; clamped to
`_EARTHDATA_RETRY_AFTER_MAX`.
"""
function _earthdata_retry_after(response)
    raw = HTTP.header(response, "Retry-After", "")
    isempty(strip(raw)) && return 0.0
    seconds = tryparse(Float64, strip(raw))
    isnothing(seconds) && return 0.0
    return clamp(seconds, 0.0, _EARTHDATA_RETRY_AFTER_MAX)
end

"""
    _earthdata_is_transient(err) -> Bool

Whether `err` is worth retrying: an [`EarthdataTransientError`](@ref) or a
connection-level failure (DNS, connect, timeout, reset, truncated body).

Deliberately narrow, and the mirror of [`_cds_is_transient`](@ref). `ArgumentError`
(401/403/404, a bad CMR query) and plain `ErrorException` (a size mismatch, a missing
download link) are *not* transient.
"""
_earthdata_is_transient(err) =
    err isa EarthdataTransientError ||
    err isa HTTP.Exceptions.ConnectError ||
    err isa HTTP.Exceptions.TimeoutError ||
    err isa HTTP.Exceptions.RequestError ||
    err isa Downloads.RequestError ||
    err isa Base.IOError ||
    err isa EOFError

"""
    _earthdata_retry(f; context, attempts, verbose, deadline) -> f()

Call `f`, retrying while it fails transiently (see [`_earthdata_is_transient`](@ref)) with
exponential backoff from `_EARTHDATA_RETRY_BASE` to `_EARTHDATA_RETRY_MAX`.

The mechanism lives in [`_retry_transient`](@ref) (`utils.jl`) and is shared with the CDS
client; this wrapper supplies only the Earthdata-specific error taxonomy, constants and
log record.
"""
_earthdata_retry(f; context::AbstractString,
                 attempts::Integer=_EARTHDATA_RETRY_ATTEMPTS, verbose::Bool=true,
                 deadline::Float64=Inf) =
    _retry_transient(f; is_transient=_earthdata_is_transient, context=context,
        attempts=attempts, base=_EARTHDATA_RETRY_BASE,
        max_backoff=_EARTHDATA_RETRY_MAX, deadline=deadline, verbose=verbose,
        retry_after = err -> err isa EarthdataTransientError ? err.retry_after : 0.0,
        on_retry = (attempt, backoff, err) ->
            @warn "Transient Earthdata failure; retrying" context attempt backoff_s = round(backoff; digits=1) error = sprint(showerror, err))

"""
    _earthdata_check_response(r, context)

Throw for an error response: [`EarthdataTransientError`](@ref) for a retryable status, an
`ArgumentError` carrying actionable advice for a permanent one. Returns `nothing` for a
success.
"""
function _earthdata_check_response(r, context::AbstractString)
    r.status < 400 && return nothing
    detail = _earthdata_body_excerpt(r)
    if r.status in _EARTHDATA_TRANSIENT_STATUSES || r.status >= 500
        throw(EarthdataTransientError(context, r.status, detail,
                                      _earthdata_retry_after(r)))
    elseif r.status == 401
        throw(ArgumentError("""
        Earthdata Login rejected the token ($(context), HTTP 401).

        $(detail)

        Tokens expire after 60 days. Check EARTHDATA_TOKEN, or list/mint one at
        $(_EDL_TOKEN_PAGE).
        """))
    elseif r.status == 403
        throw(ArgumentError("""
        Earthdata accepted the token but refused the request ($(context), HTTP 403).

        $(detail)

        This is almost always an unaccepted end-user licence or an unapproved
        application, not a bad token — log in at https://urs.earthdata.nasa.gov/ and
        approve the DAAC application, then retry. Retrying without doing so cannot help.
        """))
    else
        throw(ArgumentError("""
        Earthdata request failed ($(context), HTTP $(r.status)).

        $(detail)
        """))
    end
end

# Response bodies can be a full HTML error page; keep the message readable.
function _earthdata_body_excerpt(r; limit::Integer=800)
    body = try
        String(r.body)
    catch
        ""
    end
    s = strip(body)
    isempty(s) && return "(empty response body)"
    return length(s) <= limit ? s : string(first(s, limit), "\n… (truncated)")
end

# ------------------------------------------------------------------- granule discovery

"""
    _cmr_granules(; short_name, version, temporal=nothing, bounding_box=nothing,
                  provider=nothing, page_size=_CMR_PAGE_SIZE, verbose=true,
                  deadline=Inf) -> Vector{NamedTuple}

Search NASA's Common Metadata Repository for granules, returning
`(; id, url, bytes)` per granule: the archive filename, an HTTPS download URL, and the size
in **bytes**.

A thin shim over `EarthData.request` — pagination (`CMR-Search-After`), the UMM-G schema and
error parsing are all upstream's. `temporal` is a CMR range string
(`"2019-07-01T00:00:00Z,2019-07-01T23:59:59Z"`) and `bounding_box` is `"W,S,E,N"`.

!!! warning "`temporal` is not a nominal-date filter"
    A granule's temporal extent spans the product's full retrieval window — 16 days for
    MCD43A3 — so a *one-day* `temporal` query returns every granule whose window overlaps
    that day: 16 distinct nominal dates, not one. Callers must filter on the nominal date
    encoded in the granule id (see `_modis_granule_date`). Verified live: a one-day,
    two-tile query returned `cmr-hits: 32`.

Granules with no HTTPS `GET DATA` URL are dropped, since there is nothing to download.
"""
function _cmr_granules(; short_name::AbstractString, version::AbstractString,
                       temporal::Union{Nothing,AbstractString}=nothing,
                       bounding_box::Union{Nothing,AbstractString}=nothing,
                       provider::Union{Nothing,AbstractString}=nothing,
                       page_size::Integer=_CMR_PAGE_SIZE,
                       verbose::Bool=true, deadline::Float64=Inf,
                       requester=HTTP.request)
    1 <= page_size <= _CMR_PAGE_SIZE ||
        throw(ArgumentError("page_size must be in 1:$(_CMR_PAGE_SIZE), got $(page_size)"))

    query = Dict{String,Any}("short_name" => short_name, "version" => version,
                             "temporal" => temporal, "bounding_box" => bounding_box,
                             "provider" => provider)

    context = "CMR granule search ($(short_name))"
    out = NamedTuple{(:id, :url, :bytes),Tuple{String,String,Int}}[]
    page_num = 1
    while true
        granules = _earthdata_retry(; context="$(context) page $(page_num)",
                                    verbose=verbose, deadline=deadline) do
            # `verbose` is deliberately not forwarded: upstream passes it to HTTP.jl as
            # *curl* verbosity (a full wire dump), whereas here it gates retry logging.
            EarthData.request(EarthData.granule_url, query, EarthData.Granules.UMM_G;
                              page_num=page_num, page_size=page_size, all=false,
                              requester=_classifying_requester(requester, context))
        end

        for g in granules
            url = _granule_data_url(g)
            isnothing(url) && continue
            push!(out, (; id=_granule_id(g), url=url, bytes=_granule_bytes(g)))
        end

        # A short page is the last page. `page_size` defaults to CMR's 2000 ceiling, so this
        # loop runs once for every query this package makes (~32 hits for a tile-day).
        length(granules) < page_size && break
        page_num += 1
    end
    return out
end

# ---- local mirrors of upstream PRs 1 and 2; delete each when its PR is released ----------

# PR 1. `EarthData.request` is called rather than the exported `EarthData.granules` for one
# reason: `granules` validates its keywords against `fieldnames(GranuleRequest)`, which lists
# `polygon` but none of CMR's other spatial params, so `bounding_box=` raises
# `ArgumentError: Unknown keyword argument(s)` even though CMR itself accepts it (verified:
# HTTP 200). `request` takes the query as a plain `Dict` and applies no whitelist, so it is
# the same code path minus the validation. Everything else here — pagination, schema, error
# parsing — is upstream's.
#
# A `polygon` built from the box was the other option and is *wrong* for this use: CMR treats
# polygon edges as great-circle arcs, so at MODIS tile latitudes the enclosed area is smaller
# than the box and edge granules would be missed. `bounding_box` edges are parallels and
# meridians. Once PR 1 lands, call `EarthData.granules(; bounding_box=…)` and drop this note.
#
# `cmr_query` drops `nothing` values, so the unset keywords in `query` cost nothing.

# PR 5. `all=true` is NOT used, and must not be "simplified" back to it: upstream's
# search-after loop always leaves `page_num` in the query, and CMR rejects the combination —
# page 1 returns 200, page 2 returns **HTTP 400 `page_num is not allowed with search-after`**.
# Reproduced live. So the paging is done here with plain `page_num`, which CMR accepts (also
# verified: page 1/2 of a 32-hit query return 20 and 12 at page_size=20). At the default
# `page_size` of 2000 this is a single request anyway; the loop only exists so a query larger
# than the ceiling cannot silently truncate.

# PR 4. Upstream funnels *every* error response into `error(parse_cmr_error(r))` — a plain
# `ErrorException`, which `_earthdata_is_transient` (correctly) calls permanent. Left alone,
# a CMR 503 or 429 would abort the run instead of being retried, and the status code is gone
# by the time the exception surfaces. So the classification happens here, in the one place it
# is still visible: the requester hook upstream already provides for test injection.
#
# It only *throws* for a retryable status; a permanent error is handed back to upstream
# untouched so its `parse_cmr_error` message (which reads CMR's `errors` array) is what the
# caller sees. Retires when PR 4's transient/permanent split lands upstream.
function _classifying_requester(inner, context::AbstractString)
    return function (args...; kwargs...)
        r = inner(args...; kwargs...)
        if r.status in _EARTHDATA_TRANSIENT_STATUSES || r.status >= 500
            throw(EarthdataTransientError(context, r.status, _earthdata_body_excerpt(r),
                                          _earthdata_retry_after(r)))
        end
        return r
    end
end

"""
    _granule_id(g) -> String

The archive filename — what the legacy endpoint calls `producer_granule_id`, and what
`_modis_granule_date` / `_modis_granule_tile` parse the nominal date and tile out of.

`DataGranule.Identifiers` is consulted first, for the entry whose `IdentifierType` is
`"ProducerGranuleId"`: that is the field's authoritative home in UMM-G. `GranuleUR` is the
fallback. For MCD43A3 the two are identical (verified live), but that equality is a
per-provider convention, not a schema guarantee — reading the typed identifier means the
filter does not depend on it.
"""
function _granule_id(g)
    dg = g.DataGranule
    if !isnothing(dg) && !isnothing(dg.Identifiers)
        for ident in dg.Identifiers
            ident.IdentifierType == "ProducerGranuleId" || continue
            isnothing(ident.Identifier) && continue
            return String(ident.Identifier)
        end
    end
    return String(g.GranuleUR)
end

"""
    _granule_data_url(g) -> Union{String,Nothing}

The HTTPS download URL of a granule: the `RelatedUrls` entry whose `Type` is `"GET DATA"`.

Filtering on `Type` is what makes this exact. An MCD43A3 granule carries **eight** related
URLs, so `https_urls(g)` (unfiltered) does not identify the file — among them a DOI landing
page, the `.cmr.xml` sidecar, browse imagery and an `s3credentials` endpoint, all typed
`VIEW RELATED INFORMATION`. Note the `s3://` copy of the data is typed
`GET DATA VIA DIRECT ACCESS`, *not* `GET DATA`, so the `Type` test already excludes it; the
`https://` check is a second guard rather than the mechanism.

Upstream PR 2 replaces this with `download_url(g; type="GET DATA")`.
"""
function _granule_data_url(g)
    related = g.RelatedUrls
    isnothing(related) && return nothing
    for u in related
        u.Type == "GET DATA" || continue
        startswith(u.URL, "https://") || continue
        return String(u.URL)
    end
    return nothing
end

"""
    _granule_bytes(g) -> Int

Granule size in **bytes**, from `DataGranule.ArchiveAndDistributionInformation`, or `0` when
absent (meaning "unknown — do not verify the download against it").

`SizeInBytes` is preferred because it is unambiguous; `Size` needs `SizeUnit` to interpret.
UMM-G's own schema documentation says as much, warning that a provider reporting MegaBytes
may have used either 1000² or 1024², "and therefore there is no systematic way to know the
actual size". So this figure is only ever used with `_EARTHDATA_SIZE_TOLERANCE` slack — it
catches a truncated download, and must not be treated as an exact expected length.

Upstream PR 2 replaces this with `granule_size(g)`.
"""
function _granule_bytes(g)
    dg = g.DataGranule
    isnothing(dg) && return 0
    info = dg.ArchiveAndDistributionInformation
    (isnothing(info) || isempty(info)) && return 0
    entry = first(info)

    # Both members of the `Union{FilePackageType,FileType}` element type declare all three
    # fields, so these are plain reads rather than `hasproperty` guards.
    isnothing(entry.SizeInBytes) || return Int(entry.SizeInBytes)
    isnothing(entry.Size) && return 0
    return round(Int, Float64(entry.Size) * _size_unit_factor(entry.SizeUnit))
end

# UMM-G's `SizeUnit` enumeration. `nothing` is treated as MB because that is what the DAACs
# in use report; an unrecognised unit is a hard error rather than a silent wrong number,
# since the value feeds the truncated-download check.
function _size_unit_factor(unit)
    isnothing(unit) && return Float64(1024^2)
    u = uppercase(String(unit))
    u in ("BYTES", "B") && return 1.0
    u == "KB" && return Float64(1024)
    u == "MB" && return Float64(1024^2)
    u == "GB" && return Float64(1024^3)
    u == "TB" && return Float64(1024^4)
    throw(ArgumentError("Unrecognised UMM-G SizeUnit \"$(unit)\"; cannot convert to bytes."))
end

# ------------------------------------------------------------------------- downloading

# Tolerance on the CMR-reported size, as a fraction. `granule_size` is a rounded megabyte
# figure, so it cannot be compared exactly; this is wide enough to absorb the rounding and
# narrow enough to catch a truncated transfer (which would silently corrupt an HDF4 file,
# and GDAL's error for that is unhelpful).
const _EARTHDATA_SIZE_TOLERANCE = 0.01

"""
    _earthdata_download(url, dest; token, expected_bytes=0, verbose=true,
                        timeout=Inf, deadline=Inf) -> String

Download `url` to `dest` with an Earthdata Login bearer token, returning `dest`.

Writes to `dest * ".part"` and `mv`s on success, so an interrupted transfer never leaves a
plausible-looking truncated file in a cache. When `expected_bytes > 0` the final size is
checked against it within `_EARTHDATA_SIZE_TOLERANCE`; a short file is treated as a
transient failure and re-downloaded, because a truncated HDF4 granule fails much later and
far less legibly.

An existing `dest` is returned untouched — cache-hit handling belongs to the caller only
insofar as it chooses the path.

!!! note "The bearer header survives the redirect, and that is fine"
    LP DAAC answers a bearer-authenticated GET with a single 303 to a CloudFront URL, and
    `Downloads.download` follows it and returns the full body — verified live, with and
    without the header forwarded. There is *no* need to split this into a headerless second
    request: the presigned-URL/`Authorization` conflict that affects some S3 endpoints does
    not arise on this path. Do not add redirect handling "to be safe"; it was measured.
"""
function _earthdata_download(url::AbstractString, dest::AbstractString;
                             token::AbstractString, expected_bytes::Integer=0,
                             verbose::Bool=true, timeout::Real=Inf,
                             deadline::Float64=Inf)
    mkpath(dirname(dest))
    tmp = dest * ".part"
    headers = ["Authorization" => "Bearer $(token)"]
    limit = min(deadline, isfinite(timeout) ? time() + timeout : Inf)

    try
        _earthdata_retry(; context="download $(basename(dest))", verbose=verbose,
                         deadline=limit) do
            isfile(tmp) && rm(tmp; force=true)
            Downloads.download(url, tmp; headers=headers)
            n = filesize(tmp)
            if expected_bytes > 0 &&
               n < expected_bytes * (1 - _EARTHDATA_SIZE_TOLERANCE)
                # Retryable: a short body is a cut-off transfer, not a bad request.
                throw(EarthdataTransientError("download $(basename(dest))", 0,
                    "Got $(n) bytes, expected ≈$(expected_bytes) — transfer was truncated."))
            end
            return nothing
        end
        mv(tmp, dest; force=true)
    catch
        isfile(tmp) && rm(tmp; force=true)
        rethrow()
    end
    return dest
end
