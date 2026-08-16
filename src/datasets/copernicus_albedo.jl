"""
C3S satellite surface albedo — 10-daily gridded surface albedo derived from
Sentinel-3 OLCI + SLSTR at 300 m, from the Copernicus Climate Data Store dataset
`satellite-albedo` ("Surface albedo 10-daily gridded data from 1981 to present").

Unlike every other gridded source in this package, this product is **not** published
as a cloud-optimized store: there is no ARCO Zarr copy on ECMWF's host, no COG bucket
on AWS Open Data, and no OPeNDAP endpoint. The CDS catalogue exposes it only through
the Retrieve API, whose `satellite-albedo` process advertises `async-execute` only.

Data therefore has to be **ordered** rather than read: a job is submitted, queued,
and processed server-side (typically minutes), and the result arrives as a ZIP of
per-timestep NetCDF files. This module wraps that into an order-and-cache pipeline:

1. Work out which 10-day timesteps the requested `time_range` covers.
2. Skip any already in the local cache.
3. Order the rest — **split into chunks costing ≤ 20**, the server's hard per-request
   size limit (see below).
4. Unzip into a per-timestep cache.
5. Return a lazy `RasterSeries` over `Ti`, so pixels are only read on demand.

"Lazy" here means the returned rasters are disk-backed and the minimal subset was
requested server-side — not that the remote data is randomly accessible.

## Size limits

The CDS costing endpoint reports `cost` as the size of the **cross product**
`variable × year × month × nominal_day` — not the number of dates actually wanted —
against a `limit` of 20. Two timesteps straddling a month boundary therefore cost 4,
not 2, and a full year of one variable costs 12 × 5 = 60. Narrowing `area` reduces the
*volume* returned but **not** the cost, so a long time series is inherently many jobs:
one variable for a full year is ~3 jobs, five years ~15. Expect long runtimes for
multi-year requests.

## Layers

One ordered *variable* arrives as a NetCDF holding several *layers*: the broadband
(`_BB`), near-infrared (`_NI`) and visible (`_VI`) albedos, a `_ERR` uncertainty for
each, a `QFLAG` quality mask, and a scalar `crs`. `satellite_albedo` reads the
full-spectrum broadband layer by default; pass `layer=` to select another, and use
[`satellite_albedo_layers`](@ref) to list what a cached file contains.

## Coverage

Only the current operational Sentinel-3 era is supported: `sensor="olci_and_slstr"`,
`satellite="sentinel_3"`, `product_version="v3_1"`, 300 m, 2018–2024. The earlier
AVHRR (4 km), SPOT-VGT and PROBA-V (1 km) eras exist in the same CDS dataset but use
different grids, resolutions and nominal-day conventions; they are deliberately out of
scope here rather than silently mixed into one series.

## Prerequisite

The dataset licence must be accepted once, per account, at
https://cds.climate.copernicus.eu/datasets/satellite-albedo?tab=download
Until then CDS rejects every request with HTTP 403.
"""

using Rasters
using NCDatasets   # loads Rasters' NetCDF backend (RastersNCDatasetsExt)
using DimensionalData
using Dates
import JSON3
import SHA
import p7zip_jll

# CDS dataset (process) id.
const _ALBEDO_DATASET = "satellite-albedo"

"""
Surface-albedo variables published by the C3S product.

- `:albb_dh` — broadband directional-hemispherical (black-sky) albedo
- `:albb_bh` — broadband bi-hemispherical (white-sky) albedo
- `:alsp_dh` — spectral directional-hemispherical albedo
- `:alsp_bh` — spectral bi-hemispherical albedo

The broadband pair is what a surface energy-balance model such as GEMB wants; the
spectral pair resolves individual bands.
"""
const SATELLITE_ALBEDO_VARIABLES = (:albb_bh, :albb_dh, :alsp_bh, :alsp_dh)

# The single supported era (Sentinel-3, 300 m). Kept as separate constants rather than
# a table so adding the AVHRR/VGT eras later is an explicit, reviewed change.
const _ALBEDO_SENSOR = "olci_and_slstr"
const _ALBEDO_SATELLITE = "sentinel_3"
const _ALBEDO_VERSION = "v3_1"
const _ALBEDO_RESOLUTION = "300m"
const _ALBEDO_YEARS = 2018:2024

# Server-enforced request size limit, in cross-product units (see `_albedo_cost`).
# Confirmed against the live costing endpoint.
const _ALBEDO_COST_LIMIT = 20

# Per-job give-up time, in seconds. Far above the generic `cds_retrieve` default because
# this product queues for a long time: a 12-timestep order was measured still `running`
# after 1840 s, so a 1800 s timeout aborts jobs that would have succeeded. A timeout is
# especially expensive here — it kills the wait *before* the finished result is cached, so
# the work is lost and the order must be resubmitted for another 30–90 min of queue time.
const _ALBEDO_JOB_TIMEOUT = 10800

# Filename token per variable, as used in the C3S product filenames. Verified against a
# delivered file: `c3s_ALBB-DH_20190610000000_GLOBE_SENTINEL3_V3.1.0.area-subset.….nc`
const _ALBEDO_FILE_TOKENS = Dict(
    :albb_bh => "ALBB-BH",
    :albb_dh => "ALBB-DH",
    :alsp_bh => "ALSP-BH",
    :alsp_dh => "ALSP-DH",
)

"""
Default NetCDF layer to read for each ordered variable.

**One ordered variable is a file of several layers, not a single grid.** An `albb_dh`
file carries `AL_DH_BB` (total/broadband shortwave), `AL_DH_NI` (near-infrared) and
`AL_DH_VI` (visible), each with a matching `_ERR` uncertainty layer, plus `QFLAG` and a
scalar `crs`. The entries here name the layer that corresponds to the variable's headline
quantity — the full-spectrum broadband albedo — which is what a surface energy-balance
model wants.

This mapping matters because `Raster(path)` on such a file does **not** pick the albedo:
it selects a layer arbitrarily (in practice the scalar `crs`, giving a 1-element `Char`
raster with no `X`/`Y` dims). Layers must be named explicitly. Use `layer` to read a
different one, or [`satellite_albedo_layers`](@ref) to list what a file contains.
"""
const _ALBEDO_DEFAULT_LAYERS = Dict(
    :albb_dh => :AL_DH_BB,   # verified against a delivered file
    :albb_bh => :AL_BH_BB,   # inferred from the DH naming; not yet verified
    :alsp_dh => :AL_DH_BB,   # inferred; spectral files may name bands differently
    :alsp_bh => :AL_BH_BB,   # inferred
)

# Layers that are not geophysical fields on the X/Y grid, so never a sensible default.
const _ALBEDO_NON_DATA_LAYERS = (:crs,)

"""
    _default_albedo_cache() -> String

Default on-disk cache directory for ordered albedo files.
"""
_default_albedo_cache() = joinpath(tempdir(), "GEMB_ClimateForcing", "satellite_albedo")

"""
    _albedo_era_tag() -> String

Directory-name fragment identifying the product era, so files from different eras (or
product versions) can never collide in the cache.
"""
_albedo_era_tag() = "$(_ALBEDO_SATELLITE)_$(_ALBEDO_VERSION)_$(_ALBEDO_RESOLUTION)"

"""
    _albedo_normalize_variables(variable) -> Vector{Symbol}

Validate and canonicalise the requested variable(s).
"""
function _albedo_normalize_variables(variable)
    vars = variable isa Symbol ? [variable] : collect(variable)
    isempty(vars) && throw(ArgumentError("`variable` must name at least one albedo variable"))
    for v in vars
        v in SATELLITE_ALBEDO_VARIABLES || throw(ArgumentError(
            "Unknown albedo variable $(repr(v)). Available: $(SATELLITE_ALBEDO_VARIABLES)"))
    end
    return unique(vars)
end

"""
    _albedo_nominal_days(year, month) -> Vector{Int}

The three nominal days of the 10-daily (decadal) cycle in a given month: day 10,
day 20, and the last day of the month.

The end-of-month day varies (28/29/30/31), which `Dates.daysinmonth` handles including
leap Februaries. Verified against the live CDS constraints endpoint, which reports
`nominal_day = ["10","20","28","29","30","31"]` for the Sentinel-3 era.

```jldoctest
julia> GEMB_ClimateForcing._albedo_nominal_days(2020, 2)
3-element Vector{Int64}:
 10
 20
 29
```
"""
_albedo_nominal_days(year::Integer, month::Integer) =
    [10, 20, daysinmonth(Date(year, month))]

"""
    _albedo_timesteps(time_range; nominal_day=nothing) -> Vector{Date}

All 10-daily timestep dates falling within `time_range` (inclusive).

When `nominal_day` is given it replaces the standard `[10, 20, end-of-month]` cycle —
an escape hatch for months whose convention differs from the verified one.
"""
function _albedo_timesteps(time_range::Tuple{DateTime,DateTime};
                           nominal_day::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    t0, t1 = time_range
    t0 <= t1 || throw(ArgumentError(
        "time_range must be (start, stop) with start <= stop; got ($(t0), $(t1))"))

    dates = Date[]
    d0, d1 = Date(t0), Date(t1)
    for y in year(d0):year(d1), m in 1:12
        # Skip months entirely outside the range.
        Date(y, m, daysinmonth(Date(y, m))) < d0 && continue
        Date(y, m, 1) > d1 && continue
        days = isnothing(nominal_day) ? _albedo_nominal_days(y, m) : collect(nominal_day)
        for day in days
            day <= daysinmonth(Date(y, m)) || continue
            d = Date(y, m, day)
            d0 <= d <= d1 && push!(dates, d)
        end
    end
    return sort!(unique!(dates))
end

"""
    _albedo_validate_years(dates)

Check every timestep lies within the era's published coverage.
"""
function _albedo_validate_years(dates::AbstractVector{Date})
    isempty(dates) && throw(ArgumentError(
        "No 10-daily albedo timesteps fall within the requested time_range. Timesteps " *
        "are on day 10, day 20 and the last day of each month."))
    for d in dates
        year(d) in _ALBEDO_YEARS || throw(ArgumentError(
            "$(d) is outside the known coverage of the Sentinel-3 300 m albedo era " *
            "($(first(_ALBEDO_YEARS))–$(last(_ALBEDO_YEARS))). Earlier years are served " *
            "by the AVHRR/VGT eras, which this reader does not support."))
    end
    return nothing
end

"""
    _albedo_area(extent) -> Union{Nothing,Vector{Float64}}

Convert an `extent` to the CDS `area` field, which is ordered `[north, west, south, east]`.
Returns `nothing` for a global request.

Longitude is −180…180 here (as for the Copernicus DEM, unlike the ERA5-Land invariant
grid's 0–360°E).
"""
_albedo_area(::Nothing) = nothing
function _albedo_area(extent)
    x = extent.X
    y = extent.Y
    xmin, xmax = float(minimum(x)), float(maximum(x))
    ymin, ymax = float(minimum(y)), float(maximum(y))
    (-180.0 <= xmin && xmax <= 180.0) || throw(ArgumentError(
        "Albedo longitude must be in [-180, 180]; got X=$(x)"))
    (-90.0 <= ymin && ymax <= 90.0) || throw(ArgumentError(
        "Albedo latitude must be in [-90, 90]; got Y=$(y)"))
    return [ymax, xmin, ymin, xmax]
end

"""
    _albedo_area_tag(area) -> String

Short, stable directory name for a spatial subset, so files ordered for different
`area`s never collide (they share filenames but differ in size).
"""
_albedo_area_tag(::Nothing) = "global"
function _albedo_area_tag(area::AbstractVector{<:Real})
    rounded = round.(float.(area); digits=4)
    return first(bytes2hex(SHA.sha256(JSON3.write(rounded))), 12)
end

"""
    _albedo_cost(n_variables, dates) -> Int

The CDS request cost for a group of variables and timestep dates.

The cost is **not** `length(dates)`: a CDS request selects independent `year`, `month`
and `nominal_day` lists, and the server charges the full **cross product** of them —
including combinations that do not exist:

    cost = n_variables × |years| × |months| × |nominal_days|

This is why a single variable over one year costs 60 (12 months × 5 distinct nominal
days — a February contributes either 28 or 29, never both) rather than the 36 real
timesteps. Getting it wrong under-counts any group
straddling a month boundary, so a request that looks legal locally is rejected by the
server. Verified against the live costing endpoint across single-month, month-straddling,
full-year and multi-variable cases, so requests can be sized without a network round-trip.

```jldoctest
julia> using Dates

julia> GEMB_ClimateForcing._albedo_cost(1, [Date(2019, 6, 10), Date(2019, 6, 20)])
2

julia> GEMB_ClimateForcing._albedo_cost(1, [Date(2019, 5, 31), Date(2019, 6, 10)])
4
```
"""
function _albedo_cost(n_variables::Integer, dates::AbstractVector{Date})
    isempty(dates) && return 0
    n_years = length(unique(year(d) for d in dates))
    n_months = length(unique(month(d) for d in dates))
    n_days = length(unique(day(d) for d in dates))
    return n_variables * n_years * n_months * n_days
end

"""
    _albedo_chunk_timesteps(dates, n_variables; limit=_ALBEDO_COST_LIMIT) -> Vector{Vector{Date}}

Partition timesteps into groups each cheap enough to submit as a single job, i.e. every
group satisfies `_albedo_cost(n_variables, group) <= limit`.

Greedy over chronologically sorted dates: a date is added to the current group while it
keeps the group's cross-product cost within `limit`, otherwise a new group starts. Since
[`_albedo_cost`](@ref) charges the cross product, groups naturally break at month and
year boundaries, where the cost jumps.
"""
function _albedo_chunk_timesteps(dates::AbstractVector{Date}, n_variables::Integer;
                                 limit::Integer=_ALBEDO_COST_LIMIT)
    n_variables >= 1 || throw(ArgumentError("n_variables must be >= 1"))
    n_variables <= limit || throw(ArgumentError(
        "Cannot request $(n_variables) variables at once: the CDS size limit is " *
        "$(limit) variable × timestep combinations, so even a single timestep " *
        "exceeds it. Request fewer variables per call."))

    chunks = Vector{Vector{Date}}()
    current = Date[]
    for d in sort(collect(dates))
        if isempty(current)
            push!(current, d)
        elseif _albedo_cost(n_variables, vcat(current, d)) <= limit
            push!(current, d)
        else
            push!(chunks, current)
            current = [d]
        end
    end
    isempty(current) || push!(chunks, current)
    return chunks
end

"""
    _albedo_request(variables, dates; area=nothing) -> Dict{String,Any}

Build the CDS `inputs` dict for a group of variables and timesteps.

The request is a **cross product** of `year` × `month` × `nominal_day`, so it can cover
more timesteps than `dates` lists; extra files are discarded on the read side rather
than trying to make the request exact. Note `sensor` is a bare string while every other
field is an array of zero-padded strings.
"""
function _albedo_request(variables::AbstractVector{Symbol}, dates::AbstractVector{Date};
                         area::Union{Nothing,AbstractVector{<:Real}}=nothing)
    request = Dict{String,Any}(
        "variable" => [string(v) for v in variables],
        "satellite" => [_ALBEDO_SATELLITE],
        "sensor" => _ALBEDO_SENSOR,               # a single string, not an array
        "product_version" => [_ALBEDO_VERSION],
        "horizontal_resolution" => [_ALBEDO_RESOLUTION],
        "year" => sort(unique(string(year(d)) for d in dates)),
        "month" => sort(unique(lpad(month(d), 2, '0') for d in dates)),
        "nominal_day" => sort(unique(lpad(day(d), 2, '0') for d in dates)),
    )
    isnothing(area) || (request["area"] = collect(float.(area)))
    return request
end

"""
    _albedo_request_hash(request) -> String

Stable cache key for a request.

Uses SHA-256 over a canonical JSON form (keys sorted, array values sorted) rather than
`Base.hash`, because this key must stay identical across sessions and Julia versions for
the cache to be reused at all — unlike the DEM's VRT hash, which is regenerated on
every call and so may drift harmlessly.
"""
function _albedo_request_hash(request::AbstractDict)
    canonical = Dict{String,Any}()
    for k in sort(collect(keys(request)))
        v = request[k]
        canonical[String(k)] = v isa AbstractVector ? sort(collect(v), by=string) : v
    end
    ordered = [String(k) => canonical[String(k)] for k in sort(collect(keys(canonical)))]
    return first(bytes2hex(SHA.sha256(JSON3.write(ordered))), 16)
end

"""
    _albedo_date_from_path(path) -> Union{Nothing,Date}

Parse the timestep date from a product filename, which embeds it as `YYYYMMDD`.
Returns `nothing` when no date-like token is present, so callers can fall back to
reading the file's `time` coordinate.
"""
function _albedo_date_from_path(path::AbstractString)
    m = match(r"(\d{8})", basename(String(path)))
    m === nothing && return nothing
    s = m.captures[1]
    y, mo, d = parse(Int, s[1:4]), parse(Int, s[5:6]), parse(Int, s[7:8])
    (1 <= mo <= 12 && 1 <= d <= 31) || return nothing
    return Date(y, mo, d)
end

"""
    _albedo_variable_from_path(path) -> Union{Nothing,Symbol}

Identify which albedo variable a product file holds, from the `ALBB-DH`-style token in
its name (matched case-insensitively).
"""
function _albedo_variable_from_path(path::AbstractString)
    name = uppercase(basename(String(path)))
    for (var, token) in _ALBEDO_FILE_TOKENS
        occursin(token, name) && return var
    end
    return nothing
end

"""
    _albedo_is_zip(path) -> Bool

Whether a downloaded file is a ZIP archive, by magic bytes (`PK\\x03\\x04`). CDS returns
a ZIP for multi-file requests but may return a bare NetCDF (`CDF` or HDF5 signature) for
a single field, so the content is sniffed rather than assumed.
"""
function _albedo_is_zip(path::AbstractString)
    open(path, "r") do io
        magic = read(io, 4)
        return length(magic) == 4 && magic == UInt8[0x50, 0x4b, 0x03, 0x04]
    end
end

"""
    _albedo_extract!(archive, dest_dir) -> Vector{String}

Extract the NetCDF files from a downloaded CDS result into `dest_dir`, returning their
paths. Handles both a ZIP archive and a bare NetCDF, and walks the extraction directory
recursively since the archive may nest files inside a folder.
"""
function _albedo_extract!(archive::AbstractString, dest_dir::AbstractString)
    mkpath(dest_dir)
    if !_albedo_is_zip(archive)
        # A bare NetCDF: adopt it directly.
        target = joinpath(dest_dir, basename(archive))
        cp(archive, target; force=true)
        return [target]
    end
    p7zip_jll.p7zip() do exe
        run(pipeline(`$(exe) x -y -o$(dest_dir) $(archive)`; stdout=devnull, stderr=devnull))
    end
    files = String[]
    for (root, _, names) in walkdir(dest_dir), name in names
        endswith(lowercase(name), ".nc") && push!(files, joinpath(root, name))
    end
    isempty(files) && error("""
    No NetCDF files found after extracting the CDS result:
        $(archive)
    Extracted into $(dest_dir). The product layout may have changed.
    """)
    return files
end

"""
    _albedo_cached_files(store_dir) -> Dict{Tuple{Date,Symbol},String}

Index the per-timestep cache by `(date, variable)`.
"""
function _albedo_cached_files(store_dir::AbstractString)
    index = Dict{Tuple{Date,Symbol},String}()
    isdir(store_dir) || return index
    for name in readdir(store_dir)
        endswith(lowercase(name), ".nc") || continue
        path = joinpath(store_dir, name)
        date = _albedo_date_from_path(name)
        var = _albedo_variable_from_path(name)
        (date === nothing || var === nothing) && continue
        index[(date, var)] = path
    end
    return index
end

"""
    satellite_albedo_layers(path) -> Vector{Symbol}

The geophysical layers a downloaded C3S albedo NetCDF contains, e.g.

    [:AL_DH_BB, :AL_DH_BB_ERR, :AL_DH_NI, :AL_DH_NI_ERR, :AL_DH_VI, :AL_DH_VI_ERR, :QFLAG]

`AL_{DH|BH}_{BB|NI|VI}` is albedo over the broadband/near-infrared/visible spectrum,
`_ERR` its uncertainty, and `QFLAG` a per-pixel quality flag. The non-grid `crs` scalar
is filtered out. Pass any of these as `satellite_albedo`'s `layer` keyword.
"""
function satellite_albedo_layers(path::AbstractString)
    layers = collect(keys(RasterStack(String(path); lazy=true)))
    return [l for l in layers if !(l in _ALBEDO_NON_DATA_LAYERS)]
end

"""
    _albedo_resolve_layers(path, variable, layer) -> Vector{Symbol}

Decide which NetCDF layer(s) to read, as a vector. `layer` may be a `Symbol`, a vector of
them, or `nothing`; see [`_albedo_resolve_layer`](@ref) for the single-layer resolution
each entry goes through.
"""
function _albedo_resolve_layers(path::AbstractString, variable::Symbol, layer)
    layers = isnothing(layer) || layer isa Symbol ? [layer] : collect(layer)
    isempty(layers) && throw(ArgumentError("`layer` must name at least one layer"))
    return unique(_albedo_resolve_layer(path, variable, l) for l in layers)
end

"""
    _albedo_resolve_layer(path, variable, layer) -> Symbol

Decide which NetCDF layer to read. An explicit `layer` is validated against the file;
`nothing` falls back to the variable's headline broadband layer
([`_ALBEDO_DEFAULT_LAYERS`](@ref)), and then to the first geophysical layer if the
product's naming has changed.
"""
function _albedo_resolve_layer(path::AbstractString, variable::Symbol,
                               layer::Union{Nothing,Symbol})
    available = satellite_albedo_layers(path)
    if !isnothing(layer)
        layer in available || throw(ArgumentError(
            "Layer $(repr(layer)) is not in the $(variable) albedo file. " *
            "Available: $(join(available, ", "))"))
        return layer
    end
    preferred = get(_ALBEDO_DEFAULT_LAYERS, variable, nothing)
    (!isnothing(preferred) && preferred in available) && return preferred
    isempty(available) && error("""
    The albedo file contains no readable geophysical layer:
        $(path)
    The product layout may have changed.
    """)
    @warn """
    Expected layer $(repr(preferred)) is absent from the $(variable) albedo file; \
    falling back to $(repr(first(available))). Pass `layer` to choose explicitly.
    """ path available
    return first(available)
end

"""
    _albedo_open_layer(path, layer) -> Raster

Open one layer of a product file as a lazy 2-D (`X`, `Y`) `Raster`.

The stored files carry a singleton `time` axis; it is dropped here because time is
carried by the enclosing `RasterSeries`' `Ti` dimension instead (the same reasoning as
`_open_invariant_raster`). Packed integers are unpacked by NCDatasets via `scale_factor`
even under `lazy=true`, so values come back as fractions in 0–1.
"""
function _albedo_open_layer(path::AbstractString, layer::Symbol)
    r = Raster(String(path); name=layer, lazy=true)
    return hasdim(r, Ti) ? view(r, Ti(1)) : r
end

"""
    _albedo_check_time_coord(path, layer, expected_date; verbose)

Cross-check a file's internal `time` coordinate against the date parsed from its
filename, warning on mismatch. Run on the first file only — cheap insurance that the
filename convention is what we assume, since the `Ti` lookup depends on it.
"""
function _albedo_check_time_coord(path::AbstractString, layer::Symbol,
                                  expected_date::Date; verbose::Bool)
    verbose || return nothing
    try
        r = Raster(String(path); name=layer, lazy=true)
        hasdim(r, Ti) || return nothing
        actual = Date(first(lookup(r, Ti)))
        actual == expected_date || @warn """
        Albedo filename date disagrees with the file's time coordinate — the `Ti` lookup \
        may be wrong.
        """ path filename_date=expected_date time_coord=actual
    catch e
        @debug "Could not read time coordinate for cross-check" path exception=e
    end
    return nothing
end

"""
    _run_concurrent_jobs(f, n_jobs, max_concurrent, submit_stagger; verbose=true)

Run `f(i, submit_gate)` for `i in 1:n_jobs` concurrently, with at most `max_concurrent`
in flight and successive calls to `submit_gate()` spaced at least `submit_stagger`
seconds apart.

**Why this exists.** A CDS job spends nearly all its wall-clock queued and processed
server-side — tens of minutes, and highly variable (a 3-timestep order has been observed
finishing in ~30 min while another ran past 2 h). Ordering the jobs one after another
therefore costs the *sum* of those waits, while CDS will work on several at once.

Measured, six one-month orders submitted together: three finished at 27–32 min, the rest
at ~79–80 min, all six done in 83 min — versus ~30 min *each*, i.e. ~3 h, in series. So
the server clearly runs only ~3 of our jobs at a time and the speedup is ~2×, not linear
in `max_concurrent`; raising the cap past a handful buys little and risks the rejection
path. Concurrency, not smaller jobs, is nonetheless what makes a multi-month request
tractable.

Tasks rather than threads: the work is HTTP polling, so it is entirely I/O-bound and
`@async` imposes no thread-safety requirement on HTTP.jl or the NetCDF/`p7zip` calls.

The two knobs both exist to stay inside CDS's limits:
- `max_concurrent` bounds in-flight jobs, because CDS *rejects* (rather than queues)
  submissions past a per-dataset limit.
- `submit_stagger` spaces submissions out, because a simultaneous burst trips that
  limiter even at a count that succeeds when spread out. `f` decides when to open the
  gate, so a job served from cache does not consume a submission slot at all, and only
  the *waiting* is serialised — polling stays parallel.

Errors propagate: `@sync` rethrows, so a failed job is not silently dropped.
"""
function _run_concurrent_jobs(f, n_jobs::Integer, max_concurrent::Integer,
                              submit_stagger::Real; verbose::Bool=true)
    n_jobs >= 1 || return nothing
    limit = max(1, min(max_concurrent, n_jobs))
    slots = Base.Semaphore(limit)
    submit_lock = ReentrantLock()
    # Wall-clock of the last submission. Tracking the real time (rather than sleeping once
    # per job) means a task that queued behind the semaphore for an hour submits at once
    # instead of sitting out another `submit_stagger` seconds for nothing. Initialised in
    # the past so the first job submits immediately.
    last_submit = Ref(time() - submit_stagger)

    verbose && n_jobs > 1 &&
        @info "Running up to $(limit) CDS jobs concurrently" n_jobs = n_jobs

    submit_gate = () -> lock(submit_lock) do
        wait_s = submit_stagger - (time() - last_submit[])
        wait_s > 0 && sleep(wait_s)
        last_submit[] = time()
    end

    @sync for i in 1:n_jobs
        @async begin
            Base.acquire(slots)
            try
                f(i, submit_gate)
            finally
                Base.release(slots)
            end
        end
    end
    return nothing
end

"""
    satellite_albedo(; time_range, extent=nothing, variable=:albb_dh, layer=nothing,
                       nominal_day=nothing, token=nothing, cache_path=nothing,
                       force_download=false, verbose=true, poll_interval=10,
                       timeout=$(_ALBEDO_JOB_TIMEOUT)) -> RasterSeries

Surface albedo from the C3S `satellite-albedo` product (Sentinel-3 OLCI+SLSTR, 300 m,
10-daily), as a **lazy** `RasterSeries` over the `Ti` (time) dimension.

Data is *ordered* from the CDS Retrieve API, not streamed: each call submits one or
more jobs, waits for them, caches the resulting NetCDFs, and opens them lazily. **Expect
minutes of latency** (a single timestep took ~10 min in testing), and see the size note
below before requesting long series.

# Keywords
- `time_range::Tuple{DateTime,DateTime}`: inclusive start/stop. Timesteps fall on day 10,
  day 20 and the last day of each month; only those inside the range are returned.
- `extent`: an `Extents.Extent` or `(; X, Y)` NamedTuple in −180…180 longitude, or
  `nothing` (default) for global. Global 300 m data is ~120960 × 47040 pixels *per
  variable per timestep*, so an extent is strongly recommended.
- `variable`: one `Symbol` (default `:albb_dh`, broadband black-sky) or a vector of them.
  See [`SATELLITE_ALBEDO_VARIABLES`](@ref). One variable and one layer returns a series of
  `Raster`s; several of either return a series of `RasterStack`s.
- `layer`: which NetCDF layer(s) to read from each ordered file — one `Symbol`, a vector
  of them, or `nothing` for the variable's full-spectrum broadband albedo (`AL_DH_BB` /
  `AL_BH_BB`). Each file also holds near-infrared (`_NI`) and visible (`_VI`) bands,
  per-band `_ERR` uncertainties, and `QFLAG`; list them with
  [`satellite_albedo_layers`](@ref). Several layers of one variable come back as a
  `RasterStack` keyed by *layer* name, from the same cached file and with **no extra CDS
  job** — much cheaper than one call per layer.
- `nominal_day`: override the standard `[10, 20, end-of-month]` cycle.
- `token`: CDS API key; defaults to [`get_cds_api_key()`](@ref).
- `cache_path`: cache directory; defaults to a per-user temp directory. Cached timesteps
  are reused without submitting a job.
- `force_download`: re-order timesteps already cached.
- `poll_interval`, `timeout`: job polling cadence and give-up time, in seconds. `timeout`
  is **per job**, measured from that job's own submission. With more chunks than
  `max_concurrent_jobs` the later ones start their clock only once a slot frees, so the
  call as a whole can exceed `timeout`.
- `max_concurrent_jobs`: how many CDS jobs to keep in flight (default 6). Jobs spend
  nearly all their time waiting server-side, so this is the main throughput lever — though
  CDS appears to process only ~3 of an account's jobs at a time, so the gain saturates
  (~2× measured, not linear). CDS also rejects submissions past a per-dataset queue limit,
  hence the cap.
- `submit_stagger`: seconds between successive submissions (default 15). A simultaneous
  burst trips the queue limiter even at counts that succeed when spread out.

# Size limits
CDS charges one unit per (variable × year × month × nominal-day) combination and rejects
requests above 20. Narrowing `extent` shrinks the download but **not** the cost, so
requests are split automatically into as many jobs as needed — one variable for a year is
~3 jobs, five years ~15. Those jobs are submitted and polled **concurrently**
(`max_concurrent_jobs`), which roughly halves wall-clock versus ordering them in series,
but expect tens of minutes regardless: a single job routinely queues that long, and CDS
runs only a few of an account's jobs at once.

# Prerequisite
Accept the dataset licence once at
<https://cds.climate.copernicus.eu/datasets/satellite-albedo?tab=download#manage-licences>;
until then CDS answers every request with HTTP 403.

# Example
```julia
using Dates, Rasters
alb = satellite_albedo(;
    time_range = (DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
    extent = Extent(X = (-50.0, -45.0), Y = (66.0, 68.0)),
    variable = :albb_dh,
)
lookup(alb, Ti)                 # the 10-daily timestep dates
read(alb[1])                    # materialise the first timestep (albedo fraction, 0–1)

# The visible-band albedo from the same files instead of broadband:
vis = satellite_albedo(; time_range = (DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
                       extent = Extent(X = (-50.0, -45.0), Y = (66.0, 68.0)),
                       layer = :AL_DH_VI)

# Albedo, its uncertainty and the quality flag together — one call, one set of jobs:
qc = satellite_albedo(; time_range = (DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
                      extent = Extent(X = (-50.0, -45.0), Y = (66.0, 68.0)),
                      layer = [:AL_DH_BB, :AL_DH_BB_ERR, :QFLAG])
qc[1][:QFLAG]                   # a RasterStack keyed by layer name per timestep
```
"""
function satellite_albedo(;
    time_range::Tuple{DateTime,DateTime},
    extent=nothing,
    variable::Union{Symbol,AbstractVector{Symbol}}=:albb_dh,
    layer::Union{Nothing,Symbol,AbstractVector{Symbol}}=nothing,
    nominal_day::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    token::Union{String,Nothing}=nothing,
    cache_path::Union{String,Nothing}=nothing,
    force_download::Bool=false,
    verbose::Bool=true,
    poll_interval::Real=10,
    timeout::Real=_ALBEDO_JOB_TIMEOUT,
    max_concurrent_jobs::Integer=6,
    submit_stagger::Real=15,
)
    variables = _albedo_normalize_variables(variable)
    dates = _albedo_timesteps(time_range; nominal_day=nominal_day)
    _albedo_validate_years(dates)

    area = _albedo_area(extent)
    isnothing(area) && verbose && @warn """
    Ordering global 300 m albedo (~120960 × 47040 pixels per variable per timestep). \
    Pass `extent` to subset server-side and download far less.
    """

    cache = isnothing(cache_path) ? _default_albedo_cache() : cache_path
    store_dir = joinpath(cache, "files", _albedo_era_tag(), _albedo_area_tag(area))
    mkpath(store_dir)

    # Which (date, variable) pairs are still missing? When nothing is missing we submit
    # no job at all — this is what makes repeat calls cheap.
    cached = _albedo_cached_files(store_dir)
    wanted = [(d, v) for d in dates for v in variables]
    if force_download
        for (d, v) in wanted
            haskey(cached, (d, v)) && rm(cached[(d, v)]; force=true)
        end
        cached = _albedo_cached_files(store_dir)
    end
    missing_dates = sort(unique(d for (d, v) in wanted if !haskey(cached, (d, v))))

    if !isempty(missing_dates)
        tok = isnothing(token) ? get_cds_api_key() : token
        chunks = _albedo_chunk_timesteps(missing_dates, length(variables))
        verbose && @info "Ordering albedo from CDS" n_timesteps=length(missing_dates) n_variables=length(variables) n_jobs=length(chunks)

        n_jobs = length(chunks)
        completed = Threads.Atomic{Int}(0)

        _run_concurrent_jobs(n_jobs, max_concurrent_jobs, submit_stagger;
                             verbose=verbose) do i, submit_gate
            chunk = chunks[i]
            request = _albedo_request(variables, chunk; area=area)
            job_dir = joinpath(cache, "jobs", _albedo_request_hash(request))
            mkpath(job_dir)
            archive = joinpath(job_dir, "download")

            if force_download || !isfile(archive)
                write(joinpath(job_dir, "request.json"), JSON3.write(request))
                submit_gate()
                cds_retrieve(_ALBEDO_DATASET, request, archive;
                             token=tok, poll_interval=poll_interval,
                             timeout=timeout, verbose=verbose)
            end

            extract_dir = joinpath(job_dir, "extract")
            for file in _albedo_extract!(archive, extract_dir)
                target = joinpath(store_dir, basename(file))
                isfile(target) && !force_download ? rm(file; force=true) :
                    mv(file, target; force=true)
            end
            # The archive and scratch extraction are large and no longer needed; the
            # per-timestep store is now authoritative.
            rm(extract_dir; recursive=true, force=true)
            rm(archive; force=true)

            n = Threads.atomic_add!(completed, 1) + 1
            verbose && @info "CDS albedo job $(n)/$(n_jobs) complete" dates=chunk
            return nothing
        end
        cached = _albedo_cached_files(store_dir)
    end

    # Keep only the requested timesteps: the CDS request is a year × month × nominal_day
    # cross product, so extra files may have been returned at the range endpoints.
    available = sort(unique(d for (d, v) in keys(cached) if d in dates && v in variables))
    isempty(available) && error("""
    CDS returned no usable albedo files for the requested time range.
    Cached files were expected in:
        $(store_dir)
    This can mean the product filename convention differs from the one assumed by
    `_albedo_date_from_path` / `_albedo_variable_from_path`.
    """)

    if verbose && length(available) < length(dates)
        @warn "Some requested albedo timesteps were not returned by CDS" requested=length(dates) returned=length(available) missing_dates=setdiff(dates, available)
    end

    # Each ordered variable is a *multi-layer* file (broadband/NIR/visible albedo, their
    # uncertainties, QFLAG, and a scalar `crs`), so the layer has to be named explicitly —
    # an unnamed `Raster(path)` picks `crs` and yields a dimensionless Char. Resolve once
    # per variable from its first file; all timesteps of a variable share a layout.
    resolved = Dict(v => _albedo_resolve_layers(cached[(first(available), v)], v, layer)
                    for v in variables)
    multi_layer = length(resolved[first(variables)]) > 1

    _albedo_check_time_coord(cached[(first(available), first(variables))],
                             first(resolved[first(variables)]), first(available);
                             verbose=verbose)

    meta = Dict{String,Any}(
        "source" => "C3S satellite surface albedo (CDS dataset $(_ALBEDO_DATASET))",
        "cds_dataset" => _ALBEDO_DATASET,
        "sensor" => _ALBEDO_SENSOR,
        "satellite" => _ALBEDO_SATELLITE,
        "product_version" => _ALBEDO_VERSION,
        "resolution" => _ALBEDO_RESOLUTION,
        "variable" => length(variables) == 1 ? first(variables) : Tuple(variables),
        "layer" => let ls = length(variables) == 1 ? resolved[first(variables)] :
                            [l for v in variables for l in resolved[v]]
            length(ls) == 1 ? only(ls) : Tuple(ls)
        end,
        "extent" => isnothing(extent) ? "global" : string(extent),
        "n_timesteps" => length(available),
        "units" => "1 (dimensionless reflectance fraction, 0–1)",
        "long_name" => "Surface albedo",
    )

    # Build the lazy series: one Raster per timestep in the single-variable, single-layer
    # case, otherwise a RasterStack (mirrors `climate_model_invariant`).
    #
    # `RasterSeries` has no metadata field — `rebuild(series; metadata=...)` is accepted
    # but silently discarded — so the provenance Dict is attached to each layer instead,
    # where `view` preserves it through the lazy crop.
    layers = map(available) do d
        layer_meta = merge(meta, Dict{String,Any}("time" => d))
        # A stack is keyed by the *ordered variable* name when one layer per variable was
        # requested, so `stack[:albb_dh]` works regardless of which band was read — but by
        # *layer* name when several layers of one variable were asked for, since that is
        # what distinguishes them.
        pairs = Pair{Symbol,Any}[]
        for v in variables
            for l in resolved[v]
                r = _albedo_open_layer(cached[(d, v)], l)
                # `"file"` is the product NetCDF this layer reads. Callers needing
                # something the `Raster` does not expose — CF attributes such as the
                # `QFLAG` legend — can open it directly instead of re-deriving the cache
                # layout from private helpers. Exact per (date, variable), unlike a
                # directory scan.
                r = rebuild(r; metadata=merge(layer_meta,
                                              Dict{String,Any}("file" => cached[(d, v)],
                                                               "layer" => l)))
                key = multi_layer ? l : v
                push!(pairs, key => _crop_to_extent(r, extent))
            end
        end
        length(pairs) == 1 && return last(only(pairs))
        return RasterStack(NamedTuple(pairs); metadata=layer_meta)
    end

    return RasterSeries(layers, Ti(available))
end
