"""
Glacier bare-ice albedo at a **list of points**, from MODIS MCD43A3 (500 m, daily), in
black-sky and white-sky forms.

The statistic is identical to the C3S path's ([`compute_glacier_ice_albedo`](@ref)): per
cell and calendar year, the mean of the darkest `percentile` of valid retrievals, because
on a glacier the annual albedo minimum *is* the bare-ice state. Only the source and the
sampling geometry differ.

!!! note "Points, not a grid — and deliberately so"
    Only glacier surfaces are of interest, and a glacier point list is a vanishingly small
    fraction of any MODIS tile. Deduplicating the points to unique `(h, v, row, col)` cells
    means the accumulator is a *vector* of length `n_unique` rather than a multi-megapixel
    grid, which removes every piece of large-grid machinery the C3S path needs: no VRT
    mosaic, no reprojection, no `Mmap` scratch, no blockwise folding, no fold checkpoint.

    Do **not** add a gridded mode by warping the granules onto a lat/lon grid. An
    `:average` (or any) resampling of an individual day mixes neighbouring pixels' albedo
    *before* the percentile is taken, which biases the darkest tail toward the local mean —
    exactly the quantity being measured. Sampling the native cell is strictly more faithful.

!!! warning "Volume is the cost, and it is not reducible"
    There is no server-side subsetting: the unit of transfer is a whole ~50–95 MB granule,
    of which ~85 % is layers never read. Cost ≈ `n_tiles × n_dates × 70 MB` — a Greenland
    point list spanning 8 tiles over a full year is ~200 GB. `doy_range` is the effective
    fix, not `stride`: at high latitude there is no usable retrieval through polar night, so
    `doy_range=(120, 290)` removes ~2.5× of the download at approximately zero cost in
    surviving samples. `keep_granules=false` (the default, *inverted* from the C3S path)
    keeps peak disk at one date's tiles, because a re-download here is bandwidth-bound
    minutes rather than the C3S path's hours of queue latency.
"""

import SHA   # short stable key for the point set, so a cached sample cannot be misapplied

# Number of valid retrievals below which a cell-year is reported as NaN. Daily cadence, so
# NOT the C3S path's 10: ten usable days out of 365 is a pathological cell, and the darkest
# 5 % of ten observations is noise.
const _MODIS_ICE_MIN_SAMPLES = 30

# MCD43A3 mandatory-quality classes accepted by default. `0` is a full BRDF inversion; `1`
# (magnitude inversion) is a genuinely weaker retrieval and is opt-in.
const MCD43A3_QA_KEEP = [0]

"""
    _modis_dedup_points(lat, lon) -> (cells, cell_of_point)

Map a point list to unique MCD43A3 grid cells.

Returns the unique `(h, v, row, col)` cells in **first-encounter order** and, per input
point, the index of its cell. This is what makes multiple points in one 500 m cell cost one
sample rather than several: the accumulator is indexed by cell, and the result is scattered
back over points at the end, so the dedup is invisible to the caller except through the
reported `cell_id`.
"""
function _modis_dedup_points(lat::AbstractVector{<:Real}, lon::AbstractVector{<:Real})
    n = length(lat)
    cells = NTuple{4,Int}[]
    index = Dict{NTuple{4,Int},Int}()
    cell_of_point = Vector{Int}(undef, n)
    for i in 1:n
        cell = _modis_cell(lat[i], lon[i])
        cell_of_point[i] = get!(index, cell) do
            push!(cells, cell)
            return length(cells)
        end
    end
    return cells, cell_of_point
end

# Melt-season day-of-year windows, by hemisphere. Generous rather than tight: the low
# percentile discards bright days anyway, so the only cost of a wide window is download
# volume, while a window that is too narrow discards the darkest days themselves.
const _MODIS_MELT_SEASON_NORTH = (120, 290)   # 30 Apr – 17 Oct
const _MODIS_MELT_SEASON_SOUTH = (300, 110)   # 27 Oct – 20 Apr, wrapping New Year

"""
    _melt_season_doy_range(lat) -> Tuple{Int,Int}

Melt-season day-of-year window for a latitude: `$(_MODIS_MELT_SEASON_NORTH)` in the northern
hemisphere, `$(_MODIS_MELT_SEASON_SOUTH)` (wrapping) in the southern.

A hemisphere is a property of each *point*, not of the request, so a point list straddling
the equator has no single window. The caller resolves that by taking the union of the two,
which for a mixed list is the whole year — correct, if expensive, and better than silently
sampling the wrong season for half the points.
"""
_melt_season_doy_range(lat::Real) =
    lat >= 0 ? _MODIS_MELT_SEASON_NORTH : _MODIS_MELT_SEASON_SOUTH

"""
    _resolve_doy_range(doy_range, lat) -> Union{Nothing,Tuple{Int,Int}}

Turn a user-supplied `doy_range` into an explicit window: `nothing` and an explicit tuple
pass through, `:melt_season` resolves against the latitudes in `lat`.

For a point list confined to one hemisphere this is that hemisphere's window. For one
straddling the equator the two windows overlap only partially, so `nothing` (the whole year)
is returned rather than a window that would be wrong for some of the points.
"""
function _resolve_doy_range(doy_range, lat::AbstractVector{<:Real})
    doy_range === :melt_season || return doy_range
    north = any(>=(0), lat)
    south = any(<(0), lat)
    north && south && return nothing
    return _melt_season_doy_range(north ? 1.0 : -1.0)
end

"""
    _modis_doys(doy_range, ndays) -> Vector{Int}

Days of year selected by `doy_range` within a year of `ndays` days, in ascending order.

`lo > hi` means a window that **wraps New Year** — `(300, 110)` is 27 Oct → 20 Apr, the
southern-hemisphere melt season. Days are returned ascending (i.e. `1:hi` then `lo:ndays`)
rather than in seasonal order, because a calendar year's samples are pooled into one
statistic and the accumulator is order-independent; ascending keeps the progress bar and the
per-date cache in date order.

A wrapping window therefore samples the *tail* of one melt season and the *start* of the
next within the same calendar year. For the darkest-percentile statistic that is acceptable —
both tails are bare-ice states of the same glacier a year apart — and it avoids reworking the
whole pipeline onto melt-year buckets.
"""
function _modis_doys(doy_range::Union{Nothing,Tuple{Integer,Integer}}, ndays::Integer)
    isnothing(doy_range) && return collect(1:ndays)
    lo, hi = doy_range
    lo <= hi && return collect(lo:min(hi, ndays))
    # Wrapping: the two arcs cannot overlap, since validation enforces lo > hi.
    return vcat(1:min(hi, ndays), lo:ndays)
end

"""
    _modis_dates(year, doy_range, stride) -> Vector{Date}

The dates to sample for one calendar year: every day, restricted to `doy_range`
(inclusive day-of-year, `lo > hi` wrapping New Year, or `nothing` for the whole year),
thinned by `stride`, and clipped to the start of the record (2000-02-16).

`stride` thins the selected days *after* the window is applied, so it keeps every
`stride`-th selected day rather than every `stride`-th day of the year — for a wrapping
window those differ.

Can legitimately be empty — a `doy_range` before mid-February of 2000 selects nothing — and
the caller validates that rather than this returning a silently short list.
"""
function _modis_dates(year::Integer, doy_range::Union{Nothing,Tuple{Integer,Integer}},
                      stride::Integer)
    doys = _modis_doys(doy_range, Dates.daysinyear(year))
    isempty(doys) && return Date[]
    dates = [Date(year) + Day(d - 1) for d in doys[1:stride:end]]
    return filter(>=(_MCD43A3_FIRST_DATE), dates)
end

# ------------------------------------------------------------------- per-date sampling

"""
    _modis_sample_cells(paths, cells, cell_of_tile, layers) -> Dict{Symbol,Vector}

Sample one date's granules at the requested cells, returning per layer a vector of length
`length(cells)` of **raw** values (`Int16` albedo, `UInt8` quality), plus the quality
layers, with `missing` where no granule covered the cell.

Reads the *bounding window* of each tile's needed rows and columns once per layer rather
than one pixel at a time: a whole 2400×2400 `Int16` layer is only 11 MB, so even scattered
points are cheaper to serve from a single windowed read than from many small ones.
"""
function _modis_sample_cells(paths::Dict{Tuple{Int,Int},String},
                             cells::Vector{NTuple{4,Int}},
                             cell_of_tile::Dict{Tuple{Int,Int},Vector{Int}},
                             layers)
    quality = unique(_mcd43a3_quality_layer.(layers))
    wanted = vcat(collect(layers), quality)
    out = Dict{Symbol,Vector{Union{Missing,Int}}}(
        l => Vector{Union{Missing,Int}}(missing, length(cells)) for l in wanted)

    for (tile, path) in paths
        idx = get(cell_of_tile, tile, Int[])
        isempty(idx) && continue
        rows = extrema(cells[i][3] for i in idx)
        cols = extrema(cells[i][4] for i in idx)
        rrange, crange = rows[1]:rows[2], cols[1]:cols[2]
        for layer in wanted
            window = _modis_read_window(path, layer, rrange, crange)
            dest = out[layer]
            for i in idx
                _, _, row, col = cells[i]
                dest[i] = Int(window[row - rrange.start + 1, col - crange.start + 1])
            end
        end
    end
    return out
end

# A short, stable key for the point set, so a cached per-date sample cannot be reused for a
# different list of cells. Hashing the cells (not the input lat/lon) is what makes two
# callers whose points land on the same cells share a cache entry.
function _modis_cells_key(cells::Vector{NTuple{4,Int}})
    buf = IOBuffer()
    for c in cells
        print(buf, c[1], ",", c[2], ",", c[3], ",", c[4], ";")
    end
    return bytes2hex(SHA.sha1(take!(buf)))[1:16]
end

"""
    _modis_sample_cache_path(cache, cells_key, date) -> String

Where one date's sampled cell values are cached.

Caching the *samples* rather than checkpointing an accumulator is the right shape here: a
date's samples are a few kB, so a re-run with the same point list never re-downloads and a
run killed mid-year resumes for free — without any of the C3S path's mmap-scratch and
fold-state machinery, which exists only because a global grid cannot be held in RAM.
"""
_modis_sample_cache_path(cache::AbstractString, cells_key::AbstractString, date::Date) =
    joinpath(cache, "samples", cells_key, string(date, ".tsv"))

# Samples are stored as a tiny TSV: a header of layer names, then one line per cell with
# integer values or "NA". Plain text on purpose — these files outlive package versions, and
# a serialization format that a future `_valid_albedo` change could silently misread is a
# worse trade than a few kB of disk.
function _modis_sample_cache_write(path::AbstractString, layers::Vector{Symbol},
                                   samples::Dict{Symbol,Vector{Union{Missing,Int}}})
    mkpath(dirname(path))
    tmp = path * ".part"
    open(tmp, "w") do io
        println(io, join(layers, "\t"))
        n = length(samples[first(layers)])
        for i in 1:n
            println(io, join((ismissing(samples[l][i]) ? "NA" : string(samples[l][i])
                              for l in layers), "\t"))
        end
    end
    mv(tmp, path; force=true)
    return path
end

function _modis_sample_cache_read(path::AbstractString, layers::Vector{Symbol}, ncell::Int)
    isfile(path) || return nothing
    lines = readlines(path)
    isempty(lines) && return nothing
    header = Symbol.(split(lines[1], '\t'))
    # A cached file written for a different layer set (or a different cell count) is
    # ignored, not partially trusted — mixing the two would sample the wrong layer.
    issubset(layers, header) || return nothing
    length(lines) == ncell + 1 || return nothing
    out = Dict{Symbol,Vector{Union{Missing,Int}}}(
        l => Vector{Union{Missing,Int}}(missing, ncell) for l in layers)
    cols = Dict(l => findfirst(==(l), header) for l in layers)
    for i in 1:ncell
        fields = split(lines[i + 1], '\t')
        for l in layers
            f = fields[cols[l]]
            out[l][i] = f == "NA" ? missing : parse(Int, f)
        end
    end
    return out
end

# --------------------------------------------------------------------------- the driver

"""
    compute_glacier_ice_albedo_modis(lat, lon, years = MCD43A3_YEARS; kwargs...) -> DimStack
    compute_glacier_ice_albedo_modis(points, years = MCD43A3_YEARS; kwargs...) -> DimStack

Bare-ice albedo at each point in `lat`/`lon` (or a vector of `(lat, lon)` tuples), per
calendar year in `years`, from MODIS MCD43A3 v061 — as **black-sky** (`:albedo_bsa`,
directional-hemispherical) and **white-sky** (`:albedo_wsa`, bihemispherical) albedo.

Per cell and year, the reported value is the mean of the darkest `percentile` of that year's
valid retrievals. On a glacier the annual minimum is the bare-ice state, and averaging a low
percentile rather than taking the single minimum stops one bad retrieval from setting the
answer.

# Returns

A `DimStack` over `(Dim{:point}(1:n), Ti(...))` in **input point order**, with layers:

- `:albedo_bsa`, `:albedo_wsa` — `Float32`, `NaN32` where unresolved (all retrievals
  rejected by QC, or fewer than `min_samples` of them).
- `:n_valid_observations_bsa`, `:n_valid_observations_wsa` — `Int32`. Reported even where
  the albedo is `NaN`, so a missing value can be diagnosed: zero means QC rejected
  everything, nonzero-but-small means `min_samples` bit.
- `:latitude`, `:longitude` — `Dim{:point}`-only, the **cell centre** actually sampled (not
  the requested point), so the sampling location is auditable.
- `:cell_id` — `Dim{:point}`-only string `"h16v02_r0699_c0124"`. **Two points sharing this
  necessarily share their albedo**, because the value was derived once for the cell.

With `layers` other than the default pair, layer names follow the same scheme lowercased —
e.g. `:Albedo_BSA_vis` → `:albedo_bsa_vis` and `:n_valid_observations_bsa_vis`.

# Arguments
- `lat`, `lon`: point coordinates in degrees; longitude may be −180…180 or 0…360.
- `years`: calendar years, within [`MCD43A3_YEARS`](@ref). 2000 is partial (the record
  starts 2000-02-16), so relax `min_samples` for it.

# Keywords — the statistic
- `percentile = 0.05`: fraction of the darkest valid retrievals averaged. Unchanged from the
  C3S path even though the cadence is 36× denser: 5 % of ~365 days retains ~19 values, which
  is *more* robust than C3S's darkest 1–2, not less.
- `min_samples = $(_MODIS_ICE_MIN_SAMPLES)`: minimum valid retrievals for a cell-year to be
  reported. Higher than the C3S path's 10 because of the daily cadence.
- `albedo_range = (0.3, 1.0)`: per-observation accept range. **The 0.3 floor is
  glaciological, not physical** — exposed ice rarely goes below ~0.3 broadband, so darker
  pixels are usually rock, water, shadow or a failed inversion. Heavily dust- or
  algae-darkened ablation zones need it lowered (~0.15) or the default clips the very signal
  being measured. The upper bound is QC too: MCD43A3's valid range reaches 32.766 after
  scaling, so albedo *can* legitimately exceed 1.0 and no physical clamp is applied.

# Keywords — quality control
- `qa_keep = $(MCD43A3_QA_KEEP)`: accepted `BRDF_Albedo_Band_Mandatory_Quality_*` classes,
  as a **whitelist**. `0` is a full BRDF inversion; add `1` to admit magnitude inversions
  (weaker, but they roughly double the sample count at high latitude). Classes `2`–`7` are
  v061 Terra detector-failure cases and `255` is fill; a whitelist rejects those *and* any
  class a later product version adds, which a reject-list cannot.

  Note there is no `max_error` equivalent: MCD43A3 carries no per-pixel uncertainty layer,
  and its snow flags live in a separate granule (MCD43A2) that is out of scope. Bright,
  snow-covered days are simply discarded by the low percentile, exactly as in the C3S path.

# Keywords — download volume
- `doy_range = nothing`: inclusive day-of-year window, e.g. `(120, 290)`. **The recommended
  volume knob** — see the module note; at high latitude it removes ~2.5× of the download at
  approximately zero cost in surviving samples.

  `first > last` means a window that **wraps New Year**, which is what a southern-hemisphere
  melt season needs: `(300, 110)` is 27 Oct → 20 Apr. A wrapping window pools the tail of one
  melt season with the start of the next inside the same calendar year — acceptable for a
  darkest-percentile statistic, since both are bare-ice states of the same glacier.

  `:melt_season` picks the window from the points' latitudes:
  `$(_MODIS_MELT_SEASON_NORTH)` north of the equator, `$(_MODIS_MELT_SEASON_SOUTH)` south of
  it. A point list straddling the equator has no single melt season, so it falls back to the
  whole year rather than sampling the wrong one for half the points — split such a list by
  hemisphere and call twice.
- `stride = 1`: sample every `stride`-th day. The blunt fallback; unlike `doy_range` it
  thins the melt season too.

# Keywords — access
- `layers = MCD43A3_ALBEDO_LAYERS`: albedo layers to reduce, from [`MCD43A3_LAYERS`](@ref).
  Each is reduced independently from the same granules, so adding layers costs no extra
  download.
- `token = nothing`: Earthdata Login token; defaults to [`get_earthdata_token`](@ref).
- `cache_path = nothing`: cache root; defaults to a per-user temp directory, which the OS may
  clear. Per-date sampled values are cached under it, so a second identical call is
  near-instant and bit-identical.
- `keep_granules = false`: keep the downloaded `.hdf` granules. **Inverted from the C3S
  path** — see the module note. `true` makes a re-run with *different* QC settings free at
  the cost of ~70 MB × tiles × dates of disk.
- `force_download = false`: ignore both caches and re-fetch.
- `max_concurrent_downloads = 4`: parallel granule downloads. Bandwidth-bound, so there is
  no queue limiter to trip and no submission stagger — unlike the CDS path.
- `download_timeout = 3600`: seconds per granule.
- `progress = true`, `verbose = true`: plain-text progress bar / info logging.

# Examples
```julia
# Two points on Russell Glacier's ablation zone, melt season only.
lat = [67.09, 66.95]
lon = [-50.05, -49.90]
ice = compute_glacier_ice_albedo_modis(lat, lon, 2019:2020; doy_range=(150, 260))
ice[:albedo_bsa]                       # 2 × 2 (point × year)
ice[:cell_id]                          # which 500 m cell each point sampled

# Darkened ice, and admitting magnitude inversions to double the sample count.
ice = compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                       albedo_range=(0.15, 1.0), qa_keep=[0, 1])
```

# See Also
[`compute_glacier_ice_albedo`](@ref) — the same statistic from the C3S 300 m record, over a
gridded extent rather than a point list.
"""
function compute_glacier_ice_albedo_modis(lat::AbstractVector{<:Real},
                                          lon::AbstractVector{<:Real},
                                          years=MCD43A3_YEARS;
                                          percentile::Real=0.05,
                                          min_samples::Integer=_MODIS_ICE_MIN_SAMPLES,
                                          albedo_range::Tuple{Real,Real}=_ICE_ALBEDO_RANGE,
                                          qa_keep::AbstractVector{<:Integer}=MCD43A3_QA_KEEP,
                                          layers=MCD43A3_ALBEDO_LAYERS,
                                          doy_range::Union{Nothing,Symbol,Tuple{Integer,Integer}}=nothing,
                                          stride::Integer=1,
                                          token::Union{Nothing,AbstractString}=nothing,
                                          cache_path::Union{Nothing,AbstractString}=nothing,
                                          keep_granules::Bool=false,
                                          force_download::Bool=false,
                                          max_concurrent_downloads::Integer=4,
                                          download_timeout::Real=3600,
                                          progress::Bool=true, verbose::Bool=true)
    # ---- validation, all before any network or disk work
    length(lat) == length(lon) || throw(ArgumentError(
        "lat and lon must have the same length; got $(length(lat)) and $(length(lon))"))
    isempty(lat) && throw(ArgumentError("no points supplied"))
    all(l -> -90 <= l <= 90, lat) || throw(ArgumentError(
        "every latitude must be in −90…90"))
    all(l -> -360 <= l <= 360, lon) || throw(ArgumentError(
        "every longitude must be in −360…360"))

    # `vec(collect(...))`, not `collect`: a scalar `years=2019` collects to a *0-dimensional*
    # array, whose `Date.(...)` is also 0-d and cannot be a `Ti` lookup.
    year_list = vec(collect(years))
    isempty(year_list) && throw(ArgumentError("no years requested"))
    bad = setdiff(year_list, MCD43A3_YEARS)
    isempty(bad) || throw(ArgumentError(
        "years outside the MCD43A3 v061 record $(MCD43A3_YEARS): $(join(bad, ", "))"))

    (0 < percentile <= 1) || throw(ArgumentError(
        "percentile must be in (0, 1]; got $(percentile)"))
    min_samples >= 1 || throw(ArgumentError("min_samples must be ≥ 1; got $(min_samples)"))
    albedo_range[1] < albedo_range[2] || throw(ArgumentError(
        "albedo_range must be (low, high) with low < high; got $(albedo_range)"))
    stride >= 1 || throw(ArgumentError("stride must be ≥ 1; got $(stride)"))
    doy_range isa Symbol && doy_range !== :melt_season && throw(ArgumentError(
        "doy_range as a symbol must be :melt_season; got :$(doy_range)"))
    # Resolve `:melt_season` against the points' hemispheres before validating, so the rest
    # of the function only ever sees `nothing` or an explicit window.
    doy_range = _resolve_doy_range(doy_range, lat)
    if !isnothing(doy_range)
        lo, hi = doy_range
        # `lo > hi` is legal and means the window wraps New Year, e.g. the southern melt
        # season (300, 110). Only the bounds themselves are checked.
        all(d -> 1 <= d <= 366, doy_range) || throw(ArgumentError(
            "doy_range days must be in 1…366; got $(doy_range)"))
        lo == hi + 1 && throw(ArgumentError(
            "doy_range $(doy_range) selects the whole year; pass `nothing` instead"))
    end

    layer_list = collect(Symbol.(layers))
    isempty(layer_list) && throw(ArgumentError("no layers requested"))
    unknown = setdiff(layer_list, MCD43A3_LAYERS)
    isempty(unknown) || throw(ArgumentError(
        "not MCD43A3 albedo layers: $(join(unknown, ", ")). Valid: $(join(MCD43A3_LAYERS, ", "))"))
    isempty(qa_keep) && throw(ArgumentError(
        "qa_keep is empty, so every observation would be rejected; pass at least [0]"))

    dates_by_year = Dict(y => _modis_dates(y, doy_range, stride) for y in year_list)
    empty_years = [y for y in year_list if isempty(dates_by_year[y])]
    isempty(empty_years) || throw(ArgumentError("""
        No MCD43A3 dates for year(s) $(join(empty_years, ", ")) under \
        doy_range=$(doy_range), stride=$(stride). The record starts $(_MCD43A3_FIRST_DATE).
        """))

    _modis_assert_hdf4_driver()
    _configure_gdal_http()

    # ---- dedup: this is what makes a point list cheap
    cells, cell_of_point = _modis_dedup_points(lat, lon)
    ncell = length(cells)
    npoint = length(lat)
    cells_key = _modis_cells_key(cells)
    # Cells grouped by tile, so each granule is opened once per date and every cell it
    # covers is sampled from one windowed read.
    cell_of_tile = Dict{Tuple{Int,Int},Vector{Int}}()
    for (i, c) in enumerate(cells)
        push!(get!(cell_of_tile, (c[1], c[2]), Int[]), i)
    end
    tiles = sort!(collect(keys(cell_of_tile)))

    cache = isnothing(cache_path) ? _default_modis_cache() : cache_path
    tok = isnothing(token) ? get_earthdata_token() : String(token)
    quality_of = Dict(l => _mcd43a3_quality_layer(l) for l in layer_list)
    cache_layers = unique(vcat(layer_list, collect(values(quality_of))))

    if verbose
        @info "MCD43A3 bare-ice albedo" points = npoint unique_cells = ncell tiles = length(tiles) years = "$(first(year_list))–$(last(year_list))" layers = layer_list
        total_dates = sum(length(dates_by_year[y]) for y in year_list)
        @info "Download volume estimate" dates = total_dates granules = total_dates * length(tiles) approx_GB = round(total_dates * length(tiles) * 70 / 1024; digits=1)
    end

    # ---- results, one (cell → point) scatter per layer and year
    albedo = Dict(l => fill(NaN32, npoint, length(year_list)) for l in layer_list)
    counts = Dict(l => zeros(Int32, npoint, length(year_list)) for l in layer_list)

    total_steps = sum(length(dates_by_year[y]) for y in year_list)
    t0 = time()
    done = 0

    for (yi, y) in enumerate(year_list)
        dates = dates_by_year[y]
        # `n_expected` is the whole year's date count, so `kmax` does not depend on how the
        # year happens to be traversed.
        accs = Dict(l => _LowPercentileTopK((ncell, 1), length(dates), percentile,
                                            min_samples) for l in layer_list)
        for date in dates
            samples = force_download ? nothing :
                _modis_sample_cache_read(_modis_sample_cache_path(cache, cells_key, date),
                                         cache_layers, ncell)
            if isnothing(samples)
                paths = mcd43a3_granules(date, tiles; token=tok, cache_path=cache,
                                         force_download=force_download, verbose=verbose,
                                         timeout=download_timeout,
                                         max_concurrent_downloads=max_concurrent_downloads)
                samples = _modis_sample_cells(paths, cells, cell_of_tile, layer_list)
                _modis_sample_cache_write(
                    _modis_sample_cache_path(cache, cells_key, date), cache_layers, samples)
                if !keep_granules
                    # NCDatasets/GDAL hold each file open until its handle is collected, so
                    # the collection must precede the removal — same reason the C3S path's
                    # `discard_after_fold` does this.
                    GC.gc()
                    for p in values(paths)
                        rm(p; force=true)
                    end
                end
            end

            for l in layer_list
                # The accumulator is an (ncell, 1) "grid", so the C3S reduction kernels are
                # reused verbatim.
                masked = _valid_albedo(reshape(samples[l], ncell, 1),
                                       reshape(samples[quality_of[l]], ncell, 1), nothing;
                                       keep_values=qa_keep, scale=_MCD43A3_SCALE,
                                       albedo_range=albedo_range)
                _accumulate!(accs[l], masked)
            end

            done += 1
            progress && _report_progress("MCD43A3 $(y)", done, total_steps, t0)
        end

        for l in layer_list
            cell_albedo, cell_counts = _finalize(accs[l])
            # Scatter cell values back over points: points sharing a cell get bit-identical
            # values by construction, not by coincidence.
            for p in 1:npoint
                c = cell_of_point[p]
                albedo[l][p, yi] = cell_albedo[c, 1]
                counts[l][p, yi] = cell_counts[c, 1]
            end
        end
    end

    return _modis_ice_albedo_stack(albedo, counts, layer_list, cells, cell_of_point,
                                   year_list, percentile, min_samples, albedo_range,
                                   qa_keep, doy_range, stride)
end

function compute_glacier_ice_albedo_modis(points::AbstractVector{<:Tuple{<:Real,<:Real}},
                                          years=MCD43A3_YEARS; kwargs...)
    isempty(points) && throw(ArgumentError("no points supplied"))
    return compute_glacier_ice_albedo_modis(first.(points), last.(points), years; kwargs...)
end

"""
    _modis_ice_albedo_layer_names(layer) -> (albedo_key, count_key)

Output layer names for one MCD43A3 albedo layer: `:Albedo_BSA_shortwave` →
`(:albedo_bsa, :n_valid_observations_bsa)`, and `:Albedo_WSA_vis` →
`(:albedo_wsa_vis, :n_valid_observations_wsa_vis)`.

The broadband default drops the `_shortwave` suffix, because that is the layer GEMB
consumes and `:albedo_bsa` reads better than `:albedo_bsa_shortwave`.
"""
function _modis_ice_albedo_layer_names(layer::Symbol)
    m = match(r"^Albedo_(BSA|WSA)_(.+)$", String(layer))
    sky = lowercase(m[1])
    suffix = m[2] == "shortwave" ? "" : "_" * lowercase(m[2])
    return (Symbol("albedo_", sky, suffix),
            Symbol("n_valid_observations_", sky, suffix))
end

# `"h16v02_r0699_c0124"` — the dedup key, in a form that survives a NetCDF/CSV round-trip.
_modis_cell_id(cell::NTuple{4,Int}) =
    string("h", lpad(cell[1], 2, '0'), "v", lpad(cell[2], 2, '0'),
           "_r", lpad(cell[3], 4, '0'), "_c", lpad(cell[4], 4, '0'))

function _modis_ice_albedo_stack(albedo, counts, layer_list, cells, cell_of_point,
                                 year_list, percentile, min_samples, albedo_range,
                                 qa_keep, doy_range, stride)
    npoint = length(cell_of_point)
    point_dim = Dim{:point}(1:npoint)
    time_dim = Ti(Date.(year_list, 1, 1))
    dims2 = (point_dim, time_dim)

    sampled = [_modis_cell_center(cells[cell_of_point[p]]...) for p in 1:npoint]

    data = Pair{Symbol,Any}[]
    for l in layer_list
        akey, ckey = _modis_ice_albedo_layer_names(l)
        push!(data, akey => DimArray(albedo[l], dims2;
            metadata=Dict("units" => "1", "source_layer" => String(l),
                          "long_name" => "glacier bare-ice albedo (darkest $(round(100 * percentile; digits=1))% mean)")))
        push!(data, ckey => DimArray(counts[l], dims2;
            metadata=Dict("units" => "count", "source_layer" => String(l),
                          "long_name" => "valid MCD43A3 retrievals per point-year")))
    end
    push!(data, :latitude => DimArray([s[1] for s in sampled], (point_dim,);
        metadata=Dict("units" => "degrees_north",
                      "long_name" => "latitude of the sampled MODIS cell centre")))
    push!(data, :longitude => DimArray([s[2] for s in sampled], (point_dim,);
        metadata=Dict("units" => "degrees_east",
                      "long_name" => "longitude of the sampled MODIS cell centre")))
    push!(data, :cell_id => DimArray([_modis_cell_id(cells[cell_of_point[p]])
                                      for p in 1:npoint], (point_dim,);
        metadata=Dict("long_name" => "MCD43A3 sinusoidal grid cell (points sharing this share their albedo)")))

    metadata = Dict{String,Any}(
        "source" => "MODIS MCD43A3 v061 (500 m daily BSA/WSA albedo)",
        "n_points" => npoint,
        "n_unique_cells" => length(cells),
        "percentile" => percentile,
        "min_samples" => min_samples,
        "albedo_range" => collect(albedo_range),
        "qa_keep" => collect(qa_keep),
        "doy_range" => isnothing(doy_range) ? "all" : collect(doy_range),
        "stride" => stride,
        "years" => collect(year_list),
        "layers" => String.(layer_list),
    )
    return DimStack(NamedTuple(data); metadata=metadata)
end
