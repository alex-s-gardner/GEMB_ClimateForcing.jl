"""
MODIS/Terra+Aqua **MCD43A3 v061** — daily 500 m black-sky and white-sky albedo, on the
MODIS sinusoidal grid, from NASA's LP DAAC.

Compared with the C3S record this package also reads (`copernicus_albedo.jl`), this product
is longer (2000-02-16 → present vs 2018–2024), denser (daily vs 10-daily), coarser in space
(500 m vs 300 m), and — the reason it is here — distributes **black-sky**
(`Albedo_BSA_*`, directional-hemispherical) and **white-sky** (`Albedo_WSA_*`,
bihemispherical) albedo as *separate* layers, which GEMB's radiation scheme can use
directly.

Access is a plain authenticated HTTPS GET of a whole granule, discovered through CMR
(`earthdata.jl`). There is **no job queue and no server-side subsetting**, so the cost model
is the exact opposite of the C3S path: bandwidth-bound rather than latency-bound, and the
unit of transfer is a ~50–95 MB granule however few pixels are wanted.

!!! warning "GDAL cannot read HDF4 through `/vsicurl/`, and the failure is a false positive"
    The HDF4 library does its own POSIX I/O and cannot use GDAL's virtual filesystems.
    Opening `/vsicurl/…​.hdf` **appears to succeed** — it returns a dataset handle that only
    throws on the first metadata access — so a byte-range read cannot be substituted for the
    download, however tempting the 85 % of unread bytes makes it. Granules are downloaded to
    local disk, full stop. (Julia's `GDAL_jll` does carry the HDF4 driver; Homebrew's
    `gdalinfo` does not, which is why a shell check can disagree with what runs here.)

Grid: sinusoidal, `+R=6371007.181`, tiled `h##v##` with `h` counting east from the
antimeridian and `v` south from the north pole, 2400×2400 pixels at 463.312716528 m. The
tile geometry is closed-form (see [`_modis_tile_origin`](@ref)), so no granule needs opening
to know which cell a point falls in.

Layer values are `Int16` scaled by **0.001**, applied by us: GDAL reports the scale on the
band but a raw `ArchGDAL.read` returns the unscaled integer. Valid range 0–32766, fill
32767, so albedo may legitimately exceed 1.0 and **no physical clamp belongs here** —
`albedo_range` is quality control, not a unit fix.

Quality (`BRDF_Albedo_Band_Mandatory_Quality_*`, `UInt8`, fill 255) is a small-integer
*class*, not a bitfield: `0` = full BRDF inversion, `1` = magnitude inversion, `2`–`7` =
v061 Terra band 5/6 detector-failure cases. Hence the whitelist (`keep_values`) rather than
the C3S reject-mask, and hence no per-pixel uncertainty layer — MCD43A3 has none. Snow flags
and full per-band QA live in a *separate* granule (MCD43A2) and are out of scope.
"""

using Rasters
using ArchGDAL
using Dates

# ------------------------------------------------------------------------------ the grid

# Sinusoidal projection the MODIS land products are gridded on. `R` is the MODIS-specific
# authalic sphere radius, not WGS84's — using 6371000 puts cells ~7 m out.
const _MODIS_SINU_PROJ =
    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
const _MODIS_SPHERE_R = 6371007.181

# Grid origin (upper-left of tile h00v00) and tile span, in projection metres. The span is
# exactly 10° of arc at the equator on this sphere: 2π·R/36.
const _MODIS_UL_X = -20015109.354
const _MODIS_UL_Y = 10007554.677
const _MODIS_TILE_SPAN_M = 1111950.5196666666

# A 500 m land tile is 2400×2400, so the pixel is the tile span / 2400 — not 500 m exactly.
const _MODIS_TILE_PIXELS = 2400
const _MODIS_PIXEL_M = _MODIS_TILE_SPAN_M / _MODIS_TILE_PIXELS

# Tile index bounds (36 columns × 18 rows covers the globe).
const _MODIS_H_MAX = 35
const _MODIS_V_MAX = 17

# -------------------------------------------------------------------------- the product

const MCD43A3_SHORT_NAME = "MCD43A3"
const MCD43A3_VERSION = "061"

# The HDF-EOS grid name inside every MCD43A3 granule.
const _MCD43A3_GRID = "MOD_Grid_BRDF"

# Applied by us, not by GDAL — see the module docstring.
const _MCD43A3_SCALE = 0.001

"""
    MCD43A3_ALBEDO_LAYERS

The two broadband albedo layers reduced by default: black-sky
(`:Albedo_BSA_shortwave`, directional-hemispherical) and white-sky
(`:Albedo_WSA_shortwave`, bihemispherical).

The granule also carries `Band1`–`Band7`, `vis` and `nir` variants of each. Pass them via
the `layers` keyword of [`compute_glacier_ice_albedo_modis`](@ref) if a spectral breakdown
is wanted; every layer named there must appear in [`MCD43A3_LAYERS`](@ref).
"""
const MCD43A3_ALBEDO_LAYERS = (:Albedo_BSA_shortwave, :Albedo_WSA_shortwave)

# The band suffixes each albedo layer (and its quality companion) exists for.
const _MCD43A3_BANDS = ("Band1", "Band2", "Band3", "Band4", "Band5", "Band6", "Band7",
                        "vis", "nir", "shortwave")

"""
    MCD43A3_LAYERS

Every albedo layer an MCD43A3 granule contains: `Albedo_{BSA,WSA}_{Band1..7,vis,nir,shortwave}`,
twenty in total. Used to validate the `layers` keyword before anything is downloaded.
"""
const MCD43A3_LAYERS =
    Tuple(Symbol("Albedo_$(sky)_$(band)") for sky in ("BSA", "WSA")
          for band in _MCD43A3_BANDS)

"""
    MCD43A3_YEARS

Calendar years with usable MCD43A3 v061 coverage. The record starts 2000-02-16 (so 2000 is
partial — its `min_samples` should be relaxed accordingly) and is ongoing; the upper bound
here is conservative and only gates the `years` argument.
"""
const MCD43A3_YEARS = 2000:2025

# First day with data, for validating a `years`/`doy_range` combination that can return
# nothing.
const _MCD43A3_FIRST_DATE = Date(2000, 2, 16)

"""
    _mcd43a3_quality_layer(layer) -> Symbol

The mandatory-quality layer companion to an albedo layer. Black-sky and white-sky share one
quality band per wavelength, so both `Albedo_BSA_shortwave` and `Albedo_WSA_shortwave` map
to `BRDF_Albedo_Band_Mandatory_Quality_shortwave`.
"""
function _mcd43a3_quality_layer(layer::Symbol)
    s = String(layer)
    m = match(r"^Albedo_(?:BSA|WSA)_(.+)$", s)
    isnothing(m) && throw(ArgumentError(
        "not an MCD43A3 albedo layer: $(layer). Valid: $(join(MCD43A3_LAYERS, ", "))"))
    return Symbol("BRDF_Albedo_Band_Mandatory_Quality_$(m[1])")
end

"""
    _modis_subdataset(path, layer) -> String

GDAL subdataset spec for one layer of a local MCD43A3 granule:

    HDF4_EOS:EOS_GRID:"<path>":MOD_Grid_BRDF:<layer>

`path` must be a **local** file. Prefixing it with `/vsicurl/` yields a handle that opens
without error and then fails on first use — see the module docstring — so this function
refuses a remote path rather than letting that surface later as an unrelated GDAL message.
"""
function _modis_subdataset(path::AbstractString, layer::Union{Symbol,AbstractString})
    (startswith(path, "/vsi") || occursin("://", path)) && throw(ArgumentError("""
        MCD43A3 granules must be read from local disk, not $(path).

        GDAL's HDF4 driver cannot read through /vsicurl/ (HDF4 does its own POSIX I/O), and
        the failure is a *false positive*: the open succeeds and only the first metadata
        access throws. Download the granule first.
        """))
    return "HDF4_EOS:EOS_GRID:\"$(path)\":$(_MCD43A3_GRID):$(layer)"
end

"""
    _modis_assert_hdf4_driver()

Error with actionable advice if the running GDAL has no HDF4 driver, before any granule is
downloaded. `ArchGDAL`'s `GDAL_jll` normally provides it; a system GDAL often does not.
"""
function _modis_assert_hdf4_driver()
    ok = try
        ArchGDAL.getdriver("HDF4") !== nothing
    catch
        false
    end
    ok || error("""
    This GDAL build has no HDF4 driver, so MCD43A3 granules cannot be read.

    ArchGDAL's bundled GDAL_jll normally includes it. If a system GDAL is being picked up
    instead (e.g. via a JLL preference or LD_LIBRARY_PATH), remove that override, or
    reinstall a GDAL built with `--with-hdf4`.
    """)
    return nothing
end

# ----------------------------------------------------------------------- tile arithmetic

"""
    _modis_tile_origin(h, v) -> (ulx, uly)

Upper-left corner of MODIS sinusoidal tile `(h, v)` in projection metres:

    ulx = _MODIS_UL_X + h · _MODIS_TILE_SPAN_M
    uly = _MODIS_UL_Y − v · _MODIS_TILE_SPAN_M

Pure arithmetic — no granule is opened, which is what lets a point list be mapped to tiles
and cells before anything is downloaded. Verified against a real granule's GDAL
geotransform to under 2 mm.

```jldoctest
julia> GEMB_ClimateForcing._modis_tile_origin(0, 0)
(-2.0015109354e7, 1.0007554677e7)
```
"""
function _modis_tile_origin(h::Integer, v::Integer)
    _modis_check_tile(h, v)
    return (_modis_tile_edge_x(h), _modis_tile_edge_y(v))
end

# Grid-line positions, indexed by tile *edge* rather than tile. Every edge is computed from
# the one expression `ULX + h·T`, never as `<previous edge> + T`, so a tile's `xmax` is
# bit-identical to its neighbour's `xmin` — `ulx + T` and `ULX + (h+1)·T` round differently
# in Float64, and a gap of one ulp between tiles would put a point on the boundary into
# neither of them. Deliberately unguarded: they are called with `h+1`/`v+1` at the last
# tile, where the *edge* exists even though the tile does not.
_modis_tile_edge_x(h::Integer) = _MODIS_UL_X + h * _MODIS_TILE_SPAN_M
_modis_tile_edge_y(v::Integer) = _MODIS_UL_Y - v * _MODIS_TILE_SPAN_M

"""
    _modis_tile_bounds(h, v) -> (xmin, ymin, xmax, ymax)

Projection-metre bounding box of tile `(h, v)`. Adjacent tiles share an edge **exactly**, so
`_modis_tile_bounds(h, v)[3] === _modis_tile_bounds(h + 1, v)[1]` — see
[`_modis_tile_edge_x`](@ref) for why that needs care.
"""
function _modis_tile_bounds(h::Integer, v::Integer)
    _modis_check_tile(h, v)
    return (_modis_tile_edge_x(h), _modis_tile_edge_y(v + 1),
            _modis_tile_edge_x(h + 1), _modis_tile_edge_y(v))
end

function _modis_check_tile(h::Integer, v::Integer)
    (0 <= h <= _MODIS_H_MAX && 0 <= v <= _MODIS_V_MAX) || throw(ArgumentError(
        "MODIS tile out of range: h=$(h), v=$(v) (expected 0:$(_MODIS_H_MAX), 0:$(_MODIS_V_MAX))"))
    return nothing
end

"""
    _modis_lonlat_to_sinu(lon, lat) -> (x, y)

Forward MODIS sinusoidal projection, in metres:

    x = R · λ · cos φ
    y = R · φ

with λ, φ in radians. Closed form on purpose — validated against PROJ (through ArchGDAL) to
0.0 m at several high- and low-latitude points, so this needs no `Proj` dependency and no
per-point reprojection call.
"""
function _modis_lonlat_to_sinu(lon::Real, lat::Real)
    φ = deg2rad(lat)
    λ = deg2rad(_wrap_longitude(lon))
    return (_MODIS_SPHERE_R * λ * cos(φ), _MODIS_SPHERE_R * φ)
end

"""
    _modis_cell(lat, lon) -> (h, v, row, col)

The MCD43A3 500 m grid cell containing a point: tile indices plus **1-based** row (from the
tile's north edge) and column (from its west edge).

This is the deduplication key. Two points in the same cell necessarily share every
retrieval, so their albedo is derived once — see
[`compute_glacier_ice_albedo_modis`](@ref).

```jldoctest
julia> GEMB_ClimateForcing._modis_cell(0.0, 0.0)
(18, 9, 1, 1)
```
"""
function _modis_cell(lat::Real, lon::Real)
    -90 <= lat <= 90 ||
        throw(ArgumentError("latitude must be in −90…90, got $(lat)"))
    x, y = _modis_lonlat_to_sinu(lon, lat)

    h = floor(Int, (x - _MODIS_UL_X) / _MODIS_TILE_SPAN_M)
    v = floor(Int, (_MODIS_UL_Y - y) / _MODIS_TILE_SPAN_M)
    # A point exactly on the grid's east or south edge lands one tile past the end.
    h = clamp(h, 0, _MODIS_H_MAX)
    v = clamp(v, 0, _MODIS_V_MAX)

    ulx, uly = _modis_tile_origin(h, v)
    col = clamp(floor(Int, (x - ulx) / _MODIS_PIXEL_M) + 1, 1, _MODIS_TILE_PIXELS)
    row = clamp(floor(Int, (uly - y) / _MODIS_PIXEL_M) + 1, 1, _MODIS_TILE_PIXELS)
    return (h, v, row, col)
end

"""
    _modis_cell_center(h, v, row, col) -> (lat, lon)

Geographic centre of a grid cell, by inverting the sinusoidal projection. Reported back to
the caller so a point's *sampled* location is auditable against the location it asked for.
"""
function _modis_cell_center(h::Integer, v::Integer, row::Integer, col::Integer)
    ulx, uly = _modis_tile_origin(h, v)
    x = ulx + (col - 0.5) * _MODIS_PIXEL_M
    y = uly - (row - 0.5) * _MODIS_PIXEL_M
    φ = y / _MODIS_SPHERE_R
    # cos φ → 0 at the poles, where longitude is degenerate; guard the division.
    c = cos(φ)
    λ = abs(c) < 1e-12 ? 0.0 : x / (_MODIS_SPHERE_R * c)
    return (rad2deg(φ), rad2deg(λ))
end

# ---------------------------------------------------------------------- granule identity

"""
    _modis_granule_date(id) -> Date

Nominal acquisition date of a granule, from the `A<YYYYDDD>` field of its
`producer_granule_id`.

!!! warning "This, never `time_start`"
    A granule's CMR `time_start`/`time_end` span the product's full **16-day retrieval
    window**, so a one-day `temporal` query returns granules for 16 distinct nominal dates
    (verified live: `cmr-hits: 32` for one day and two tiles). Only the granule id carries
    the nominal date.

```jldoctest
julia> GEMB_ClimateForcing._modis_granule_date("MCD43A3.A2019182.h17v02.061.2019191033013.hdf")
2019-07-01
```
"""
function _modis_granule_date(id::AbstractString)
    m = match(r"\.A(\d{4})(\d{3})\.", id)
    isnothing(m) && throw(ArgumentError(
        "no A<YYYYDDD> nominal-date field in MODIS granule id: $(id)"))
    return Date(parse(Int, m[1])) + Day(parse(Int, m[2]) - 1)
end

"""
    _modis_granule_tile(id) -> (h, v)

Tile indices of a granule, from the `h##v##` field of its `producer_granule_id`.
"""
function _modis_granule_tile(id::AbstractString)
    m = match(r"\.h(\d{2})v(\d{2})\.", id)
    isnothing(m) && throw(ArgumentError(
        "no h##v## tile field in MODIS granule id: $(id)"))
    h, v = parse(Int, m[1]), parse(Int, m[2])
    _modis_check_tile(h, v)
    return (h, v)
end

"""
    _modis_tile_bbox(h, v) -> String

CMR `bounding_box` ("W,S,E,N") covering tile `(h, v)`, from the geographic corners of its
projected extent.

The sinusoidal tile is not a lat/lon rectangle, so the box is the *bounding* box of its
corners and mid-edge points — deliberately generous. Over-selecting is harmless (granules
are filtered by their id's tile field afterwards); under-selecting would silently drop a
tile.
"""
function _modis_tile_bbox(h::Integer, v::Integer)
    xmin, ymin, xmax, ymax = _modis_tile_bounds(h, v)
    lats = Float64[]
    lons = Float64[]
    # Sample the edges rather than just the corners: the east/west edges bow with cos φ, so
    # the extreme longitude sits at whichever edge is closest to the equator.
    for y in range(ymin, ymax; length=5), x in (xmin, xmax)
        φ = y / _MODIS_SPHERE_R
        abs(φ) > π / 2 && continue
        c = cos(φ)
        λ = abs(c) < 1e-12 ? 0.0 : x / (_MODIS_SPHERE_R * c)
        push!(lats, rad2deg(φ))
        push!(lons, clamp(rad2deg(λ), -180.0, 180.0))
    end
    isempty(lats) && throw(ArgumentError("tile h=$(h), v=$(v) has no valid geographic extent"))
    return string(minimum(lons), ",", minimum(lats), ",", maximum(lons), ",", maximum(lats))
end

# ----------------------------------------------------------------------------- retrieval

_default_modis_cache() =
    joinpath(tempdir(), "GEMB_ClimateForcing", "$(MCD43A3_SHORT_NAME).$(MCD43A3_VERSION)")

"""
    mcd43a3_granules(date, tiles; token=nothing, cache_path=nothing, force_download=false,
                     verbose=true, timeout=Inf, deadline=Inf,
                     max_concurrent_downloads=4) -> Dict{Tuple{Int,Int},String}

Local paths to the MCD43A3 granules for one **nominal** `date` and the requested
`tiles::AbstractVector{Tuple{Int,Int}}`, downloading whatever is not already cached.

Tiles with no granule for that date are simply absent from the result — a gap in the record
is normal and is not an error.

CMR is queried once per tile (the tile's bounding box), and results are filtered on the
granule id's nominal date *and* tile, because `temporal` alone returns 16 dates. Downloads
run concurrently through [`_run_concurrent_jobs`](@ref) with no submission stagger: unlike
the CDS path there is no queue limiter to trip, only bandwidth to share.
"""
function mcd43a3_granules(date::Date, tiles::AbstractVector{<:Tuple{Integer,Integer}};
                          token::Union{Nothing,AbstractString}=nothing,
                          cache_path::Union{Nothing,AbstractString}=nothing,
                          force_download::Bool=false, verbose::Bool=true,
                          timeout::Real=Inf, deadline::Float64=Inf,
                          max_concurrent_downloads::Integer=4)
    isempty(tiles) && return Dict{Tuple{Int,Int},String}()
    cache = isnothing(cache_path) ? _default_modis_cache() : cache_path
    mkpath(cache)

    tile_list = collect(Tuple{Int,Int}.(tiles))
    # (id, url, bytes) per tile that has a granule, resolved from CMR unless the file is
    # already on disk (in which case no query is needed at all).
    wanted = Dict{Tuple{Int,Int},NamedTuple{(:path, :url, :bytes),Tuple{String,String,Int}}}()
    found = Dict{Tuple{Int,Int},String}()

    for tile in tile_list
        cached = _modis_cached_granule(cache, date, tile)
        if !force_download && !isnothing(cached)
            found[tile] = cached
            continue
        end
        temporal = _modis_cmr_temporal(date)
        granules = _cmr_granules(; short_name=MCD43A3_SHORT_NAME, version=MCD43A3_VERSION,
                                 temporal=temporal, bounding_box=_modis_tile_bbox(tile...),
                                 verbose=verbose, deadline=deadline)
        hit = nothing
        for g in granules
            # BOTH filters are required: `temporal` returns 16 nominal dates, and the tile
            # bounding box overlaps neighbouring tiles.
            (_modis_granule_date(g.id) == date && _modis_granule_tile(g.id) == tile) ||
                continue
            hit = g
            break
        end
        isnothing(hit) && continue
        wanted[tile] = (; path=joinpath(cache, _modis_granule_filename(hit.id)),
                        url=hit.url, bytes=hit.bytes)
    end

    if !isempty(wanted)
        pending = collect(wanted)
        tok = isnothing(token) ? get_earthdata_token() : String(token)
        results = Vector{Union{Nothing,Pair{Tuple{Int,Int},String}}}(nothing,
                                                                    length(pending))
        _run_concurrent_jobs(length(pending), max_concurrent_downloads, 0;
                             verbose=false) do i, _gate
            tile, spec = pending[i]
            path = spec.path
            if !force_download && isfile(path) && (spec.bytes == 0 ||
                filesize(path) >= spec.bytes * (1 - _EARTHDATA_SIZE_TOLERANCE))
                results[i] = tile => path
                return nothing
            end
            verbose && @info "Downloading MCD43A3 granule" date tile size_mb = round(spec.bytes / 1024^2; digits=1)
            _earthdata_download(spec.url, path; token=tok, expected_bytes=spec.bytes,
                                verbose=verbose, timeout=timeout, deadline=deadline)
            results[i] = tile => path
            return nothing
        end
        for r in results
            isnothing(r) || (found[r.first] = r.second)
        end
    end
    return found
end

# The archive filename, from a producer_granule_id that may or may not already carry the
# extension.
_modis_granule_filename(id::AbstractString) =
    endswith(id, ".hdf") ? String(id) : string(id, ".hdf")

# A granule's production-timestamp field is unpredictable, so a cached file is found by
# globbing on the fields that *are* known: product, nominal date and tile.
function _modis_cached_granule(cache::AbstractString, date::Date, tile::Tuple{Integer,Integer})
    isdir(cache) || return nothing
    stem = string(MCD43A3_SHORT_NAME, ".A", year(date), lpad(dayofyear(date), 3, '0'),
                  ".h", lpad(tile[1], 2, '0'), "v", lpad(tile[2], 2, '0'), ".")
    for f in readdir(cache)
        startswith(f, stem) && endswith(f, ".hdf") && return joinpath(cache, f)
    end
    return nothing
end

"""
    _modis_cmr_temporal(date) -> String

CMR `temporal` range for a nominal `date`.

A single day is enough — a granule's window *contains* its nominal date — and the extra 15
nominal dates it drags in are filtered out by [`_modis_granule_date`](@ref). Narrowing this
further is impossible; widening it only costs more filtering.
"""
_modis_cmr_temporal(date::Date) =
    string(date, "T00:00:00Z,", date, "T23:59:59Z")

"""
    _modis_read_window(path, layer, rows, cols) -> Matrix

Read the `(rows, cols)` window of one granule layer, as a `Matrix` indexed
`[row_offset, col_offset]` — i.e. row-major in the caller's terms, matching
[`_modis_cell`](@ref)'s `(row, col)`.

`ArchGDAL.read(ds, band, rows, cols)` takes the ranges in that order but returns an
**x-major** array — `length(cols) × length(rows)`, matching GDAL's own `[x, y]` storage.
The `permutedims` here is what makes the result agree with `(row, col)` indexing, and a test
pins the convention on a real granule, because getting it backwards silently samples the
wrong pixels rather than erroring.
"""
function _modis_read_window(path::AbstractString, layer::Union{Symbol,AbstractString},
                            rows::AbstractUnitRange, cols::AbstractUnitRange)
    spec = _modis_subdataset(path, layer)
    return ArchGDAL.read(spec) do ds
        permutedims(ArchGDAL.read(ds, 1, rows, cols), (2, 1))
    end
end
