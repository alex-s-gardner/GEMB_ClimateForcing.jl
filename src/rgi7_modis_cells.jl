"""
RGI 7.0 glacier outlines, rasterized onto the **native MCD43A3 500 m sinusoidal grid**.

Two halves live here. The *burn kernel* (`_rgi7_*`) maps glacier polygons to
`(h, v, row, col)` cells and is driven by `data/make_rgi7_modis_cells.jl`, which is run by
hand when RGI is updated. The *loader* reads the resulting vendored tables
(`data/rgi7_modis_cells.csv.gz`, `data/rgi7_glaciers.csv.gz`) and hands the cell list to
[`compute_glacier_ice_albedo_modis`](@ref).

Rasterizing *into* the product's own grid, rather than warping the product onto a lat/lon
grid, is the whole point: an `:average` resample of an individual day mixes neighbouring
pixels' albedo *before* the darkest-percentile reduction, biasing the very tail being
measured. See the note in `glacier_ice_albedo_modis.jl`.

!!! note "Two rules that shape the artifact, both deliberate"
    **Cells are assigned centre-in-polygon, and a contested cell goes to the *smaller*
    glacier.** Both are the opposite of the obvious choice, for reasons recorded at
    [`_rgi7_burn_tile!`](@ref) and [`_rgi7_tile_cells`](@ref):

    - `ALL_TOUCHED` is rejected on *physical* grounds, not to save space. A perimeter cell is
      majority off-glacier — rock, moraine, water, shadow — all **darker** than ice. So it
      does not merely add noisy cells, it adds cells biased low in exactly the quantity being
      measured, and a darkest-5 % reduction amplifies that bias rather than averaging it away.
    - Smallest-area-wins, because a 500-cell glacier that concedes one cell loses 0.2 %, while
      a one-cell glacier that concedes its only cell **disappears from the product**.

!!! warning "Glaciers smaller than a cell get one forced cell, and it is biased"
    83.6 % of RGI 7.0 is under 1 km² against a 0.2146 km² cell, so centre-in-polygon alone
    drops a large fraction of the inventory *by count*. Such glaciers are represented by the
    single cell containing their representative point, flagged `forced`. That cell is majority
    non-glacier, so its albedo carries the same dark bias `ALL_TOUCHED` would have — hence the
    flag travels all the way into the output product. Filter on it before averaging.

    The forced cell is **added, never substituted**. Letting it displace the cell's existing
    owner under smallest-wins was tried and is wrong: on Iceland it displaced 13 owners and
    left 8 glaciers (1.4 %) with nothing, which is the very erasure smallest-wins exists to
    prevent. So a cell may carry more than one glacier, and only ever when a glacier would
    otherwise have none — two sub-cell glaciers in one 463 m pixel really do share that
    pixel's albedo. See [`RGI7ModisCells`](@ref) for what that means when sampling.
"""

using Rasters      # rasterize!, Raster
using ArchGDAL     # vector side only here: geometry walking, pointonsurface, envelope
# Sampled/Regular/Intervals/Center and the order types are NOT re-exported by Rasters or by
# DimensionalData itself — they live in this submodule. The burn target has to be built from
# them explicitly rather than from a range, so that its coordinates are bit-identical to
# `_modis_cell_center`'s; see `_rgi7_tile_dims`.
using DimensionalData.Lookups

# Sentinel for "no glacier in this cell". `typemax` rather than 0 so the `minimum` reducer
# needs no special-casing: `min(typemax, rank) == rank` for any real rank.
const _RGI7_NO_GLACIER = typemax(Int32)

# Area of one 500 m cell, km². The pixel is the tile span / 2400, not 500 m exactly, so this
# is 0.21464 km² — used only for the generator's area sanity check.
const _RGI7_CELL_AREA_KM2 = (_MODIS_PIXEL_M / 1000)^2

# --------------------------------------------------------------------------- the burn grid

"""
    _rgi7_tile_dims(h, v) -> (X, Y)

Dimensions of a MODIS tile as a Rasters burn target: 2400×2400 **cell centres** in
projection metres, `X` ascending east and `Y` descending south (row 1 is the tile's north
edge, as in the granule).

The centres are built by the *same expression* [`_modis_cell_center`](@ref) uses, and the
origin comes from [`_modis_tile_origin`](@ref) — i.e. through `_modis_tile_edge_x/_y`, never
`ulx + T`. That is what makes the round trip in [`_rgi7_tile_cells`](@ref) hold by
construction rather than by two implementations happening to agree: Rasters evaluates its
inside test at `first(lookup) + (i-1)·step`, which for these lookups *is* the cell centre.

An explicit `Vector` rather than a range, so every value is bit-identical to
`_modis_cell_center`'s arithmetic; `Regular` so `step` is exactly `±_MODIS_PIXEL_M`; and
`Intervals(Center())` so the `Touches` bounds Rasters derives are cell *edges* and its
`maybeshiftlocus(Center(), ·)` is a no-op.
"""
function _rgi7_tile_dims(h::Integer, v::Integer)
    ulx, uly = _modis_tile_origin(h, v)
    xs = [ulx + (c - 0.5) * _MODIS_PIXEL_M for c in 1:_MODIS_TILE_PIXELS]
    ys = [uly - (r - 0.5) * _MODIS_PIXEL_M for r in 1:_MODIS_TILE_PIXELS]
    return (X(Sampled(xs; order=ForwardOrdered(), span=Regular(_MODIS_PIXEL_M),
                      sampling=Intervals(Center()))),
            Y(Sampled(ys; order=ReverseOrdered(), span=Regular(-_MODIS_PIXEL_M),
                      sampling=Intervals(Center()))))
end

"""
    _rgi7_new_tile_grid() -> Matrix{Int32}

A fresh 2400×2400 burn grid, filled with [`_RGI7_NO_GLACIER`](@ref).

Indexed `[col, row]`, because [`_rgi7_tile_dims`](@ref) puts `X` first. 23 MB each, and the
generator holds one per touched tile (103 of them, ~2.4 GB) so that the global `minimum`
reduction can stay order-independent while regions are streamed one at a time.
"""
_rgi7_new_tile_grid() = fill(_RGI7_NO_GLACIER, _MODIS_TILE_PIXELS, _MODIS_TILE_PIXELS)

# ----------------------------------------------------------------- WGS84 → sinusoidal

"""
    _rgi7_sinu_ring(ring, rgi_id) -> Vector{Tuple{Float64,Float64}}

One polygon ring reprojected to MODIS sinusoidal metres with
[`_modis_lonlat_to_sinu`](@ref).

The closed form is used rather than PROJ (through `ArchGDAL.createcoordtrans`) on purpose:
it is the *identical* function [`_modis_cell`](@ref) uses, so the burn and the cell lookup
cannot disagree. Using PROJ here and the closed form there would introduce a second
implementation that only probably agrees. `test_rgi7_modis_cells.jl` pins the closed form
against PROJ to well under a millimetre, which is the right place for that comparison.

Throws on a ring straddling the antimeridian. `x = R·λ·cos φ` with `λ ∈ (−180°, 180°]` sends
179.9° and −179.9° to opposite ends of the grid, so such a ring silently becomes a
planet-spanning polygon that would burn millions of spurious cells. RGI 7.0 is not believed
to contain one, but "not believed to" is not a guard.
"""
function _rgi7_sinu_ring(ring, rgi_id::AbstractString)
    n = ArchGDAL.ngeom(ring)
    n >= 4 || throw(ArgumentError(
        "$(rgi_id): ring has $(n) vertices, need at least 4 for a closed polygon"))

    lons = Vector{Float64}(undef, n)
    lats = Vector{Float64}(undef, n)
    for j in 1:n
        lons[j] = ArchGDAL.getx(ring, j - 1)
        lats[j] = ArchGDAL.gety(ring, j - 1)
    end

    lo, hi = extrema(lons)
    hi - lo > 180 && error(
        "$(rgi_id): outline spans $(round(hi - lo; digits=1))° of longitude " *
        "($(round(lo; digits=3))…$(round(hi; digits=3))), which the sinusoidal " *
        "projection cannot represent as one ring. Split the outline at the antimeridian.")

    return [_modis_lonlat_to_sinu(lons[j], lats[j]) for j in 1:n]
end

"""
    _rgi7_sinu_geometry(geom, rgi_id) -> ArchGDAL.IGeometry

An RGI 7.0 outline (Polygon or MultiPolygon, EPSG:4326) rebuilt in MODIS sinusoidal metres.

Rebuilt rather than transformed in place because `ArchGDAL.reproject` mutates and would pull
PROJ into the hot path; see [`_rgi7_sinu_ring`](@ref) for why the closed form is preferred.
Interior rings are carried through in order, so nunataks stay holes — Rasters' even-odd
crossing test then leaves them unburned, which `test_rgi7_modis_cells.jl` pins against a real
outline rather than trusting.
"""
function _rgi7_sinu_geometry(geom, rgi_id::AbstractString)
    gt = ArchGDAL.getgeomtype(geom)
    if gt == ArchGDAL.wkbPolygon || gt == ArchGDAL.wkbPolygon25D
        rings = [_rgi7_sinu_ring(ArchGDAL.getgeom(geom, i - 1), rgi_id)
                 for i in 1:ArchGDAL.ngeom(geom)]
        return ArchGDAL.createpolygon(rings)
    elseif gt == ArchGDAL.wkbMultiPolygon || gt == ArchGDAL.wkbMultiPolygon25D
        out = ArchGDAL.createmultipolygon()
        for i in 1:ArchGDAL.ngeom(geom)
            ArchGDAL.addgeom!(out, _rgi7_sinu_geometry(ArchGDAL.getgeom(geom, i - 1), rgi_id))
        end
        return out
    end
    throw(ArgumentError("$(rgi_id): unsupported geometry type $(gt); " *
                        "expected Polygon or MultiPolygon"))
end

"""
    _rgi7_candidate_tiles(sinu_geom) -> Vector{Tuple{Int,Int}}

MODIS tiles whose extent the sinusoidal envelope of a glacier overlaps.

The envelope is taken **after** reprojection, in metres. Taking it in WGS84 and reprojecting
the corners would be wrong: sinusoidal does not preserve axis-aligned boxes, so a lon/lat box
maps to a curved quadrilateral whose metre-space envelope is larger *and* differently placed.

No ±1 padding: the envelope bounds the geometry exactly, and Rasters' `Touches` selector
clamps a selector that runs past a tile edge, so a glacier straddling a boundary is handled
by each tile it genuinely reaches and by no others.
"""
function _rgi7_candidate_tiles(sinu_geom)
    env = ArchGDAL.envelope(sinu_geom)
    h0 = clamp(floor(Int, (env.MinX - _MODIS_UL_X) / _MODIS_TILE_SPAN_M), 0, _MODIS_H_MAX)
    h1 = clamp(floor(Int, (env.MaxX - _MODIS_UL_X) / _MODIS_TILE_SPAN_M), 0, _MODIS_H_MAX)
    # y is inverted: MaxY is the *north* edge, which is the *lowest* v.
    v0 = clamp(floor(Int, (_MODIS_UL_Y - env.MaxY) / _MODIS_TILE_SPAN_M), 0, _MODIS_V_MAX)
    v1 = clamp(floor(Int, (_MODIS_UL_Y - env.MinY) / _MODIS_TILE_SPAN_M), 0, _MODIS_V_MAX)
    return [(h, v) for h in h0:h1 for v in v0:v1]
end

"""
    _rgi7_representative_cell(geom, rgi_id) -> (h, v, row, col)

The cell containing a glacier's representative point, for the forced-cell rule.

`ArchGDAL.pointonsurface` rather than the centroid, because a centroid is only guaranteed
inside a *convex* polygon: a crescent-shaped or multi-lobed glacier — of which RGI has many —
has its centroid off the ice, which would assign the forced cell to a neighbouring valley.
PointOnSurface is guaranteed interior.

Takes the **original WGS84** geometry, since [`_modis_cell`](@ref) wants lon/lat.
"""
function _rgi7_representative_cell(geom, rgi_id::AbstractString)
    pt = ArchGDAL.pointonsurface(geom)
    lon, lat = ArchGDAL.getx(pt, 0), ArchGDAL.gety(pt, 0)
    isfinite(lon) && isfinite(lat) || error(
        "$(rgi_id): representative point is not finite ($(lon), $(lat))")
    return _modis_cell(lat, lon)
end

# -------------------------------------------------------------------------------- the burn

"""
    _rgi7_burn_tile!(grid, h, v, sinu_geoms, ranks) -> grid

Burn glacier **area ranks** into one tile's grid, keeping the lowest rank per cell.

`ranks` is each geometry's position in a global ascending-area ordering, so `minimum`
resolves a contested cell to the **smallest** glacier. Two properties follow, and both are
the reason for this design rather than a `fill`-last or largest-wins rule:

1. **A one-cell glacier cannot be erased.** Conceding its only cell would delete it from the
   product entirely, whereas a 500-cell glacier concedes one and loses 0.2 %.
2. **The artifact is order-independent.** `min` is commutative, so the vendored file is
   bit-identical however the 19 regional shapefiles are streamed, whether one region is
   re-run, and whether the burn is threaded. A last-wins rule would make a *vendored data
   file* a function of directory iteration order.

`boundary=:center` is the centre-in-polygon rule; see this file's header for why
`ALL_TOUCHED` is rejected on physical rather than budgetary grounds. `verbose`/`progress` are
off because Rasters warns once per non-burning geometry and there are 274,531 of them.

If the `reducer`/`init` plumbing ever regresses, the fallback is `Rasters.boolmask!` per
geometry plus an explicit `min`-merge over `findall` — noted so the next reader need not
rediscover it.
"""
function _rgi7_burn_tile!(grid::AbstractMatrix{Int32}, h::Integer, v::Integer,
                          sinu_geoms, ranks::AbstractVector{Int32})
    isempty(sinu_geoms) && return grid
    length(sinu_geoms) == length(ranks) || throw(ArgumentError(
        "geometry/rank length mismatch: $(length(sinu_geoms)) vs $(length(ranks))"))

    ras = Raster(grid, _rgi7_tile_dims(h, v); missingval=_RGI7_NO_GLACIER)
    Rasters.rasterize!(minimum, ras, sinu_geoms;
                       fill=ranks, init=_RGI7_NO_GLACIER, boundary=:center,
                       verbose=false, progress=false, threaded=false)
    return grid
end

"""
    _rgi7_tile_cells(grid, h, v) -> Vector{Tuple{Int,Int,Int,Int,Int}}

Extract the burned cells of one tile as `(h, v, row, col, rank)`, in canonical
`(row, col)` order.

`grid` is indexed `[col, row]` — `X` is the first dimension of
[`_rgi7_tile_dims`](@ref) — and getting that backwards would silently transpose every
glacier, so `test_rgi7_modis_cells.jl` pins it.

Every emitted cell is checked to satisfy

    _modis_cell(_modis_cell_center(h, v, row, col)...) === (h, v, row, col)

and the generator errors out if one does not. This is not defensive padding: the public
albedo API takes lon/lat, so the vendored cells are fed back through as *centres*, and the
identity is the assumption the whole pipeline rests on. It can genuinely fail near the poles,
where `cos φ → 0` makes the recovered longitude exceed ±180° and `_wrap_longitude` folds it —
those are the sinusoidal grid's off-Earth corner cells, which real glacier geometry should
never reach.
"""
function _rgi7_tile_cells(grid::AbstractMatrix{Int32}, h::Integer, v::Integer)
    out = NTuple{5,Int}[]
    for row in 1:_MODIS_TILE_PIXELS, col in 1:_MODIS_TILE_PIXELS
        rank = grid[col, row]
        rank == _RGI7_NO_GLACIER && continue
        back = _modis_cell(_modis_cell_center(h, v, row, col)...)
        back === (Int(h), Int(v), row, col) || error(
            "cell round trip failed at h=$(h) v=$(v) row=$(row) col=$(col): " *
            "centre maps back to $(back). This is an off-Earth corner of the sinusoidal " *
            "grid and must not have been burned.")
        push!(out, (Int(h), Int(v), row, col, Int(rank)))
    end
    return out
end

# =========================================================================== the vendored
# tables, and the loader the albedo driver actually calls

const _RGI7_CELLS_CSV    = joinpath(@__DIR__, "..", "data", "rgi7_modis_cells.csv.gz")
const _RGI7_GLACIERS_CSV = joinpath(@__DIR__, "..", "data", "rgi7_glaciers.csv.gz")

# Row counts of the shipped tables. Hints for `sizehint!` and a `@warn` tripwire only, so a
# regenerated file of any size still loads.
const _RGI7_N_GLACIERS = 274_531
const _RGI7_AREA_KM2 = 706_744.0

"""
    RGI7ModisCells

RGI 7.0 glaciers rasterized onto the MCD43A3 500 m grid, as read from the vendored tables.

Two groups of parallel vectors. **Per row**, in canonical `(h, v, row, col, glacier)` order:
`h`, `v`, `row`, `col`, and `glacier` (a 1-based index into the second group). **Per
glacier**, sorted by `rgi_id`: `rgi_id`, `area_km2`, `n_cells`, and `forced`.

!!! note "Rows are unique as (cell, glacier) pairs, not as cells"
    A cell carries more than one glacier in exactly one situation: when a glacier smaller than
    a 463 m pixel would otherwise have no cell at all (see `forced` below). Two sub-cell
    glaciers in one pixel genuinely share that pixel's albedo, so both are recorded rather
    than one being deleted. Every glacier therefore has `n_cells >= 1`.

    The practical consequence: **the albedo point list is
    [`rgi7_modis_unique_cells`](@ref), not the raw columns**, or `length` would over-count the
    pixels to download. Sharing is rare — measured at well under 1 % of rows.

Longitude and latitude are deliberately *not* stored: they are recoverable exactly from
`(h, v, row, col)` through [`_modis_cell_center`](@ref), and storing them would triple the
file. Use [`rgi7_modis_cell_points`](@ref) to materialise them.

!!! warning "`forced` marks a biased cell"
    `forced[g]` means glacier `g` is smaller than one 0.2146 km² cell and contains no cell
    centre, so it is represented by the single cell holding its representative point. That
    cell is majority *not* glacier — rock, moraine or water, all darker than ice — so its
    bare-ice albedo is biased low. This affects roughly a third of RGI 7.0 by count but only
    ~2 % by area. **Filter on it before computing any aggregate.**
"""
struct RGI7ModisCells
    # per cell, sorted by (h, v, row, col)
    h::Vector{UInt8}
    v::Vector{UInt8}
    row::Vector{UInt16}
    col::Vector{UInt16}
    glacier::Vector{Int32}
    # per glacier, sorted by rgi_id
    rgi_id::Vector{String}
    area_km2::Vector{Float32}
    n_cells::Vector{Int32}
    forced::Vector{Bool}
end

"""
    length(cells::RGI7ModisCells) -> Int

Number of `(cell, glacier)` rows.

*Not* the number of glaciers — use [`n_glaciers`](@ref) — and *not* the number of distinct
MODIS pixels either, since a pixel shared by two sub-cell glaciers appears twice. For the
download-relevant count use `length(rgi7_modis_unique_cells(cells))`.
"""
Base.length(t::RGI7ModisCells) = length(t.h)

"""
    n_glaciers(cells::RGI7ModisCells) -> Int

Number of RGI 7.0 glaciers represented.
"""
n_glaciers(t::RGI7ModisCells) = length(t.rgi_id)

function Base.show(io::IO, t::RGI7ModisCells)
    print(io, "RGI7ModisCells(", length(t), " rows over ", n_glaciers(t),
          " RGI 7.0 glaciers, ", length(unique(zip(t.h, t.v))), " tiles, ",
          count(t.forced), " forced)")
end

# Parsed once per session. 3.3M rows of text is a few seconds, and the albedo driver asks for
# the table repeatedly (once per year per hemisphere).
const _RGI7_CELLS_CACHE = Ref{Union{Nothing,RGI7ModisCells}}(nothing)

"""
    rgi7_modis_cells() -> RGI7ModisCells

Load (and cache) the vendored RGI 7.0 → MCD43A3 cell list.

Reads `data/rgi7_modis_cells.csv.gz` and `data/rgi7_glaciers.csv.gz`, both produced by
`data/make_rgi7_modis_cells.jl`. No network access.

# Example
```julia
cells = rgi7_modis_cells()
north, south = rgi7_hemisphere_split(cells)
lat, lon = rgi7_modis_cell_points(cells; index = north)
```
"""
function rgi7_modis_cells()
    cached = _RGI7_CELLS_CACHE[]
    cached === nothing || return cached

    for p in (_RGI7_CELLS_CSV, _RGI7_GLACIERS_CSV)
        isfile(p) || error("RGI7 cell table not found at $(p); generate it with " *
                           "data/make_rgi7_modis_cells.jl")
    end

    rgi_id = String[]
    area_km2 = Float32[]
    n_cells = Int32[]
    forced = Bool[]
    for vec in (rgi_id, area_km2, n_cells, forced)
        sizehint!(vec, _RGI7_N_GLACIERS)
    end
    open(GzipDecompressorStream, _RGI7_GLACIERS_CSV) do io
        for line in eachline(io)
            (isempty(line) || startswith(line, '#')) && continue
            startswith(line, "glacier_index,") && continue
            f = split(line, ',')
            length(f) == 5 || error("malformed row in $(_RGI7_GLACIERS_CSV): $line")
            # f[1] is the row number, carried in the file so the join key is explicit but
            # redundant on read: position in these vectors *is* the glacier index.
            push!(rgi_id, f[2])
            push!(area_km2, parse(Float32, f[3]))
            push!(n_cells, parse(Int32, f[4]))
            push!(forced, f[5] == "1")
        end
    end

    ng = length(rgi_id)
    # rgi7_glacier_cells binary-searches this, so an out-of-order regenerated file must fail
    # here rather than quietly return the wrong glacier.
    issorted(rgi_id) || error("$(_RGI7_GLACIERS_CSV) is not sorted by rgi_id; regenerate " *
                              "it with data/make_rgi7_modis_cells.jl")
    all(>=(Int32(1)), n_cells) || error("$(_RGI7_GLACIERS_CSV) has a glacier with no cell; " *
                                       "the generator's forced pass guarantees at least one, " *
                                       "so regenerate both tables")

    nc = Int(sum(n_cells))
    h = Vector{UInt8}(undef, nc)
    v = Vector{UInt8}(undef, nc)
    row = Vector{UInt16}(undef, nc)
    col = Vector{UInt16}(undef, nc)
    glacier = Vector{Int32}(undef, nc)
    i = 0
    open(GzipDecompressorStream, _RGI7_CELLS_CSV) do io
        for line in eachline(io)
            (isempty(line) || startswith(line, '#')) && continue
            startswith(line, "h,") && continue
            f = split(line, ',')
            length(f) == 5 || error("malformed row in $(_RGI7_CELLS_CSV): $line")
            i += 1
            i <= nc || error("$(_RGI7_CELLS_CSV) has more rows than the glacier table's " *
                             "n_cells sum ($(nc)); the two tables are out of sync, " *
                             "regenerate both with data/make_rgi7_modis_cells.jl")
            h[i] = parse(UInt8, f[1])
            v[i] = parse(UInt8, f[2])
            row[i] = parse(UInt16, f[3])
            col[i] = parse(UInt16, f[4])
            glacier[i] = parse(Int32, f[5])
        end
    end
    i == nc || error("$(_RGI7_CELLS_CSV) has $(i) rows but the glacier table's n_cells sum " *
                     "is $(nc); regenerate both with data/make_rgi7_modis_cells.jl")

    ng == _RGI7_N_GLACIERS ||
        @warn "unexpected RGI7 glacier count (expected $(_RGI7_N_GLACIERS) for v7.0)" ng

    tbl = RGI7ModisCells(h, v, row, col, glacier, rgi_id, area_km2, n_cells, forced)
    _RGI7_CELLS_CACHE[] = tbl
    return tbl
end

"""
    rgi7_modis_cell_tuples(cells; index = eachindex(cells)) -> Vector{NTuple{4,Int}}

The cell list as `(h, v, row, col)` tuples, ready for the cells-first method of
[`compute_glacier_ice_albedo_modis`](@ref).

Preferred over [`rgi7_modis_cell_points`](@ref) for that call: it skips the lon/lat round trip
and the 3.3M-entry deduplication `Dict` the point methods have to build.
"""
function rgi7_modis_cell_tuples(t::RGI7ModisCells; index = eachindex(t.h))
    return [(Int(t.h[i]), Int(t.v[i]), Int(t.row[i]), Int(t.col[i])) for i in index]
end

"""
    rgi7_modis_unique_cells(cells; index = eachindex(cells)) -> Vector{NTuple{4,Int}}

The **distinct** MODIS cells, in canonical order — the actual point list to hand to
[`compute_glacier_ice_albedo_modis`](@ref).

Use this rather than [`rgi7_modis_cell_tuples`](@ref) whenever the cells are being *sampled*:
a pixel shared by two sub-cell glaciers appears twice in the raw table, and sampling it twice
would double-count it in the accumulator and in the progress total. Cheap, because the table is
already sorted, so duplicates are adjacent.
"""
function rgi7_modis_unique_cells(t::RGI7ModisCells; index = eachindex(t.h))
    out = NTuple{4,Int}[]
    for i in index
        c = (Int(t.h[i]), Int(t.v[i]), Int(t.row[i]), Int(t.col[i]))
        (isempty(out) || out[end] != c) && push!(out, c)
    end
    return out
end

"""
    rgi7_modis_cell_points(cells; index = eachindex(cells)) -> (lat, lon)

Cell-**centre** latitude and longitude, for the lon/lat methods of
[`compute_glacier_ice_albedo_modis`](@ref) and for writing the output geometry.

Recomputed from `(h, v, row, col)` rather than stored, via [`_modis_cell_center`](@ref). Every
vendored cell is verified at generation time to satisfy
`_modis_cell(_modis_cell_center(cell...)...) === cell`, so feeding these points back in
reproduces exactly the same cells.
"""
function rgi7_modis_cell_points(t::RGI7ModisCells; index = eachindex(t.h))
    n = length(index)
    lat = Vector{Float64}(undef, n)
    lon = Vector{Float64}(undef, n)
    for (k, i) in enumerate(index)
        lat[k], lon[k] = _modis_cell_center(Int(t.h[i]), Int(t.v[i]),
                                            Int(t.row[i]), Int(t.col[i]))
    end
    return (lat, lon)
end

"""
    rgi7_glacier_cells(cells, rgi_id) -> Vector{Int}

Indices of the cells belonging to one glacier.

`rgi_id` is found by binary search, but the cells themselves are stored in canonical
`(h, v, row, col)` order rather than grouped by glacier, so a glacier's cells are **not
contiguous** and this scans the `glacier` column. That is a few milliseconds over 3.3M
`Int32`s — cheap enough not to justify carrying a second index in the vendored file.
"""
function rgi7_glacier_cells(t::RGI7ModisCells, rgi_id::AbstractString)
    gi = searchsortedfirst(t.rgi_id, rgi_id)
    (gi <= length(t.rgi_id) && t.rgi_id[gi] == rgi_id) ||
        throw(ArgumentError("glacier $(rgi_id) is not in the RGI 7.0 cell table"))
    return findall(==(Int32(gi)), t.glacier)
end

"""
    rgi7_hemisphere_split(cells) -> (north, south)

Partition the cell indices by hemisphere, as `v <= 8` and `v >= 9`.

An exact test, not an approximation: `_MODIS_UL_Y == 9 · _MODIS_TILE_SPAN_M` to the bit, so
`_modis_tile_edge_y(9)` is exactly `0.0` and no MODIS tile straddles the equator. That is why
this needs no coordinate computation at all.

The split is **required** before calling [`compute_glacier_ice_albedo_modis`](@ref) on a
global list. `_resolve_doy_range(:melt_season, lat)` returns `nothing` — the whole year — for
any list spanning both hemispheres, which more than doubles the download for no benefit. Call
once per half with an explicit `doy_range` instead.
"""
function rgi7_hemisphere_split(t::RGI7ModisCells)
    north = findall(<=(UInt8(8)), t.v)
    south = findall(>=(UInt8(9)), t.v)
    return (north, south)
end
