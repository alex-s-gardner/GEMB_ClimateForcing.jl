"""
Glacier bare-ice albedo sampled at a GeoInterface geometry.

`bare_ice_albedo(geometry, reducer)` reads the pooled tables written by
`data/run_rgi7_pooled_albedo.jl` — one row per MCD43A3 cell, holding the darkest-`percentile`
mean of every valid retrieval in the record — and returns the cells a point, line or polygon
selects.

The geometry is reprojected into the MODIS sinusoidal grid and burned there. Nothing is
resampled: the value reported for a cell is the value computed for that cell.
"""

using Rasters      # boolmask
import GeoInterface as GI

# Years the pooled record spans. Not a parameter of `bare_ice_albedo` — the product holds one
# value per cell over all of them — but recorded because it is what "the record" means, and a
# re-fold over a different span is a different statistic. MCD43A3 begins 2000-02-16, which is
# inside the southern melt-season window and before the northern one, so 2000 is a full season
# on both halves.
#
# Every year holds its full complement of dates, but yield falls off markedly: the fraction of
# cells resolving in a *single* year runs from 20.4 % (2002) to 9.1 % (2024) and 6.5 % (2025).
# That decline is retrieval quality, not missing input, and pooling across the record rather
# than reporting per-year values is the direct response to it.
const BARE_ICE_ALBEDO_YEARS = 2000:2025

# ------------------------------------------------------------------ geometry → sinusoidal

"""
    _bia_shape(geom) -> Symbol

`:point`, `:line` or `:polygon` — which burn `Rasters.boolmask` should apply to `geom`.

Derived from the GeoInterface trait rather than left to `boolmask`'s own inference so that a
geometry it cannot classify fails here, with a message naming the trait.
"""
function _bia_shape(geom)
    trait = GI.geomtrait(geom)
    trait isa GI.AbstractPointTrait && return :point
    trait isa GI.AbstractMultiPointTrait && return :point
    trait isa GI.AbstractCurveTrait && return :line
    trait isa GI.AbstractMultiCurveTrait && return :line
    trait isa GI.AbstractPolygonTrait && return :polygon
    trait isa GI.AbstractMultiPolygonTrait && return :polygon
    throw(ArgumentError(
        "cannot sample a $(trait) geometry; pass a point, multipoint, line, ring, polygon " *
        "or multipolygon, or set `shape` explicitly"))
end

"""
    _bia_project(geom) -> geom

`geom` with every vertex mapped from (lon, lat) degrees into MODIS sinusoidal metres by
[`_modis_lonlat_to_sinu`](@ref), rebuilt as the same GeoInterface geometry type.

The closed form is used rather than PROJ for the same reason `_rgi7_sinu_ring` does: it is the
*identical* projection [`_modis_cell`](@ref) and the vendored cell table were built with, so a
burned cell and a looked-up cell cannot disagree. Reprojecting the geometry — rather than
warping the albedo onto a lon/lat grid — is what keeps the reported value equal to the value
computed for that cell.

Edges are projected vertex-to-vertex, so a long edge chords across the curvature the
sinusoidal projection introduces. At glacier scale the error is far below the 463 m cell; for
edges spanning degrees, densify the geometry before calling.
"""
function _bia_project(geom)
    trait = GI.geomtrait(geom)
    if trait isa GI.AbstractPointTrait
        return GI.Point(_modis_lonlat_to_sinu(GI.x(geom), GI.y(geom)))
    elseif trait isa GI.AbstractMultiPointTrait
        return GI.MultiPoint([_modis_lonlat_to_sinu(GI.x(p), GI.y(p)) for p in GI.getpoint(geom)])
    elseif trait isa GI.LinearRingTrait
        return GI.LinearRing([_modis_lonlat_to_sinu(GI.x(p), GI.y(p)) for p in GI.getpoint(geom)])
    elseif trait isa GI.AbstractCurveTrait
        return GI.LineString([_modis_lonlat_to_sinu(GI.x(p), GI.y(p)) for p in GI.getpoint(geom)])
    elseif trait isa GI.AbstractPolygonTrait
        return GI.Polygon([_bia_project(r) for r in GI.getring(geom)])
    elseif trait isa GI.AbstractMultiCurveTrait
        return GI.MultiLineString([_bia_project(g) for g in GI.getgeom(geom)])
    elseif trait isa GI.AbstractMultiPolygonTrait
        return GI.MultiPolygon([_bia_project(g) for g in GI.getgeom(geom)])
    end
    throw(ArgumentError("cannot project a $(trait) geometry"))
end

"""
    _bia_points(geom)

Every vertex of `geom`, as an iterable of GeoInterface points.

`GI.getpoint` has no method for `PointTrait` — a point is not a collection of points — so a
bare point is wrapped in a 1-tuple. Without this, sampling at a single point throws a
`MethodError` from inside GeoInterface's fallbacks.
"""
_bia_points(geom) = GI.geomtrait(geom) isa GI.AbstractPointTrait ? (geom,) : GI.getpoint(geom)

"""
    _bia_check_longitude(geom)

Reject a geometry spanning more than 180° of longitude.

The sinusoidal projection maps each parallel to a line segment centred on `lon_0 = 0`, so an
antimeridian-crossing ring projects to a shape sweeping the whole globe rather than the narrow
one intended, and would silently select cells across every intervening tile. Split such a
geometry at the antimeridian and call once per part. Same guard, same reason, as
`_rgi7_sinu_ring`.
"""
function _bia_check_longitude(geom)
    lo, hi = Inf, -Inf
    for p in _bia_points(geom)
        x = GI.x(p)
        lo = min(lo, x)
        hi = max(hi, x)
    end
    isfinite(lo) || throw(ArgumentError("geometry has no vertices"))
    hi - lo > 180 && throw(ArgumentError(
        "geometry spans $(round(hi - lo; digits=1))° of longitude " *
        "($(round(lo; digits=3))…$(round(hi; digits=3))), which the sinusoidal projection " *
        "cannot represent as one shape. Split it at the antimeridian."))
    return nothing
end

# ----------------------------------------------------------------------- the burn

"""
    _bia_tile_windows(sinu_geom; pad=1) -> Vector{Tuple{Int,Int,UnitRange{Int},UnitRange{Int}}}

`(h, v, rows, cols)` for every MODIS tile the projected geometry's bounding box touches.

The window is the geometry's bounding box, not the whole tile, so the burn target stays small
for a glacier-scale geometry (a 2400×2400 tile mask is 5.8 M cells). `pad` widens it by whole
cells so that `boundary = :touches` — which selects cells the geometry merely clips — cannot
be cut off at the window edge.
"""
function _bia_tile_windows(sinu_geom; pad::Integer = 1)
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    for p in _bia_points(sinu_geom)
        x, y = GI.x(p), GI.y(p)
        xmin, xmax = min(xmin, x), max(xmax, x)
        ymin, ymax = min(ymin, y), max(ymax, y)
    end

    hlo = clamp(floor(Int, (xmin - _MODIS_UL_X) / _MODIS_TILE_SPAN_M), 0, _MODIS_H_MAX)
    hhi = clamp(floor(Int, (xmax - _MODIS_UL_X) / _MODIS_TILE_SPAN_M), 0, _MODIS_H_MAX)
    vlo = clamp(floor(Int, (_MODIS_UL_Y - ymax) / _MODIS_TILE_SPAN_M), 0, _MODIS_V_MAX)
    vhi = clamp(floor(Int, (_MODIS_UL_Y - ymin) / _MODIS_TILE_SPAN_M), 0, _MODIS_V_MAX)

    out = Tuple{Int,Int,UnitRange{Int},UnitRange{Int}}[]
    for v in vlo:vhi, h in hlo:hhi
        ulx, uly = _modis_tile_origin(h, v)
        clo = floor(Int, (xmin - ulx) / _MODIS_PIXEL_M) + 1 - pad
        chi = floor(Int, (xmax - ulx) / _MODIS_PIXEL_M) + 1 + pad
        rlo = floor(Int, (uly - ymax) / _MODIS_PIXEL_M) + 1 - pad
        rhi = floor(Int, (uly - ymin) / _MODIS_PIXEL_M) + 1 + pad
        cols = clamp(clo, 1, _MODIS_TILE_PIXELS):clamp(chi, 1, _MODIS_TILE_PIXELS)
        rows = clamp(rlo, 1, _MODIS_TILE_PIXELS):clamp(rhi, 1, _MODIS_TILE_PIXELS)
        (isempty(cols) || isempty(rows)) && continue
        push!(out, (h, v, rows, cols))
    end
    return out
end

"""
    _bia_window_dims(h, v, rows, cols) -> (X, Y)

Burn-target dimensions for one tile window, sliced out of [`_rgi7_tile_dims`](@ref).

Taken from that function rather than recomputed so the coordinates are **bit-identical** to
the cell centres the vendored table was burned against, and so `Rasters`' inside test — which
it evaluates at the lookup values — lands on the same cell centres. Building them from
`_MODIS_UL_X + (h * 2400 + col - 0.5) * pixel` instead would reintroduce exactly the `ulx + T`
rounding that [`_modis_tile_edge_x`](@ref) exists to avoid.
"""
function _bia_window_dims(h::Integer, v::Integer, rows::AbstractUnitRange,
                          cols::AbstractUnitRange)
    dx, dy = _rgi7_tile_dims(h, v)
    xs = parent(parent(dx))   # Dimension -> Sampled lookup -> the centre Vector itself
    ys = parent(parent(dy))
    return (X(Sampled(xs[cols]; order = ForwardOrdered(),
                      span = Regular(_MODIS_PIXEL_M), sampling = Intervals(Center()))),
            Y(Sampled(ys[rows]; order = ReverseOrdered(),
                      span = Regular(-_MODIS_PIXEL_M), sampling = Intervals(Center()))))
end

"""
    _bia_burn(sinu_geom, shape, boundary) -> Vector{NTuple{4,Int}}

The `(h, v, row, col)` cells a projected geometry selects, in canonical sorted order.

`Rasters.boolmask` does the burn, per tile window. It — not `Rasters.extract` — is the
primitive here for two reasons. `extract` **throws on a multipoint geometry** in Rasters
0.15.0 (`TypeError: non-boolean (Int64) used in boolean context`, from the `::Bool` assertion
in its `AbstractMultiPointTrait` method), whereas `boolmask` handles every trait; and
`extract` cannot index a 3-D `(X, Y, Ti)` raster at all, failing in the `Extractor`
constructor. Where both work they agree cell-for-cell, which a test pins.

Cells are deduplicated by construction — a mask holds each cell once — so two points landing
in one 463 m cell select it once, matching `_modis_dedup_points`.
"""
function _bia_burn(sinu_geom, shape::Symbol, boundary::Symbol)
    cells = NTuple{4,Int}[]
    for (h, v, rows, cols) in _bia_tile_windows(sinu_geom)
        dims = _bia_window_dims(h, v, rows, cols)
        mask = shape === :polygon ? boolmask(sinu_geom; to = dims, shape, boundary) :
                                    boolmask(sinu_geom; to = dims, shape)
        for I in findall(parent(mask))
            push!(cells, (h, v, rows[I[2]], cols[I[1]]))
        end
    end
    sort!(cells)
    return cells
end

# ------------------------------------------------------------------- reading the values

"""
    _bia_next_uint(line, i) -> (value, next_i)

Parse one unsigned decimal field of `line` starting at byte `i`, returning it and the byte
after the delimiter.

Exists so the per-year files can be filtered **without `split`**. Every row of a hemisphere
file is a candidate but only a handful match a glacier-scale geometry, and `split(line, ',')`
allocates a `Vector` of eight `SubString`s for each — 65 M rows for a 25-year query, which
measured at 38.5 GiB allocated and 59 % GC time. Reading `h` and `v` this way rejects a
foreign tile after ~4 bytes and allocates nothing.

Valid only for the non-negative integer key columns (`h`, `v`, `row`, `col`); the albedo
columns are floats and are parsed by `split` on the rows that survive, where the cost is
irrelevant.
"""
@inline function _bia_next_uint(buf::AbstractVector{UInt8}, i::Int, stop::Int)
    v = 0
    while i <= stop
        c = buf[i]
        c == UInt8(',') && break
        v = 10 * v + Int(c - 0x30)
        i += 1
    end
    return v, i + 1
end

"""
    _bia_hemisphere(v) -> Symbol

Which per-year file holds tile row `v`.

`_MODIS_UL_Y == 9 * _MODIS_TILE_SPAN_M` to the bit, so `_modis_tile_edge_y(9) === 0.0` and no
tile straddles the equator: `v <= 8` *is* the northern hemisphere.
"""
_bia_hemisphere(v::Integer) = v <= 8 ? :north : :south


"""
    _bia_read(cells, dir) -> NamedTuple

Pooled values for `cells`, as one vector per column of the pooled tables.

Only the hemisphere files the cells actually occupy are opened, so a northern geometry never
touches a southern file. Rows are located in `cells` by `searchsortedfirst`, not by position:
the files and `cells` share a canonical sort, but a file holds every glacierized cell on its
half of the planet while `cells` holds a handful, so nothing about the order can be assumed.
"""
function _bia_read(cells::Vector{NTuple{4,Int}}, dir::AbstractString)
    ncell = length(cells)
    hemis = unique(_bia_hemisphere(c[2]) for c in cells)
    for hemi in hemis
        path = bare_ice_albedo_path(dir, hemi)
        isfile(path) || throw(ArgumentError(
            "no pooled bare-ice albedo file at $(path). Generate it with " *
            "`julia --project=. data/run_rgi7_pooled_albedo.jl $(hemi)` — that re-folds the " *
            "cached samples and downloads nothing."))
    end

    # The tiles the cells occupy, for the cheap prefilter in the read loop below.
    tiles = Set((c[1], c[2]) for c in cells)
    bsa, wsa = fill(NaN32, ncell), fill(NaN32, ncell)
    nb, nw = zeros(Int32, ncell), zeros(Int32, ncell)
    kb, kw = zeros(Int32, ncell), zeros(Int32, ncell)
    # Whether a cell appears in the product at all, which is *not* the same as having a
    # resolved albedo: a glacierized cell too cloudy to resolve is present with a NaN, while a
    # cell off the RGI 7.0 outlines is absent entirely. Separating the two is what lets a
    # caller tell "not glacier" from "glacier, no data".
    found = falses(ncell)

    for hemi in hemis
        # Decompressed whole rather than iterated with `eachline`, which allocates a `String`
        # per row: a hemisphere file is 2.58 M rows, so line-based reading built millions of
        # throwaway strings and spent most of its time in GC.
        buf = open(io -> read(GzipDecompressorStream(io)), bare_ice_albedo_path(dir, hemi))
        n = length(buf)
        i = 1
        while i <= n
            e = something(findnext(isequal(UInt8('\n')), buf, i), n + 1)
            stop = e - 1
            c1 = buf[i]
            if c1 != UInt8('#') && c1 != UInt8('h')
                h, j = _bia_next_uint(buf, i, stop)
                v, j = _bia_next_uint(buf, j, stop)
                if (h, v) in tiles
                    row, j = _bia_next_uint(buf, j, stop)
                    col, j = _bia_next_uint(buf, j, stop)
                    cell = (h, v, row, col)
                    k = searchsortedfirst(cells, cell)
                    if k <= ncell && cells[k] == cell
                        # Only the surviving rows pay for string handling, which for a
                        # glacier-scale geometry is a handful out of millions.
                        f = split(String(@view buf[j:stop]), ',')
                        length(f) == 6 || error(
                            "$(bare_ice_albedo_path(dir, hemi)) row for cell $(cell) has " *
                            "$(length(f) + 4) columns, expected 10. The file was written by a " *
                            "different version of data/run_rgi7_pooled_albedo.jl.")
                        found[k] = true
                        bsa[k] = parse(Float32, f[1])
                        wsa[k] = parse(Float32, f[2])
                        nb[k] = parse(Int32, f[3])
                        nw[k] = parse(Int32, f[4])
                        kb[k] = parse(Int32, f[5])
                        kw[k] = parse(Int32, f[6])
                    end
                end
            end
            i = e + 1
        end
    end
    return (; bsa, wsa, nb, nw, kb, kw, found)
end

# --------------------------------------------------------------------------- public API

"""
    bare_ice_albedo_path(dir, hemisphere) -> String

Path of one pooled hemisphere table, matching what `data/run_rgi7_pooled_albedo.jl` writes.
"""
bare_ice_albedo_path(dir::AbstractString, hemisphere::Symbol) =
    joinpath(dir, "rgi7_ice_albedo_pooled_$(hemisphere).csv.gz")

"""
    bare_ice_albedo(geometry, reducer = nothing; kwargs...) -> DimStack or NamedTuple

Glacier bare-ice albedo at the MCD43A3 cells a geometry selects.

Each cell carries **one** value: the mean of the darkest `percentile` (5 %) of *every* valid
retrieval over the whole 2000–2025 record, with the number of retrievals behind it. Pooled over
the record, not per year — year-to-year coverage is far too variable for an annual value to be
comparable, and `:n_valid_*` is what qualifies the albedo.

`geometry` is any GeoInterface geometry in **(lon, lat) degrees** — point, multipoint, line,
ring, polygon or multipolygon. It is reprojected into the MODIS sinusoidal grid and burned
there by `Rasters.boolmask`, so the value reported for a cell is the value computed for that
cell: nothing is interpolated, warped or resampled.

Cell selection is nearest-cell throughout. A point (or each point of a multipoint) selects the
cell containing it; a line selects every cell it passes through; a polygon selects every cell
whose *centre* it covers, or a different rule via `boundary`.

`reducer` controls the shape of the result:

- `nothing` (default) — a `DimStack` over `Dim{:cell}`, one entry per selected cell.
- a function (`median`, `mean`, `minimum`, `x -> quantile(x, 0.1)`, …) — a `NamedTuple` of
  scalars, the reducer applied across the selected cells. It never sees a `NaN`: unresolved
  cells are dropped first, which is what makes `n_cells_bsa` meaningful.

# Keywords
- `dir = joinpath(@__DIR__, "..", "data")`: where the pooled tables live.
- `boundary = :center`: for polygons, `:center` (cell centre inside), `:touches` (any overlap)
  or `:inside` (wholly within). `:center` matches how the vendored cell table itself was
  burned; `:touches` adds perimeter cells that are majority off-glacier and therefore biased
  dark, which the darkest-percentile statistic amplifies.
- `shape`: override the `:point` / `:line` / `:polygon` burn inferred from the geometry type.
- `min_cells = 1`: with a `reducer`, the fewest resolved cells needed before reducing at all.

No retrieval-count threshold is applied. `:n_valid_bsa` / `:n_valid_wsa` report how many
retrievals passed QA and the albedo range, and `:k_used_bsa` / `:k_used_wsa` how many of the
darkest were averaged, so filtering on sample count is the caller's to do.

# Returns
With `reducer === nothing`, a `DimStack` over `Dim{:cell}` carrying `:albedo_bsa` (black-sky)
and `:albedo_wsa` (white-sky), their `:n_valid_*` and `:k_used_*` counts, `:latitude`,
`:longitude`, `:cell_id`, and `:in_product` (false means the cell is off-glacier, not merely
cloudy). Albedo is `NaN32` where no retrieval survived QC.

[`bare_ice_albedo_points`](@ref) turns the result into the cell centres as geometries, for
sampling another dataset at these same cells.

# Example
```julia
using GEMB_ClimateForcing, Statistics
import GeoInterface as GI

a = bare_ice_albedo(GI.Point(-49.5, 69.2))
a[:albedo_bsa][1]                    # the cell's bare-ice albedo
a[:n_valid_bsa][1]                   # retrievals behind it

poly = GI.Polygon([GI.LinearRing([(-50.0, 69.0), (-49.0, 69.0),
                                  (-49.0, 69.5), (-50.0, 69.5), (-50.0, 69.0)])])
bare_ice_albedo(poly, median)        # (albedo_bsa = 0.44f0, n_cells_bsa = 812, …)
```

!!! warning "Albedo is not clamped to 1"
    MCD43A3's valid range reaches 32766, so a legitimate retrieval can exceed 1.0. The QC range
    (`POOLED_ICE_ALBEDO_RANGE`, `(0.25, 1.0)`) is a *glaciological* filter, not a physical
    clamp. Values are unsuitable as-is for GEMB's `albedo_ice`, which asserts
    `0.2 ≤ albedo_ice ≤ 0.6`.

See also [`pool_ice_albedo_from_cache`](@ref), which computes this product, and
[`compute_glacier_ice_albedo_modis`](@ref) for the per-year statistic at an arbitrary point list.
"""
function bare_ice_albedo(geometry, reducer = nothing;
                         dir::AbstractString = joinpath(@__DIR__, "..", "data"),
                         boundary::Symbol = :center,
                         shape::Union{Symbol,Nothing} = nothing,
                         min_cells::Integer = 1)
    GI.isgeometry(geometry) || throw(ArgumentError(
        "expected a GeoInterface geometry, got a $(typeof(geometry)). Wrap coordinates with " *
        "e.g. `GeoInterface.Point(lon, lat)` or `GeoInterface.Polygon([GeoInterface.LinearRing(pts)])`."))
    boundary in (:center, :touches, :inside) || throw(ArgumentError(
        "boundary must be :center, :touches or :inside, got :$(boundary)"))
    min_cells >= 1 || throw(ArgumentError("min_cells must be >= 1, got $(min_cells)"))

    burn = something(shape, _bia_shape(geometry))
    _bia_check_longitude(geometry)
    cells = _bia_burn(_bia_project(geometry), burn, boundary)
    isempty(cells) && throw(ArgumentError(
        "the geometry selects no MCD43A3 cell. A geometry smaller than one 463 m cell " *
        "selects nothing unless it covers a cell centre; use `boundary = :touches`, or pass " *
        "the cell centre as a point."))

    v = _bia_read(cells, dir)
    if !any(v.found)
        @warn("none of the $(length(cells)) selected cells is in the RGI 7.0 bare-ice albedo " *
              "product, so every value is NaN — the geometry is probably off-glacier",
              cells = length(cells))
    end

    if reducer !== nothing
        out = Pair{Symbol,Any}[]
        for (alb, nval, kused, label) in ((v.bsa, v.nb, v.kb, "bsa"), (v.wsa, v.nw, v.kw, "wsa"))
            keep = findall(!isnan, alb)
            red = length(keep) >= min_cells ? Float32(reducer(alb[keep])) : NaN32
            push!(out, Symbol("albedo_", label) => red)
            push!(out, Symbol("n_cells_", label) => length(keep))
            push!(out, Symbol("n_valid_", label) => sum(Int, nval[keep]; init = 0))
            push!(out, Symbol("k_used_", label) => sum(Int, kused[keep]; init = 0))
        end
        push!(out, :n_cells => length(cells))
        push!(out, :n_cells_in_product => count(v.found))
        push!(out, :reduction => _rgi7_reduction_name(reducer))
        return NamedTuple(out)
    end

    cell_dim = Dim{:cell}(1:length(cells))
    alb_meta(label) = Dict{String,Any}(
        "units" => "1",
        "long_name" => "bare-ice albedo, darkest-percentile mean over the whole record ($(label))",
        "source" => "MODIS MCD43A3 v061, pooled 2000-2025",
        "note" => "NOT clamped to 1; pooled over all dates, not a per-year value")
    out = Pair{Symbol,Any}[]
    for (alb, nval, kused, label) in ((v.bsa, v.nb, v.kb, "bsa"), (v.wsa, v.nw, v.kw, "wsa"))
        push!(out, Symbol("albedo_", label) => DimArray(alb, (cell_dim,);
            metadata = alb_meta(label)))
        push!(out, Symbol("n_valid_", label) => DimArray(nval, (cell_dim,);
            metadata = Dict("units" => "count",
                            "long_name" => "retrievals passing QA and the albedo range")))
        push!(out, Symbol("k_used_", label) => DimArray(kused, (cell_dim,);
            metadata = Dict("units" => "count",
                            "long_name" => "darkest retrievals averaged into the albedo")))
    end
    centres = [_modis_cell_center(c...) for c in cells]
    push!(out, :latitude => DimArray([c[1] for c in centres], (cell_dim,);
        metadata = Dict("units" => "degrees_north",
                        "long_name" => "latitude of the MODIS cell centre")))
    push!(out, :longitude => DimArray([c[2] for c in centres], (cell_dim,);
        metadata = Dict("units" => "degrees_east",
                        "long_name" => "longitude of the MODIS cell centre")))
    push!(out, :cell_id => DimArray(_modis_cell_id.(cells), (cell_dim,);
        metadata = Dict("long_name" => "MCD43A3 sinusoidal grid cell")))
    push!(out, :in_product => DimArray(collect(v.found), (cell_dim,);
        metadata = Dict("long_name" => "cell is glacierized in RGI 7.0 and present in the " *
                                       "product (false means off-glacier, not cloudy)")))

    return DimStack(NamedTuple(out); metadata = Dict{String,Any}(
        "source" => "MODIS MCD43A3 v061 bare-ice albedo pooled over 2000-2025, at RGI 7.0 cells",
        "n_cells" => length(cells),
        "n_cells_in_product" => count(v.found),
        "shape" => String(burn),
        "boundary" => String(boundary),
        "grid" => _MODIS_SINU_PROJ,
        "note" => "Cell values are native-grid; the geometry was reprojected, not the albedo.",
    ))
end

"""
    bare_ice_albedo_points(a) -> Vector{GeoInterface.Point}

The MCD43A3 cell centres of a [`bare_ice_albedo`](@ref) result, as `(lon, lat)` points in
EPSG:4326 and in the same order as its `Dim{:cell}` axis.

This is what samples another dataset at the same cells — `Rasters.extract` and everything else
taking GeoInterface input dispatch on geometries, so no coordinate handling is needed here:

```julia
a = bare_ice_albedo(poly)
# Crop before extracting. Sampling a lazy global mosaic costs one HTTP read per point; cropping
# to the geometry's window first turns 258 points from ~90 s into milliseconds.
window = read(view(dem, X = (-141.06 .. -140.94), Y = (60.44 .. 60.56)))
elevation = extract(window, bare_ice_albedo_points(a); geometry = false)
```

A plain `Vector` of `Point`s, deliberately, rather than a `MultiPoint` or a `DimArray`. Rasters
0.15.0 reads a `DimArray` as a *table* and rejects it for having no `geometry` column, and
`extract` throws outright on a `MultiPoint` — the reason [`_bia_burn`](@ref) uses `boolmask`
instead.
"""
bare_ice_albedo_points(a) =
    [GI.Point(lon, lat) for (lon, lat) in zip(parent(a[:longitude]), parent(a[:latitude]))]
