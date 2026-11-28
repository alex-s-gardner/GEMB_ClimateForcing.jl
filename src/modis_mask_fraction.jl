"""
    modis_mask_fraction(cells, mask, classes; factor = 10, buffer = 0, crs = Rasters.crs(mask)) -> Vector{Float32}
    modis_mask_fraction(lat, lon, mask, classes; kwargs...) -> Vector{Float32}

Share (0–1) of each MCD43A3 500 m cell covered by the categories `classes` of a categorical
`mask` raster on a regular grid, for example BedMachine's ice mask (`classes = (2, 3)`, grounded
and floating ice).

The cell is split into `factor × factor` sub-cells in the MODIS sinusoidal grid, the sub-cell
geometry of the RGI 7.0 fraction behind [`modis_glacier_fraction`](@ref), and each sub-cell centre
takes the value of the mask pixel containing it. A sub-cell centre outside the mask, or on its
missing value, counts as not in `classes`. The sinusoidal grid's polar cells reach past ±180° of
longitude, off the Earth's surface; the share is taken over the sub-cells on the surface. With
`buffer` [m] the sampled square is the cell widened by `buffer` on every side at the same sub-cell
spacing, so a fraction of 1 means the mask holds `classes` everywhere within `buffer` of the cell:
use it where the mask pixel is as coarse as the cell, which makes a within-cell fraction
meaningless.

`cells` are `(h, v, row, col)` tuples, or the cells containing points `lat`, `lon` (degrees).
`crs` is the mask's coordinate reference system; it must be given when the raster carries none.
"""
function modis_mask_fraction(cells::AbstractVector{<:NTuple{4,Integer}}, mask::Rasters.AbstractRaster, classes;
                             factor::Integer = 10, buffer::Real = 0, crs = Rasters.crs(mask))
    factor >= 1 || throw(ArgumentError("factor must be >= 1, got $(factor)"))
    buffer >= 0 || throw(ArgumentError("buffer must be >= 0 m, got $(buffer)"))
    crs === nothing && throw(ArgumentError(
        "the mask raster carries no CRS; pass it, e.g. `crs = Rasters.EPSG(3031)`"))
    n = ceil(Int, factor * (1 + 2 * buffer / _MODIS_PIXEL_M))
    d = (_MODIS_PIXEL_M + 2 * buffer) / n
    # Index of the mask pixel holding a coordinate, from the regular grid's first pixel centre and step.
    xs, ys = lookup(mask, X), lookup(mask, Y)
    c1x, c1y = _first_centre(xs), _first_centre(ys)
    index(l, c1, v) = (i = round(Int, (v - c1) / step(l)) + 1; firstindex(l) <= i <= lastindex(l) ? i : nothing)
    pts = Vector{Tuple{Float64,Float64}}(undef, n * n)
    ArchGDAL.importEPSG(4326; order = :trad) do src
        ArchGDAL.importCRS(crs; order = :trad) do dst
            ArchGDAL.createcoordtrans(src, dst) do tr
                map(cells) do (h, v, row, col)
                    ulx, uly = _modis_tile_origin(h, v)
                    x0 = ulx + (col - 1) * _MODIS_PIXEL_M - buffer
                    y0 = uly - (row - 1) * _MODIS_PIXEL_M + buffer
                    for j in 1:n, i in 1:n
                        lon, lat = _modis_sinu_to_lonlat(x0 + (i - 0.5) * d, y0 - (j - 0.5) * d)
                        # Near the poles a sinusoidal cell reaches past ±180° of longitude, or to ±90° of
                        # latitude at the grid's edge, off the Earth's surface; such a sub-cell is left out.
                        if abs(lon) > 180 || abs(lat) >= 90
                            pts[(j - 1) * n + i] = (NaN, NaN)
                            continue
                        end
                        p = ArchGDAL.transform!(ArchGDAL.createpoint(lon, lat), tr)
                        pts[(j - 1) * n + i] = (ArchGDAL.getx(p, 0), ArchGDAL.gety(p, 0))
                    end
                    ij = [isnan(x) ? (nothing, nothing) : (index(xs, c1x, x), index(ys, c1y, y)) for (x, y) in pts]
                    ok = [i !== nothing && j !== nothing for (i, j) in ij]
                    any(ok) || return 0.0f0
                    ir = extrema(first(ij[k]) for k in eachindex(ij) if ok[k])
                    jr = extrema(last(ij[k]) for k in eachindex(ij) if ok[k])
                    w = mask[X(ir[1]:ir[2]), Y(jr[1]:jr[2])]       # one read of the cell's window
                    mv = missingval(mask)
                    hits = count(eachindex(ij)) do k
                        ok[k] || return false
                        val = w[ij[k][1] - ir[1] + 1, ij[k][2] - jr[1] + 1]
                        !ismissing(val) && val !== mv && val in classes
                    end
                    Float32(hits / count(p -> !isnan(first(p)), pts))
                end
            end
        end
    end
end
function modis_mask_fraction(lat::AbstractVector{<:Real}, lon::AbstractVector{<:Real}, mask::Rasters.AbstractRaster,
                             classes; kwargs...)
    length(lat) == length(lon) || throw(DimensionMismatch(
        "lat and lon must have the same length, got $(length(lat)) and $(length(lon))"))
    return modis_mask_fraction([_modis_cell(la, lo) for (la, lo) in zip(lat, lon)], mask, classes; kwargs...)
end

"""
    _first_centre(l) -> Real

Centre of the first pixel of a regular lookup: its first value is the pixel centre for points or
centred intervals, and an edge for intervals with a start or end locus.
"""
function _first_centre(l)
    Rasters.DimensionalData.Lookups.isregular(l) || throw(ArgumentError("the mask must be on a regular grid"))
    s = Rasters.DimensionalData.Lookups.sampling(l)
    s isa Rasters.DimensionalData.Lookups.Points && return first(l)
    loc = Rasters.DimensionalData.Lookups.locus(s)
    return loc isa Rasters.DimensionalData.Lookups.Start ? first(l) + step(l) / 2 :
           loc isa Rasters.DimensionalData.Lookups.End ? first(l) - step(l) / 2 : first(l)
end
