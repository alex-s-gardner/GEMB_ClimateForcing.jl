# Glacier fraction and glacier elevation of every MCD43A3 500 m pixel in the RGI 7.0 sample cache.
#
# Each MODIS sinusoidal tile that holds a cache pixel is rasterized at FACTOR × the MODIS resolution
# (default 10, ~46 m sub-cells) from the RGI 7.0 outlines reprojected to MODIS sinusoidal metres, and the
# sub-cells are averaged to each 500 m pixel: glacier_fraction is the share of a pixel's sub-cells whose
# centres lie inside an outline. The glacier elevation is the mean Copernicus GLO-30 DEM height over
# those glacier sub-cells, sampled at sub-cell centres.
#
# Writes, per hemisphere and in canonical (h, v, row, col) order, one row per cache pixel:
#   <out>/glacier_fraction_<hemisphere>.csv
#       h, v, row, col, lat, lon, glacier_fraction, elevation_glacier
#   <out>/rgi7_modis_glacier_fraction_<hemisphere>.csv.gz   (the table vendored in data/)
#       h, v, row, col, glacier_percent   (glacier_fraction × 100, rounded; exact at FACTOR = 10)
#
#   julia -t 32 --project=<GEMB_ClimateForcing> data/make_rgi7_modis_glacier_fraction.jl <rgi7 gpkg> <dem root> <out dir>
using GEMB_ClimateForcing, Rasters, ArchGDAL, CodecZlib, Dates, Statistics, Printf
using DimensionalData.Lookups
const G = GEMB_ClimateForcing

rgi_path, dem_root, out = ARGS
mkpath(out)
const FACTOR = parse(Int, get(ENV, "FACTOR", "10"))
say(x...) = (println(x...); flush(stdout))

# --- RGI 7.0 outlines in MODIS sinusoidal metres, grouped by the MODIS tile they reach -------------
t0 = time()
by_tile = Dict{Tuple{Int,Int},Vector{ArchGDAL.IGeometry}}()
ArchGDAL.read(rgi_path) do ds
    layer = ArchGDAL.getlayer(ds, 0)
    for feat in layer
        g = ArchGDAL.getgeom(feat)
        id = string(ArchGDAL.getfield(feat, "rgi_id"))
        sg = G._rgi7_sinu_geometry(g, id)
        for ht in G._rgi7_candidate_tiles(sg)
            push!(get!(by_tile, ht, ArchGDAL.IGeometry[]), sg)
        end
    end
end
say("outlines grouped into ", length(by_tile), " MODIS tiles in ", round(time() - t0), " s")

index_tiles = Set{String}(strip(l) for l in eachline(joinpath(dem_root, "tileList.txt")) if !isempty(strip(l)))

# Sub-cell grid of one MODIS tile: FACTOR×FACTOR sub-cells per 500 m pixel.
function subgrid_dims(h, v)
    ulx, uly = G._modis_tile_origin(h, v)
    n = G._MODIS_TILE_PIXELS * FACTOR
    d = G._MODIS_PIXEL_M / FACTOR
    xs = [ulx + (c - 0.5) * d for c in 1:n]
    ys = [uly - (r - 0.5) * d for r in 1:n]
    return (X(Sampled(xs; order = ForwardOrdered(), span = Regular(d), sampling = Intervals(Center()))),
            Y(Sampled(ys; order = ReverseOrdered(), span = Regular(-d), sampling = Intervals(Center()))))
end

t = G.rgi7_modis_cells()
north, south = G.rgi7_hemisphere_split(t)
for (hemi, idx) in ((:north, north), (:south, south))
    cells = G.rgi7_modis_unique_cells(t; index = idx)
    frac = fill(NaN32, length(cells)); zgl = fill(NaN32, length(cells))
    lat = fill(NaN, length(cells)); lon = fill(NaN, length(cells))
    tiles_here = sort(unique((c[1], c[2]) for c in cells))
    only_tiles = get(ENV, "ONLY_TILES", "")
    if !isempty(only_tiles)
        want = Set(Tuple(parse.(Int, split(x, ':'))) for x in split(only_tiles, ','))
        tiles_here = filter(in(want), tiles_here)
        isempty(tiles_here) && continue
    end
    pos_of_tile = Dict{Tuple{Int,Int},Vector{Int}}()
    for (i, c) in enumerate(cells); push!(get!(pos_of_tile, (c[1], c[2]), Int[]), i); end
    say(hemi, ": ", length(cells), " pixels in ", length(tiles_here), " MODIS tiles")
    lk = ReentrantLock()
    Threads.@threads for ht in tiles_here
        h, v = ht
        geoms = get(by_tile, ht, ArchGDAL.IGeometry[])
        dims = subgrid_dims(h, v)
        grid = Raster(falses(length(dims[1]), length(dims[2])), dims; missingval = false)
        isempty(geoms) || Rasters.rasterize!(last, grid, geoms; fill = true, boundary = :center,
                                             verbose = false, progress = false, threaded = false)
        idxs = pos_of_tile[ht]
        # DEM for this tile's pixels: the lon/lat envelope of their pixel centres, padded.
        lonlat = [G._modis_cell_center(cells[i]...) for i in idxs]       # (lat, lon)
        la = first.(lonlat); lo = last.(lonlat)
        ids = G._copernicus_dem_tiles_for_extent(minimum(lo) - 0.02, maximum(lo) + 0.02,
                                                 minimum(la) - 0.02, maximum(la) + 0.02, index_tiles)
        dem = nothing; dlon = 0.0; dlat = 0.0; lon0 = 0.0; lat0 = 0.0
        if !isempty(ids)
            vrt = G._copernicus_dem_write_vrt(ids, joinpath(out, "dem_$(h)_$(v).vrt"),
                                              [joinpath(dem_root, id, id * ".tif") for id in ids])
            r = read(view(Raster(vrt; lazy = true), X = (minimum(lo) - 0.02) .. (maximum(lo) + 0.02),
                          Y = (minimum(la) - 0.02) .. (maximum(la) + 0.02)))
            dem = parent(r); lx = lookup(r, X); ly = lookup(r, Y)
            lon0, lat0 = first(lx), first(ly); dlon = step(lx); dlat = step(ly)
        end
        demz(lon_s, lat_s) = begin
            i = round(Int, (lon_s - lon0) / dlon) + 1; j = round(Int, (lat_s - lat0) / dlat) + 1
            (1 <= i <= size(dem, 1) && 1 <= j <= size(dem, 2)) ? dem[i, j] : missing
        end
        xs = lookup(dims[1]); ys = lookup(dims[2])
        for (q, i) in enumerate(idxs)
            _, _, r_, c_ = cells[i]
            ic = ((c_ - 1) * FACTOR + 1):(c_ * FACTOR); ir = ((r_ - 1) * FACTOR + 1):(r_ * FACTOR)
            sub = view(grid, ic, ir)
            n_gl = count(sub)
            zsum = 0.0; zn = 0
            if dem !== nothing && n_gl > 0
                for (jj, j) in enumerate(ir), (kk, k) in enumerate(ic)
                    sub[kk, jj] || continue
                    lon_s, lat_s = G._modis_sinu_to_lonlat(xs[k], ys[j])
                    z = demz(lon_s, lat_s)
                    (ismissing(z) || !isfinite(z) || z < -500) && continue
                    zsum += z; zn += 1
                end
            end
            frac[i] = n_gl / FACTOR^2; zgl[i] = zn > 0 ? zsum / zn : NaN32
            lat[i], lon[i] = lonlat[q]
        end
        lock(lk) do; say(@sprintf("  tile h%02dv%02d: %d pixels, %d outlines", h, v, length(idxs), length(geoms))); end
    end
    open(joinpath(out, "glacier_fraction_$(hemi).csv"), "w") do io
        println(io, "h,v,row,col,lat,lon,glacier_fraction,elevation_glacier")
        for (i, c) in enumerate(cells)
            isempty(only_tiles) || (c[1], c[2]) in tiles_here || continue
            @printf(io, "%d,%d,%d,%d,%.5f,%.5f,%.4f,%.1f\n", c..., lat[i], lon[i], frac[i], zgl[i])
        end
    end
    isempty(only_tiles) && open(GzipCompressorStream,
                                joinpath(out, "rgi7_modis_glacier_fraction_$(hemi).csv.gz"), "w") do io
        println(io, "h,v,row,col,glacier_percent")
        for (i, c) in enumerate(cells)
            @printf(io, "%d,%d,%d,%d,%d\n", c..., round(Int, 100 * frac[i]))
        end
    end
    say(hemi, " written: ", count(≥(0.95f0), frac), " of ", length(cells), " pixels ≥ 95 % glacier")
end
say("GLACIER_FRACTION_DONE in ", round(time() - t0), " s")
