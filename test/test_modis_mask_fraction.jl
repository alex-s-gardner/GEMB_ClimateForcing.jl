using Test
using GEMB_ClimateForcing
using Rasters
using ArchGDAL
using Rasters.Lookups: Points

const G = GEMB_ClimateForcing

@testset "modis_mask_fraction" begin
    # A synthetic categorical mask in MODIS sinusoidal metres around one cell: class 2 west of the cell's
    # centre, class 0 east of it, so the cell is half covered.
    cell = (17, 2, 100, 100)
    ulx, uly = G._modis_tile_origin(cell[1], cell[2])
    x0 = ulx + (cell[4] - 1) * G._MODIS_PIXEL_M; y0 = uly - (cell[3] - 1) * G._MODIS_PIXEL_M
    d = 10.0
    xs = (x0 - 2000):d:(x0 + 2500); ys = (y0 + 2000):-d:(y0 - 2500)
    xc = x0 + G._MODIS_PIXEL_M / 2
    vals = Int8[x < xc ? 2 : 0 for x in xs, y in ys]
    sinu = Rasters.ProjString(
        "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=$(G._MODIS_SPHERE_R) +units=m +no_defs")
    m = Raster(vals, (X(xs; sampling = Points()), Y(ys; sampling = Points())); crs = sinu)

    @test modis_mask_fraction([cell], m, (2, 3)) ≈ Float32[0.5] atol = 0.01
    @test modis_mask_fraction([cell], m, (0,)) ≈ Float32[0.5] atol = 0.01
    @test modis_mask_fraction([cell], m, (2, 3); factor = 20) ≈ Float32[0.5] atol = 0.01
    # The same cell located from its centre point.
    lat, lon = G._modis_cell_center(cell...)
    @test modis_mask_fraction([lat], [lon], m, (2, 3)) == modis_mask_fraction([cell], m, (2, 3))
    # A buffer widens the square symmetrically, so the share stays at one half here.
    @test modis_mask_fraction([cell], m, (2, 3); buffer = 500) ≈ Float32[0.5] atol = 0.02
    # A cell wholly inside class 2, and one off the raster (counts as not in `classes`).
    @test modis_mask_fraction([(17, 2, 100, 97)], m, (2,)) == Float32[1.0]
    @test modis_mask_fraction([(17, 2, 100, 200)], m, (2,)) == Float32[0.0]

    # The same mask with interval sampling whose lookup holds pixel edges, not centres.
    using Rasters.Lookups: Intervals, Start
    me = Raster(vals, (X(xs .- d / 2; sampling = Intervals(Start())), Y(ys .+ d / 2; sampling = Intervals(Start())));
                crs = sinu)
    @test modis_mask_fraction([cell], me, (2, 3)) == modis_mask_fraction([cell], m, (2, 3))

    @test_throws "the mask raster carries no CRS" modis_mask_fraction(
        [cell], Raster(vals, (X(xs; sampling = Points()), Y(ys; sampling = Points()))), (2,))
    @test_throws "factor must be >= 1" modis_mask_fraction([cell], m, (2,); factor = 0)
    @test_throws "buffer must be >= 0 m" modis_mask_fraction([cell], m, (2,); buffer = -1)
    @test_throws "lat and lon must have the same length" modis_mask_fraction([lat], [lon, lon], m, (2,))
end
