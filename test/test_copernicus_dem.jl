#!/usr/bin/env julia

"""
Tests for the Copernicus GLO-30 (30 m) DEM reader, reached through
`climate_model_invariant(model=:copernicus_dem_30m, extent=…)`.

Offline tests (tile-ID formatting, extent→tile enumeration, argument validation,
model registration) always run. Integration tests actually read tiles over HTTP
from the AWS Open Data bucket and are opt-in: set `GEMB_TEST_COPERNICUS_DEM=1` to
enable them (they need network; only small windows are fetched via GDAL /vsicurl/,
not whole tiles).
"""

using Test
using GEMB_ClimateForcing
using Rasters
using DimensionalData

using GEMB_ClimateForcing: _copernicus_dem_tile_id, _copernicus_dem_tile_url,
    _copernicus_dem_vsicurl, _copernicus_dem_extent, _copernicus_dem_tiles_for_extent,
    _copernicus_dem_ncols, _copernicus_dem_parse_corner, _copernicus_dem_write_vrt,
    _INVARIANT_EXTENT_MODELS, _COPERNICUS_DEM_30M_BASE

@testset "Copernicus GLO-30 DEM" begin

    @testset "Model registration" begin
        @test :copernicus_dem_30m in _INVARIANT_EXTENT_MODELS
        @test startswith(_COPERNICUS_DEM_30M_BASE, "https://")
    end

    @testset "Tile-ID formatting" begin
        # SW corner encoding: N/S 2-digit lat, E/W 3-digit lon, each + "_00".
        @test _copernicus_dem_tile_id(59, -150) == "Copernicus_DSM_COG_10_N59_00_W150_00_DEM"
        @test _copernicus_dem_tile_id(0, 6)     == "Copernicus_DSM_COG_10_N00_00_E006_00_DEM"
        @test _copernicus_dem_tile_id(-34, 18)  == "Copernicus_DSM_COG_10_S34_00_E018_00_DEM"
        @test _copernicus_dem_tile_id(-1, -1)   == "Copernicus_DSM_COG_10_S01_00_W001_00_DEM"

        id = _copernicus_dem_tile_id(59, -150)
        @test _copernicus_dem_tile_url(id) ==
            "$(_COPERNICUS_DEM_30M_BASE)/$(id)/$(id).tif"
        @test _copernicus_dem_vsicurl(id) == "/vsicurl/" * _copernicus_dem_tile_url(id)

        # Round-trip: parse the SW corner back out of a tile id.
        @test _copernicus_dem_parse_corner(id) == (59, -150)
        @test _copernicus_dem_parse_corner("Copernicus_DSM_COG_10_S34_00_E018_00_DEM") == (-34, 18)
    end

    @testset "Latitude-banded column counts" begin
        # Longitude is sub-sampled poleward; rows are always 3600. Verified against
        # real tile headers across every band.
        @test _copernicus_dem_ncols(0)   == 3600
        @test _copernicus_dem_ncols(49)  == 3600
        @test _copernicus_dem_ncols(50)  == 2400
        @test _copernicus_dem_ncols(59)  == 2400
        @test _copernicus_dem_ncols(65)  == 1800
        @test _copernicus_dem_ncols(72)  == 1200
        @test _copernicus_dem_ncols(82)  == 720
        @test _copernicus_dem_ncols(87)  == 360
        # Southern hemisphere: equatorward edge is |lat_sw|-1.
        @test _copernicus_dem_ncols(-50) == 3600   # S50 tile spans [-50,-49): band <50
        @test _copernicus_dem_ncols(-51) == 2400
    end

    @testset "Analytical VRT authoring (offline)" begin
        # Build a VRT across the 50°N band boundary (N49=3600 cols, N50=2400 cols) with
        # no network access, and check the XML geometry is internally consistent.
        ids = ["Copernicus_DSM_COG_10_N49_00_W123_00_DEM",
               "Copernicus_DSM_COG_10_N50_00_W123_00_DEM"]
        path = tempname() * ".vrt"
        _copernicus_dem_write_vrt(ids, path)
        xml = read(path, String)
        # Common finest grid: 1° wide × 2° tall at 3600/° → 3600 × 7200.
        @test occursin("""rasterXSize="3600" rasterYSize="7200\"""", xml)
        @test occursin("<GeoTransform>-123.0, ", xml)          # west edge
        # Each source keeps its native column count in SourceProperties/SrcRect …
        @test occursin("""RasterXSize="3600\"""", xml)
        @test occursin("""RasterXSize="2400\"""", xml)
        # … but is placed on the 3600-col destination grid.
        @test occursin("""<DstRect xOff="0" yOff="0" xSize="3600" ySize="3600"/>""", xml)
        @test occursin("/vsicurl/", xml)
        rm(path; force=true)
    end

    @testset "Extent normalization" begin
        # Extents.Extent (re-exported by Rasters) and NamedTuple both accepted.
        @test _copernicus_dem_extent(Extent(X=(-150.0, -149.0), Y=(59.0, 60.0))) ==
            (-150.0, -149.0, 59.0, 60.0)
        @test _copernicus_dem_extent((X=(-150.0, -149.0), Y=(59.0, 60.0))) ==
            (-150.0, -149.0, 59.0, 60.0)
        # Out-of-range boxes are rejected.
        @test_throws ArgumentError _copernicus_dem_extent((X=(-181.0, -179.0), Y=(0.0, 1.0)))
        @test_throws ArgumentError _copernicus_dem_extent((X=(0.0, 1.0), Y=(80.0, 95.0)))
    end

    @testset "Extent → covering tiles" begin
        idx = Set([
            "Copernicus_DSM_COG_10_N59_00_W150_00_DEM",
            "Copernicus_DSM_COG_10_N59_00_W149_00_DEM",
        ])
        # A box inside [-150,-149) resolves to the W150 tile (SW-corner anchored).
        @test _copernicus_dem_tiles_for_extent(-149.6, -149.4, 59.2, 59.4, idx) ==
            ["Copernicus_DSM_COG_10_N59_00_W150_00_DEM"]
        # A box spanning the boundary picks up both tiles.
        @test Set(_copernicus_dem_tiles_for_extent(-150.5, -148.5, 59.2, 59.4, idx)) == idx
        # An all-ocean box (no tiles in the index) errors clearly.
        @test_throws ArgumentError _copernicus_dem_tiles_for_extent(-10.5, -10.4, 20.2, 20.4, idx)
    end

    @testset "Argument validation (offline)" begin
        # Unsupported model.
        @test_throws ArgumentError climate_model_invariant(model=:not_a_model,
                                        extent=Extent(X=(0,1), Y=(0,1)))
        # DEM is extent-based and does not take a `parameter`.
        @test_throws ArgumentError climate_model_invariant(model=:copernicus_dem_30m,
                                        parameter=:z, extent=Extent(X=(0,1), Y=(0,1)))
    end

    # ------------------------------------------------------------------
    # Integration: real /vsicurl/ reads (opt-in, no token required).
    # ------------------------------------------------------------------
    run_integration = get(ENV, "GEMB_TEST_COPERNICUS_DEM", "") in ("1", "true", "TRUE")

    if run_integration
        println("\n" * "="^70)
        println("Running Copernicus DEM /vsicurl/ tests (network; small windows only)")
        println("="^70)

        cache = mktempdir()

        @testset "Single tile → lazy Raster" begin
            # A land box inside one tile in the Kenai Peninsula, Alaska.
            dem = climate_model_invariant(model=:copernicus_dem_30m,
                      extent=Extent(X=(-149.9, -149.7), Y=(59.6, 59.8)),
                      cache_path=cache)

            @test dem isa Rasters.AbstractRaster
            @test hasdim(dem, X)
            @test hasdim(dem, Y)
            @test metadata(dem)["n_tiles"] == 1
            @test metadata(dem)["tiles"] == ["Copernicus_DSM_COG_10_N59_00_W150_00_DEM"]
            @test metadata(dem)["resolution_m"] == 30

            # Lazy until read: the parent is a disk/DiskArray view, not a dense Array.
            @test !(parent(dem) isa Array)

            # Tile index landed in the cache.
            @test isfile(joinpath(cache, "tileList.txt"))

            elev = read(dem)
            vals = collect(skipmissing(elev))
            @test !isempty(vals)
            # Plausible Kenai mountain elevations: sea level to a few thousand metres.
            @test minimum(vals) > -50.0
            @test maximum(vals) > 200.0     # genuine relief, not a flat/ocean patch
            @test maximum(vals) < 4000.0
        end

        @testset "Tile index cached (no re-download)" begin
            f = joinpath(cache, "tileList.txt")
            @test isfile(f)
            mtime_before = mtime(f)
            _ = climate_model_invariant(model=:copernicus_dem_30m,
                    extent=Extent(X=(-149.9, -149.7), Y=(59.6, 59.8)), cache_path=cache)
            @test mtime(f) == mtime_before
        end

        @testset "Multi-tile mosaic → single lazy Raster (analytical VRT)" begin
            dem = climate_model_invariant(model=:copernicus_dem_30m,
                      extent=Extent(X=(-150.2, -148.8), Y=(59.5, 59.9)),
                      cache_path=cache)
            @test dem isa Rasters.AbstractRaster
            @test metadata(dem)["n_tiles"] >= 2
            @test !(parent(dem) isa Array)

            vals = collect(skipmissing(read(dem)))
            @test !isempty(vals)
            @test maximum(vals) > 200.0
        end

        @testset "Mosaic values match direct tile reads" begin
            # Cross-check a window near the 50°N band boundary: the analytical VRT must
            # reproduce a direct single-tile read exactly (same coverage/values).
            box = Extent(X=(-122.9, -122.7), Y=(50.3, 50.5))
            viavrt = climate_model_invariant(model=:copernicus_dem_30m,
                         extent=Extent(X=(-123.0, -121.0), Y=(50.0, 51.0)), cache_path=cache)
            single = climate_model_invariant(model=:copernicus_dem_30m,
                         extent=box, cache_path=cache)   # a single N50 tile
            v_all = collect(skipmissing(read(viavrt[X=(-122.9 .. -122.7), Y=(50.3 .. 50.5)])))
            v_one = collect(skipmissing(read(single)))
            @test !isempty(v_one)
            @test extrema(v_all) == extrema(v_one)
        end

        rm(cache; recursive=true, force=true)
    else
        println("\n" * "="^70)
        println("⚠️  Skipping Copernicus DEM /vsicurl/ tests")
        println("Set GEMB_TEST_COPERNICUS_DEM=1 to run them (network).")
        println("="^70)
    end

end
