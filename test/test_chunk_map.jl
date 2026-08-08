@testset "climate_chunk_map" begin

    @testset "Exports & symbols" begin
        @test isdefined(GEMB_ClimateForcing, :climate_chunk_map)
    end

    @testset "Argument validation (offline)" begin
        # Missing token
        @test_throws ArgumentError climate_chunk_map(:era5land)

        # Unsupported dataset
        @test_throws ArgumentError climate_chunk_map(:not_a_dataset; token="x")

        # Invalid chunk_strategy (caught by climate_chunk_map before dispatch)
        @test_throws ArgumentError climate_chunk_map(:era5land; chunk_strategy=:bad, token="x")

        # Future datasets not yet implemented
        @test_throws ArgumentError climate_chunk_map(:era5;   token="x")
        @test_throws ArgumentError climate_chunk_map(:merra2; token="x")

        # Unknown variable_group (dispatched through to era5_land_chunk_map)
        @test_throws ArgumentError climate_chunk_map(
            :era5land; token="x", variable_group="not-a-group"
        )
    end

    @testset "Chunk ID formula (unit test, no network)" begin
        # Synthetic grid: 7 lons, 5 lats; chunks of size 3×2
        lon_chunk_size = 3
        lat_chunk_size = 2
        n_lon = 7
        n_lat = 5
        n_lat_chunks = cld(n_lat, lat_chunk_size)   # = 3
        n_lon_chunks = cld(n_lon, lon_chunk_size)   # = 3

        chunk_ids = Matrix{Int64}(undef, n_lon, n_lat)
        for j in 1:n_lat, i in 1:n_lon
            lon_cid = (i - 1) ÷ lon_chunk_size
            lat_cid = (j - 1) ÷ lat_chunk_size
            chunk_ids[i, j] = lon_cid * n_lat_chunks + lat_cid
        end

        # Cells within the same chunk share an ID
        @test chunk_ids[1,1] == chunk_ids[2,1] == chunk_ids[1,2] == chunk_ids[2,2]

        # Different lon chunks differ
        @test chunk_ids[1,1] != chunk_ids[4,1]

        # Different lat chunks differ
        @test chunk_ids[1,1] != chunk_ids[1,3]

        # Partial final chunks: col 7 (lon_cid=2), row 5 (lat_cid=2) → id = 2*3+2 = 8
        @test chunk_ids[7,5] == 8

        # All IDs are non-negative
        @test all(chunk_ids .>= 0)

        # Number of unique IDs equals the product of chunk counts
        @test length(unique(chunk_ids)) == n_lon_chunks * n_lat_chunks
    end

    # Integration tests — opt-in via GEMB_TEST_CHUNK_MAP=1 + valid CDS_API_KEY
    run_integration = get(ENV, "GEMB_TEST_CHUNK_MAP", "") in ("1", "true", "TRUE")
    token = run_integration ? get(ENV, "CDS_API_KEY", nothing) : nothing

    if run_integration && !isnothing(token) && !isempty(strip(token))
        @testset "ERA5-Land geo chunk map (integration)" begin
            raster = climate_chunk_map(:era5land; chunk_strategy=:geo, token=token)

            @test raster isa Rasters.Raster{Int64}
            @test hasdim(raster, Rasters.X)
            @test hasdim(raster, Rasters.Y)
            @test eltype(raster) == Int64
            @test all(raster .>= 0)

            # ERA5-Land global 0.1° grid: 3600 lon × 1801 lat
            @test size(raster, Rasters.X) == 3600
            @test size(raster, Rasters.Y) == 1801

            # Metadata present and consistent
            meta = metadata(raster)
            @test haskey(meta, "chunk_strategy")
            @test haskey(meta, "lon_chunk_size")
            @test haskey(meta, "lat_chunk_size")
            @test haskey(meta, "n_lon_chunks")
            @test haskey(meta, "n_lat_chunks")
            @test haskey(meta, "total_spatial_chunks")
            @test meta["chunk_strategy"] == "geo"
            @test meta["total_spatial_chunks"] == meta["n_lon_chunks"] * meta["n_lat_chunks"]
        end

        @testset "ERA5-Land geo vs time produce different chunk patterns" begin
            geo  = climate_chunk_map(:era5land; chunk_strategy=:geo,  token=token)
            time = climate_chunk_map(:era5land; chunk_strategy=:time, token=token)

            # Different strategies should yield different chunk shapes
            @test metadata(geo)["lon_chunk_size"] != metadata(time)["lon_chunk_size"] ||
                  metadata(geo)["lat_chunk_size"] != metadata(time)["lat_chunk_size"]

            # Chunk ID matrices should differ
            @test geo != time
        end
    else
        println("  (Skipping ERA5-Land chunk map integration tests — set GEMB_TEST_CHUNK_MAP=1 and CDS_API_KEY to enable)")
    end

end
