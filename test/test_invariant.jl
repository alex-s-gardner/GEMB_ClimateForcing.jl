#!/usr/bin/env julia

"""
Tests for `climate_model_invariant` — loading time-invariant climate-model
parameters (land-sea mask, geopotential, vegetation, …) as lazy Rasters.

Offline tests (registry structure, argument validation, exports) always run.
Integration tests actually download files from ECMWF and are opt-in: set
`GEMB_TEST_INVARIANT=1` to enable them (they need network and fetch ~50 MB
NetCDF files, cached to a temp dir).
"""

using Test
using GEMB_ClimateForcing
using Rasters
using DimensionalData

using GEMB_ClimateForcing: ERA5_LAND_INVARIANT_PARAMETERS,
    _INVARIANT_REGISTRY, _INVARIANT_BASE_URL, _default_invariant_cache

@testset "Climate model invariant" begin

    @testset "Exports & symbols" begin
        @test isdefined(GEMB_ClimateForcing, :climate_model_invariant)
        @test isdefined(GEMB_ClimateForcing, :ERA5_LAND_INVARIANT_PARAMETERS)
    end

    @testset "ERA5-Land parameter registry" begin
        params = ERA5_LAND_INVARIANT_PARAMETERS

        # The 10 documented invariant parameters (GRIB shortNames).
        expected = Set([:lsm, :z, :cl, :dl, :cvl, :cvh, :glm, :slt, :tvl, :tvh])
        @test Set(keys(params)) == expected

        # Every entry carries a NetCDF filename, an integer paramId, and a description.
        for (k, info) in params
            @test info.file isa String
            @test endswith(info.file, ".nc")
            @test info.paramId isa Int
            @test info.paramId > 0
            @test info.description isa String && !isempty(info.description)
        end

        # A couple of known paramIds (ECMWF GRIB) as a spot check.
        @test params[:lsm].paramId == 172
        @test params[:z].paramId == 129

        # Model is wired into the registry with a base URL.
        @test haskey(_INVARIANT_REGISTRY, :era5_land)
        @test _INVARIANT_REGISTRY[:era5_land] === params
        @test haskey(_INVARIANT_BASE_URL, :era5_land)
        @test startswith(_INVARIANT_BASE_URL[:era5_land], "https://")
    end

    @testset "Argument validation (offline)" begin
        # Unsupported model.
        @test_throws ArgumentError climate_model_invariant(model=:not_a_model)

        # Unknown parameter for a supported model.
        @test_throws ArgumentError climate_model_invariant(model=:era5_land, parameter=:not_a_param)

        # Default cache path is model-specific and under a temp location.
        cache = _default_invariant_cache(:era5_land)
        @test occursin("era5_land", cache)
        @test occursin("invariant", cache)
    end

    # ------------------------------------------------------------------
    # Integration: real downloads (opt-in, no CDS token required).
    # ------------------------------------------------------------------
    run_integration = get(ENV, "GEMB_TEST_INVARIANT", "") in ("1", "true", "TRUE")

    if run_integration
        println("\n" * "="^70)
        println("Running invariant-parameter download tests (network, ~50 MB/file)")
        println("="^70)

        # Isolated cache dir so the test is self-contained and cleanable.
        cache = mktempdir()

        @testset "Single parameter → lazy Raster (:lsm)" begin
            lsm = climate_model_invariant(parameter=:lsm, cache_path=cache)

            @test lsm isa Raster
            # 2-D spatial raster (time axis dropped), ERA5-Land 0.1° global grid.
            @test hasdim(lsm, X)
            @test hasdim(lsm, Y)
            @test !hasdim(lsm, Ti)
            @test size(lsm) == (3600, 1801)

            # Lazy: data is disk-backed until read/cropped.
            @test occursin("FileArray", string(typeof(parent(lsm))))

            # File landed in the cache.
            @test isfile(joinpath(cache, ERA5_LAND_INVARIANT_PARAMETERS[:lsm].file))

            # Crop to a box around Iceland (0–360 lon convention) and read.
            iceland = read(lsm[X = 335 .. 347, Y = 63 .. 67])
            @test size(iceland, X) > 0
            @test size(iceland, Y) > 0

            vals = collect(skipmissing(iceland))
            @test !isempty(vals)
            # Land-sea mask is a fraction in [0, 1] …
            @test all(0.0 .<= vals .<= 1.0)
            # … and Iceland's coastline gives genuinely fractional cells (not just 0/1).
            @test any(v -> 0.0 < v < 1.0, vals)
        end

        @testset "Caching (no re-download)" begin
            f = joinpath(cache, ERA5_LAND_INVARIANT_PARAMETERS[:lsm].file)
            @test isfile(f)
            mtime_before = mtime(f)
            lsm2 = climate_model_invariant(parameter=:lsm, cache_path=cache)
            @test lsm2 isa Raster
            # Cached file was reused, not rewritten.
            @test mtime(f) == mtime_before
        end

        @testset "Geopotential → orography (:z)" begin
            z = climate_model_invariant(parameter=:z, cache_path=cache)
            @test z isa Raster
            @test size(z) == (3600, 1801)

            # z is geopotential (m² s⁻²); orography (m) = z / g.
            oro = read((z ./ 9.80665)[X = 335 .. 347, Y = 63 .. 67])
            elev = collect(skipmissing(oro))
            @test !isempty(elev)
            # Iceland surface elevations: from ~sea level to well under 3000 m.
            @test minimum(elev) > -200.0
            @test maximum(elev) < 3000.0
        end

        @testset "All parameters → lazy RasterStack" begin
            st = climate_model_invariant(cache_path=cache)
            @test st isa RasterStack
            @test Set(keys(st)) == Set(keys(ERA5_LAND_INVARIANT_PARAMETERS))
            # Layers remain lazy (disk-backed) in the stack.
            @test occursin("FileArray", string(typeof(parent(st[:z]))))
        end

        rm(cache; recursive=true, force=true)
    else
        println("\n" * "="^70)
        println("⚠️  Skipping invariant-parameter download tests")
        println("Set GEMB_TEST_INVARIANT=1 to run them (network, ~50 MB per file).")
        println("="^70)
    end

end
