#!/usr/bin/env julia

"""
ERA5-Land ingest: unit conversions, nearest-cell selection, record bounds and validation guards.
The offline tests run on hand-built arrays. The live tests read the ARCO stores and run only when
`CDS_API_KEY` is set; `GEMB_TEST_ERA5_CACHE` names a writable Zarr cache directory for them
(default: a temporary directory).
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData
using Statistics

using GEMB_ClimateForcing: _era5_land_conversions, _nearest_longitude_index, find_nearest_index,
    _record_slice, _check_cached_record, validate_climate_forcing_units,
    SHORTWAVE_NEGATIVE_TOLERANCE, PRECIPITATION_NEGATIVE_TOLERANCE

@testset "ERA5-Land ingest" begin
    @testset "unit conversions" begin
        c = _era5_land_conversions(Float32[250, 260, 270], Float32[240, 250, 260],
                                   Float32[70000, 70000, 70000], [0.0, 1e-3, 2.5e-3],
                                   [3.0, -5.0, 0.0], [4.0, 12.0, 0.0],
                                   3600 .* [0.0, 100.0, 350.5], fill(3600 * 300.0, 3))
        @test c.shortwave_downward ≈ [0.0, 100.0, 350.5] rtol = 1e-12
        @test c.longwave_downward ≈ fill(300.0, 3) rtol = 1e-12
        @test c.precipitation ≈ [0.0, 1.0, 2.5] rtol = 1e-12
        @test c.wind_speed ≈ [5.0, 13.0, 0.0] rtol = 1e-12
        @test c.temperature_air == [250.0, 260.0, 270.0]
        @test c.vapor_pressure ≈ dewpoint_to_vapor_pressure([240.0, 250.0, 260.0])
        @test all(v -> v isa Vector{Float64}, values(c))

        # De-accumulation noise down to the tolerance is floored at zero; anything more negative
        # fails, naming the layer.
        noisy = _era5_land_conversions([260.0], [250.0], [7e4], [-0.5e-3 * 1e-3], [1.0], [0.0],
                                       [-3600 * 0.05], [3600 * 300.0])
        @test noisy.shortwave_downward == [0.0]
        @test noisy.precipitation == [0.0]
        @test_throws "shortwave_downward" _era5_land_conversions(
            [260.0], [250.0], [7e4], [0.0], [1.0], [0.0],
            [-3600 * 2SHORTWAVE_NEGATIVE_TOLERANCE], [3600 * 300.0])
        @test_throws "precipitation" _era5_land_conversions(
            [260.0], [250.0], [7e4], [-2PRECIPITATION_NEGATIVE_TOLERANCE * 1e-3], [1.0], [0.0],
            [0.0], [3600 * 300.0])
        # NaN passes through for the caller to judge.
        @test isnan(only(_era5_land_conversions([NaN], [NaN], [NaN], [NaN], [NaN], [NaN],
                                                [NaN], [NaN]).shortwave_downward))
    end

    @testset "nearest longitude across the seam" begin
        # The ARCO longitude axis, −179.9…180.0 ascending, with the stored round-off.
        lon = collect(range(-179.9, 180.0; length = 3600)) .+ 1e-11
        @test lon[_nearest_longitude_index(lon, 210.0)] ≈ -150.0 atol = 1e-9
        @test lon[_nearest_longitude_index(lon, -179.97)] ≈ 180.0 atol = 1e-9
        @test lon[_nearest_longitude_index(lon, 359.97)] ≈ 0.0 atol = 1e-9
        @test lon[_nearest_longitude_index(lon, -38.46)] ≈ -38.5 atol = 1e-9
        # A 0–359.9 axis, as the invariant grids store it.
        inv = collect(0.0:0.1:359.9)
        @test inv[_nearest_longitude_index(inv, -0.04)] == 0.0
        @test inv[_nearest_longitude_index(inv, -0.1)] ≈ 359.9
        # Latitude needs no wrap.
        lat = collect(range(-90.0, 90.0; length = 1801))
        @test lat[find_nearest_index(lat, 46.53)] ≈ 46.5
    end

    @testset "time range must lie within the record" begin
        epoch = DateTime(1970)
        hours = collect(0.0:23.0)                    # 1970-01-01T00 to T23
        @test _record_slice(hours, epoch, (DateTime(1970, 1, 1, 2), DateTime(1970, 1, 1, 5))) == 3:6
        @test _record_slice(hours, epoch, (epoch, DateTime(1970, 1, 1, 23))) == 1:24
        @test_throws "not within the ERA5-Land record" _record_slice(
            hours, epoch, (DateTime(1970, 1, 1, 2), DateTime(1970, 1, 2, 0)))
        @test_throws ArgumentError _record_slice(hours, epoch, (epoch - Hour(1), epoch + Hour(3)))
    end

    @testset "a cached record that cannot answer the request" begin
        epoch = DateTime(1970)
        hours = collect(0.0:23.0)
        arrays = ntuple(_ -> zeros(1, 1, 24), 8)
        within = (epoch, DateTime(1970, 1, 1, 23))
        @test _check_cached_record("/cache", hours, epoch, within, arrays) === nothing
        # Past the cached end: the error says how to refresh, and that chunks must go too.
        past = (epoch, DateTime(1970, 1, 2, 3))
        @test_throws "after the record end in the cached metadata" _check_cached_record(
            "/cache", hours, epoch, past, arrays)
        @test_throws "every chunk file" _check_cached_record("/cache", hours, epoch, past, arrays)
        # Groups cached at different times disagree on the record length.
        mixed = (ntuple(_ -> zeros(1, 1, 24), 7)..., zeros(1, 1, 20))
        @test_throws "different record lengths" _check_cached_record(
            "/cache", hours, epoch, within, mixed)
    end

    @testset "non-finite values fail validation" begin
        ti = Ti(DateTime(2020, 1, 1) .+ Hour.(0:2))
        good = (temperature_air = [270.0, 271.0, 272.0], pressure_air = fill(9e4, 3),
                vapor_pressure = fill(400.0, 3), wind_speed = fill(4.0, 3),
                precipitation = fill(0.1, 3), shortwave_downward = fill(100.0, 3),
                longwave_downward = fill(280.0, 3))
        stack(nt) = DimStack(map(v -> DimArray(v, ti), nt))
        @test validate_climate_forcing_units(stack(good))
        for name in keys(good), bad in (NaN, Inf)
            v = copy(good[name])
            v[2] = bad
            @test_throws ArgumentError validate_climate_forcing_units(stack(merge(good, (; name => v))))
            @test_throws "non-finite" validate_climate_forcing_units(stack(merge(good, (; name => v))))
        end
    end

    token = get(ENV, "CDS_API_KEY", "")
    if isempty(strip(token))
        println("⚠️  Skipping live ERA5-Land ingest tests (CDS_API_KEY not set)")
    else
        cache = get(ENV, "GEMB_TEST_ERA5_CACHE", mktempdir())
        @testset "live: hourly axis and a leap year" begin
            s = climate_forcing(:era5land, 46.5, 8.0; token, cache_path = cache,
                                time_range = (DateTime(2020, 1, 1), DateTime(2020, 12, 31, 23)))
            t = collect(dims(s, Ti))
            @test length(t) == 8784
            @test all(==(Hour(1)), diff(t))
            @test_throws ArgumentError climate_forcing(:era5land, 46.5, 8.0; token,
                cache_path = cache, time_range = (DateTime(1900, 1, 1), DateTime(1900, 1, 2)))
            # Past the end of the cached record: refused with refresh instructions.
            @test_throws "cached metadata" climate_forcing(:era5land, 46.5, 8.0; token,
                cache_path = cache, time_range = (DateTime(2020, 1, 1), now(UTC) + Day(30)))
        end

        @testset "live: a cell outside the land mask" begin
            s = climate_forcing(:era5land, 0.0, -25.0; token, cache_path = cache,
                                time_range = (DateTime(2020, 1, 1), DateTime(2020, 1, 2)))
            @test DimensionalData.metadata(s)["outside_land_mask"]
            @test all(isnan, s[:temperature_air])
        end
    end
end
