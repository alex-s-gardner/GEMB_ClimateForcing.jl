using Test
using GEMB_ClimateForcing
using GEMB
using Dates
using Statistics
using DimensionalData

println("="^70)
println("GEMB_ClimateForcing.jl Test Suite")
println("="^70)

@testset "GEMB_ClimateForcing.jl" begin

    @testset "Package Loading" begin
        @test isdefined(GEMB_ClimateForcing, :climate_forcing)
        println("  ✓ Package loaded successfully")
    end

    # Include unit validation tests
    include("test_unit_validation.jl")

    @testset "Input Validation" begin
        # Missing time_range
        @test_throws ArgumentError climate_forcing(:era5land, 72.0, -38.0)

        # Invalid time range (start >= end)
        @test_throws ArgumentError climate_forcing(
            :era5land, 72.0, -38.0;
            time_range=(DateTime(2020,12,31), DateTime(2020,1,1))
        )

        # Invalid latitude
        @test_throws ArgumentError climate_forcing(
            :era5land, 95.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )
        @test_throws ArgumentError climate_forcing(
            :era5land, -95.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )

        # Invalid longitude
        @test_throws ArgumentError climate_forcing(
            :era5land, 72.0, 400.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )

        # Invalid chunk_strategy
        @test_throws ArgumentError climate_forcing(
            :era5land, 72.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1)),
            chunk_strategy=:invalid
        )

        # Unsupported dataset
        @test_throws ArgumentError climate_forcing(
            :unsupported, 72.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )

        # Future datasets (not yet implemented)
        @test_throws ArgumentError climate_forcing(
            :era5, 72.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )
        @test_throws ArgumentError climate_forcing(
            :merra2, 72.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )
    end

    @testset "ERA5-Land Loading (Integration)" begin
        # Only run if CDS_API_KEY is available
        token = get(ENV, "CDS_API_KEY", nothing)

        if !isnothing(token)
            println("\n" * "="^70)
            println("Running ERA5-Land Data Retrieval Tests")
            println("(This will download actual data from ECMWF servers)")
            println("="^70)

            # Test with small time window
            lat, lon = 72.58, -38.46  # Summit, Greenland
            time_range = (DateTime(2020, 1, 1), DateTime(2020, 1, 2))  # 1 day

            @testset "Load Small Dataset" begin
                cf = climate_forcing(
                    :era5land, lat, lon;
                    time_range=time_range,
                    token=token
                )

                # Check type
                @test cf isa GEMB.ClimateForcing

                # Check all fields are present
                @test hasfield(typeof(cf), :temperature_air)
                @test hasfield(typeof(cf), :pressure_air)
                @test hasfield(typeof(cf), :precipitation)
                @test hasfield(typeof(cf), :wind_speed)
                @test hasfield(typeof(cf), :shortwave_downward)
                @test hasfield(typeof(cf), :longwave_downward)
                @test hasfield(typeof(cf), :vapor_pressure)

                # Check dimensions
                n_steps = length(cf.temperature_air)
                @test n_steps >= 24  # At least 24 hourly steps for 1 day
                @test length(cf.pressure_air) == n_steps
                @test length(cf.precipitation) == n_steps
                @test length(cf.wind_speed) == n_steps
                @test length(cf.shortwave_downward) == n_steps
                @test length(cf.longwave_downward) == n_steps
                @test length(cf.vapor_pressure) == n_steps

                # Check physical ranges
                @test all(cf.temperature_air .> 200)  # > -73°C
                @test all(cf.temperature_air .< 320)  # < 47°C
                @test all(cf.pressure_air .> 0)
                @test all(cf.pressure_air .< 150000)
                @test all(cf.precipitation .>= 0)
                @test all(cf.wind_speed .>= 0)
                @test all(cf.vapor_pressure .>= 0)

                # Check metadata
                @test cf.time_step > 0
                @test cf.temperature_air_mean > 200
                @test cf.wind_speed_mean >= 0
                @test cf.precipitation_mean >= 0
                @test cf.temperature_observation_height == 2.0
                @test cf.wind_observation_height == 10.0
            end

            @testset "Compatible with GEMB" begin
                # Verify ClimateForcing works with GEMB
                cf = climate_forcing(
                    :era5land, lat, lon;
                    time_range=time_range,
                    token=token
                )

                mp = GEMB.ModelParameters(output_frequency=:last)
                profile = GEMB.initialize_profile(mp, cf)

                # Profile can be either NamedTuple or DimStack depending on GEMB version
                @test profile isa Union{NamedTuple, DimStack}

                # Test with single timestep
                output = GEMB.gemb(profile, cf, mp)
                @test output isa DimStack
            end

        else
            println("\n" * "="^70)
            println("⚠️  Skipping ERA5-Land Data Retrieval Tests")
            println("="^70)
            println("CDS_API_KEY not set.")
            println("\nTo test actual data retrieval:")
            println("  1. Get a free API key from https://cds.climate.copernicus.eu/")
            println("  2. Set environment variable: export CDS_API_KEY='your-token-here'")
            println("  3. Run: julia --project=. -e 'using Pkg; Pkg.test()'")
            println("\nOr run standalone test:")
            println("  julia --project=. test/test_data_retrieval.jl")
            println("="^70)
        end
    end

    @testset "ERA5-Land Error Handling" begin
        # Missing token
        @test_throws ArgumentError climate_forcing(
            :era5land, 72.0, -38.0;
            time_range=(DateTime(2020,1,1), DateTime(2020,2,1))
        )

        # Invalid token (only test if CDS_API_KEY not set, to avoid rate limiting)
        if isnothing(get(ENV, "CDS_API_KEY", nothing))
            @test_throws Exception climate_forcing(
                :era5land, 72.0, -38.0;
                time_range=(DateTime(2020,1,1), DateTime(2020,1,2)),
                token="invalid_token_12345"
            )
        end
    end

end
