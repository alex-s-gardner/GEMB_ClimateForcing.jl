using Test
using GEMB_ClimateForcing
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
                stack = climate_forcing(
                    :era5land, lat, lon;
                    time_range=time_range,
                    token=token
                )

                # Check type
                @test stack isa DimStack

                # Check all required variables are present
                @test haskey(stack, :temperature_air)
                @test haskey(stack, :pressure_air)
                @test haskey(stack, :precipitation)
                @test haskey(stack, :wind_speed)
                @test haskey(stack, :shortwave_downward)
                @test haskey(stack, :longwave_downward)
                @test haskey(stack, :vapor_pressure)

                # Check time dimension
                @test hasdim(stack, Ti)
                time_dim = dims(stack, Ti)
                n_steps = length(time_dim)
                @test n_steps >= 24  # At least 24 hourly steps for 1 day

                # Check dimensions match
                @test length(stack[:temperature_air]) == n_steps
                @test length(stack[:pressure_air]) == n_steps
                @test length(stack[:precipitation]) == n_steps
                @test length(stack[:wind_speed]) == n_steps
                @test length(stack[:shortwave_downward]) == n_steps
                @test length(stack[:longwave_downward]) == n_steps
                @test length(stack[:vapor_pressure]) == n_steps

                # Check physical ranges
                @test all(stack[:temperature_air] .> 200)  # > -73°C
                @test all(stack[:temperature_air] .< 320)  # < 47°C
                @test all(stack[:pressure_air] .> 0)
                @test all(stack[:pressure_air] .< 150000)
                @test all(stack[:precipitation] .>= 0)
                @test all(stack[:wind_speed] .>= 0)
                @test all(stack[:vapor_pressure] .>= 0)

                # Check metadata
                meta = metadata(stack)
                @test haskey(meta, "latitude")
                @test haskey(meta, "longitude")
                @test haskey(meta, "temperature_air_mean")
                @test haskey(meta, "wind_speed_mean")
                @test haskey(meta, "precipitation_mean")
                @test haskey(meta, "temperature_observation_height")
                @test haskey(meta, "wind_observation_height")
                @test meta["temperature_observation_height"] == 2.0
                @test meta["wind_observation_height"] == 10.0
            end

            @testset "Parallel Loading Consistency" begin
                # Load same data twice (relies on caching for speed)
                lat, lon = 72.58, -38.46
                time_range = (DateTime(2020, 1, 1), DateTime(2020, 1, 2))

                stack1 = climate_forcing(:era5land, lat, lon; time_range=time_range, token=token)
                stack2 = climate_forcing(:era5land, lat, lon; time_range=time_range, token=token)

                # Verify identical results (parallel loading should be deterministic)
                @test stack1[:temperature_air] == stack2[:temperature_air]
                @test stack1[:pressure_air] == stack2[:pressure_air]
                @test stack1[:wind_speed] == stack2[:wind_speed]
                @test stack1[:precipitation] == stack2[:precipitation]
                @test stack1[:shortwave_downward] == stack2[:shortwave_downward]
                @test stack1[:longwave_downward] == stack2[:longwave_downward]
                @test stack1[:vapor_pressure] == stack2[:vapor_pressure]

                println("  ✓ Parallel loading produces consistent results")
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
