#!/usr/bin/env julia

"""
Test unit validation for climate forcing data
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData

# Import the validation function
using GEMB_ClimateForcing: validate_climate_forcing_units

# Helper to create DimStack for testing
function make_test_stack(;
    temperature_air = [273.0, 274.0, 275.0],
    pressure_air = [101300.0, 101300.0, 101300.0],
    precipitation = [0.0, 0.5, 1.0],
    wind_speed = [5.0, 6.0, 7.0],
    shortwave_downward = [200.0, 250.0, 300.0],
    longwave_downward = [300.0, 310.0, 320.0],
    vapor_pressure = [1000.0, 1100.0, 1200.0],
    metadata = Dict{String,Any}()
)
    time_dim = Ti([DateTime(2020,1,1,i) for i in 1:length(temperature_air)])
    return DimStack((
        temperature_air = DimArray(temperature_air, (time_dim,)),
        pressure_air = DimArray(pressure_air, (time_dim,)),
        precipitation = DimArray(precipitation, (time_dim,)),
        wind_speed = DimArray(wind_speed, (time_dim,)),
        shortwave_downward = DimArray(shortwave_downward, (time_dim,)),
        longwave_downward = DimArray(longwave_downward, (time_dim,)),
        vapor_pressure = DimArray(vapor_pressure, (time_dim,))
    ); metadata)
end

@testset "Unit Validation" begin
    # Valid data
    valid_stack = make_test_stack()

    @testset "Valid units pass" begin
        @test validate_climate_forcing_units(valid_stack) == true
    end

    @testset "Temperature out of range" begin
        # Temperature in Celsius instead of Kelvin
        bad_temp = make_test_stack(temperature_air=[20.0, 25.0, 30.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_temp)

        # Extreme temperatures
        bad_temp2 = make_test_stack(temperature_air=[100.0, 150.0, 200.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_temp2)
    end

    @testset "Pressure out of range" begin
        # Pressure too low (possibly in hPa instead of Pa)
        bad_pressure = make_test_stack(pressure_air=[1013.0, 1013.0, 1013.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_pressure)

        # Pressure too high
        bad_pressure2 = make_test_stack(pressure_air=[150000.0, 150000.0, 150000.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_pressure2)
    end

    @testset "Precipitation errors" begin
        # Negative precipitation
        bad_precip = make_test_stack(precipitation=[-0.1, 0.5, 1.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_precip)

        # Precipitation in meters instead of kg/m² (too large)
        bad_precip2 = make_test_stack(precipitation=[0.001, 0.002, 0.003] .* 1e6)
        @test_throws ArgumentError validate_climate_forcing_units(bad_precip2)
    end

    @testset "Wind speed errors" begin
        # Negative wind speed
        bad_wind = make_test_stack(wind_speed=[-5.0, 6.0, 7.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_wind)

        # Unrealistic wind speed
        bad_wind2 = make_test_stack(wind_speed=[150.0, 160.0, 170.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_wind2)
    end

    @testset "Shortwave radiation errors" begin
        # Negative shortwave (nighttime should be 0, not negative)
        bad_sw = make_test_stack(shortwave_downward=[-10.0, 200.0, 300.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_sw)

        # Shortwave in J/m² instead of W/m² (too large)
        bad_sw2 = make_test_stack(shortwave_downward=[200.0, 250.0, 300.0] .* 3600)
        @test_throws ArgumentError validate_climate_forcing_units(bad_sw2)
    end

    @testset "Longwave radiation errors" begin
        # The J/m² to W/m² conversion applied twice: ÷3600² leaves a value near zero.
        bad_lw = make_test_stack(longwave_downward=[0.1, 0.2, 0.3])
        @test_throws "conversion was applied twice" validate_climate_forcing_units(bad_lw)

        # Longwave in J/m² instead of W/m² (too large)
        bad_lw2 = make_test_stack(longwave_downward=[300.0, 310.0, 320.0] .* 3600)
        @test_throws "expected ≤ 500 W/m²" validate_climate_forcing_units(bad_lw2)
    end

    # LW = ε·σ·T⁴, so what counts as a plausible irradiance depends on the air temperature. A fixed
    # W/m² floor cannot express that, and gets both ends wrong.
    @testset "Longwave is bounded through emissivity, not irradiance" begin
        # A dry clear sky at 220 K radiating 49 W/m² is ε = 0.37 — ordinary for a polar or
        # high-altitude record, and exactly what a 50 W/m² floor rejected.
        cold_dry = make_test_stack(temperature_air=[220.0, 222.0, 224.0],
                                   longwave_downward=[49.0, 52.0, 55.0],
                                   vapor_pressure=[1.0, 2.0, 3.0])
        @test validate_climate_forcing_units(cold_dry) == true

        # The same irradiance at 290 K is ε = 0.15, below the dry clear-sky limit of the
        # parameterization: no atmosphere is that transparent. A 50 W/m² floor accepted this.
        warm_dark = make_test_stack(temperature_air=[288.0, 289.0, 290.0],
                                    longwave_downward=[60.0, 61.0, 62.0])
        @test_throws "bulk emissivity" validate_climate_forcing_units(warm_dark)

        # ε above 1 is real: a surface inversion with a cloud base warmer than the 2 m temperature.
        inversion = make_test_stack(temperature_air=[220.0, 221.0, 222.0],
                                    longwave_downward=[160.0, 162.0, 164.0])
        @test validate_climate_forcing_units(inversion) == true
    end

    # The ceiling detects metres read as kg/m², a factor of 1000. A deliberate rescaling is the
    # experiment, so re-testing its product against the source bound would reject valid runs.
    @testset "Precipitation ceiling applies to ingested data only" begin
        scaled = make_test_stack(precipitation=[0.0, 80.0, 120.0],
                                 metadata=Dict{String,Any}("precipitation_scaling" => 4.0))
        @test validate_climate_forcing_units(scaled) == true

        # Identity scaling means nothing was rescaled, so the ceiling still holds.
        unscaled = make_test_stack(precipitation=[0.0, 80.0, 120.0],
                                   metadata=Dict{String,Any}("precipitation_scaling" => 1.0))
        @test_throws "expected ≤ 100 kg/m²/hr" validate_climate_forcing_units(unscaled)

        # Non-negativity is an invariant and survives scaling.
        @test_throws "expected ≥ 0 kg/m²" validate_climate_forcing_units(
            make_test_stack(precipitation=[-0.1, 0.5, 1.0],
                            metadata=Dict{String,Any}("precipitation_scaling" => 4.0)))
    end

    @testset "Vapor pressure errors" begin
        # Negative vapor pressure
        bad_vp = make_test_stack(vapor_pressure=[-100.0, 1000.0, 1100.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_vp)

        # Vapor pressure too high (exceeds saturation)
        bad_vp2 = make_test_stack(vapor_pressure=[15000.0, 16000.0, 17000.0])
        @test_throws ArgumentError validate_climate_forcing_units(bad_vp2)
    end

    @testset "Boundary conditions" begin
        # Test exactly at boundaries (should pass)
        boundary_stack = make_test_stack(
            temperature_air = [180.0, 273.15, 330.0],
            pressure_air = [30000.0, 101300.0, 110000.0],
            precipitation = [0.0, 50.0, 100.0],
            wind_speed = [0.0, 10.0, 100.0],
            shortwave_downward = [0.0, 500.0, 1500.0],
            longwave_downward = [50.0, 300.0, 500.0],
            vapor_pressure = [0.0, 1000.0, 10000.0]
        )
        @test validate_climate_forcing_units(boundary_stack) == true
    end
end

println("\n✓ All unit validation tests passed!")
