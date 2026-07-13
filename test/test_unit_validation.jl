#!/usr/bin/env julia

"""
Test unit validation for climate forcing data
"""

using Test
using GEMB_ClimateForcing
using Dates

# Import the validation function
using GEMB_ClimateForcing: validate_climate_forcing_units

@testset "Unit Validation" begin
    # Valid data
    valid_data = (
        temperature_air = [273.0, 274.0, 275.0],
        pressure_air = [101300.0, 101300.0, 101300.0],
        precipitation = [0.0, 0.5, 1.0],
        wind_speed = [5.0, 6.0, 7.0],
        shortwave_downward = [200.0, 250.0, 300.0],
        longwave_downward = [300.0, 310.0, 320.0],
        vapor_pressure = [1000.0, 1100.0, 1200.0]
    )

    @testset "Valid units pass" begin
        @test validate_climate_forcing_units(;valid_data...) == true
    end

    @testset "Temperature out of range" begin
        # Temperature in Celsius instead of Kelvin
        bad_temp = (valid_data..., temperature_air=[20.0, 25.0, 30.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_temp...)

        # Extreme temperatures
        bad_temp2 = (valid_data..., temperature_air=[100.0, 150.0, 200.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_temp2...)
    end

    @testset "Pressure out of range" begin
        # Pressure too low (possibly in hPa instead of Pa)
        bad_pressure = (valid_data..., pressure_air=[1013.0, 1013.0, 1013.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_pressure...)

        # Pressure too high
        bad_pressure2 = (valid_data..., pressure_air=[150000.0, 150000.0, 150000.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_pressure2...)
    end

    @testset "Precipitation errors" begin
        # Negative precipitation
        bad_precip = (valid_data..., precipitation=[-0.1, 0.5, 1.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_precip...)

        # Precipitation in meters instead of kg/m² (too large)
        bad_precip2 = (valid_data..., precipitation=[0.001, 0.002, 0.003] .* 1e6)
        @test_throws ArgumentError validate_climate_forcing_units(;bad_precip2...)
    end

    @testset "Wind speed errors" begin
        # Negative wind speed
        bad_wind = (valid_data..., wind_speed=[-5.0, 6.0, 7.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_wind...)

        # Unrealistic wind speed
        bad_wind2 = (valid_data..., wind_speed=[150.0, 160.0, 170.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_wind2...)
    end

    @testset "Shortwave radiation errors" begin
        # Negative shortwave (nighttime should be 0, not negative)
        bad_sw = (valid_data..., shortwave_downward=[-10.0, 200.0, 300.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_sw...)

        # Shortwave in J/m² instead of W/m² (too large)
        bad_sw2 = (valid_data..., shortwave_downward=[200.0, 250.0, 300.0] .* 3600)
        @test_throws ArgumentError validate_climate_forcing_units(;bad_sw2...)
    end

    @testset "Longwave radiation errors" begin
        # Longwave too low (possibly in J/m² instead of W/m²)
        bad_lw = (valid_data..., longwave_downward=[0.1, 0.2, 0.3])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_lw...)

        # Longwave in J/m² instead of W/m² (too large)
        bad_lw2 = (valid_data..., longwave_downward=[300.0, 310.0, 320.0] .* 3600)
        @test_throws ArgumentError validate_climate_forcing_units(;bad_lw2...)
    end

    @testset "Vapor pressure errors" begin
        # Negative vapor pressure
        bad_vp = (valid_data..., vapor_pressure=[-100.0, 1000.0, 1100.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_vp...)

        # Vapor pressure too high (exceeds saturation)
        bad_vp2 = (valid_data..., vapor_pressure=[15000.0, 16000.0, 17000.0])
        @test_throws ArgumentError validate_climate_forcing_units(;bad_vp2...)
    end

    @testset "Boundary conditions" begin
        # Test exactly at boundaries (should pass)
        boundary_data = (
            temperature_air = [180.0, 273.15, 330.0],
            pressure_air = [30000.0, 101300.0, 110000.0],
            precipitation = [0.0, 50.0, 100.0],
            wind_speed = [0.0, 10.0, 100.0],
            shortwave_downward = [0.0, 500.0, 1500.0],
            longwave_downward = [50.0, 300.0, 500.0],
            vapor_pressure = [0.0, 1000.0, 10000.0]
        )
        @test validate_climate_forcing_units(;boundary_data...) == true
    end
end

println("\n✓ All unit validation tests passed!")
