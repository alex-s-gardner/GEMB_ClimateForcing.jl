using Test
using GEMB_ClimateForcing

@testset "Vapor Pressure Conversions" begin
    @testset "vapor_pressure_to_relative_humidity" begin
        # Test example from MATLAB documentation
        vapor_pressure = 313.9  # Pa
        temperature_air = 265.3  # K
        rh = vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)

        @test rh ≈ 92.7913 atol=0.001

        # Test edge cases
        @test vapor_pressure_to_relative_humidity(0.0, 273.15) == 0.0
        @test vapor_pressure_to_relative_humidity(611.0, 273.15) ≈ 100.0 atol=1.0

        # Test clamping
        # Very high vapor pressure should clamp to 100%
        @test vapor_pressure_to_relative_humidity(10000.0, 273.15) == 100.0

        # Note: MATLAB version doesn't have explicit error checking,
        # so Julia version also doesn't throw on negative inputs

        # Test array broadcasting
        vp_array = [313.9, 611.0, 0.0]
        temp_array = [265.3, 273.15, 273.15]
        rh_array = vapor_pressure_to_relative_humidity(vp_array, temp_array)

        @test length(rh_array) == 3
        @test rh_array[1] ≈ 92.7913 atol=0.001
        @test rh_array[2] ≈ 100.0 atol=1.0
        @test rh_array[3] == 0.0
    end

    @testset "relative_humidity_to_vapor_pressure" begin
        # Test roundtrip conversion
        temperature_air = 265.3  # K
        rh_in = 92.7913  # %
        vp = relative_humidity_to_vapor_pressure(temperature_air, rh_in)

        @test vp ≈ 313.9 atol=0.1

        # Test at freezing point with 100% RH
        vp_freeze = relative_humidity_to_vapor_pressure(273.15, 100.0)
        @test vp_freeze ≈ 611.0 atol=1.0

        # Test at 0% RH
        vp_zero = relative_humidity_to_vapor_pressure(273.15, 0.0)
        @test vp_zero == 0.0

        # Test array broadcasting
        temp_array = [265.3, 273.15, 280.0]
        rh_array = [92.7913, 100.0, 50.0]
        vp_array = relative_humidity_to_vapor_pressure(temp_array, rh_array)

        @test length(vp_array) == 3
        @test vp_array[1] ≈ 313.9 atol=0.1
        @test vp_array[2] ≈ 611.0 atol=1.0

        # Test roundtrip: RH -> VP -> RH
        rh_roundtrip = vapor_pressure_to_relative_humidity(vp_array[1], temp_array[1])
        @test rh_roundtrip ≈ rh_array[1] atol=0.01
    end

    @testset "dewpoint_to_vapor_pressure" begin
        # Test at freezing point (dewpoint = temp means 100% RH)
        vp_freeze = dewpoint_to_vapor_pressure(273.15)
        @test vp_freeze ≈ 611.0 atol=1.0

        # Test typical values
        vp_1 = dewpoint_to_vapor_pressure(263.15)  # -10°C dewpoint
        @test vp_1 > 0.0
        @test vp_1 < vp_freeze

        # Test monotonicity: higher dewpoint -> higher vapor pressure
        vp_2 = dewpoint_to_vapor_pressure(268.15)
        @test vp_2 > vp_1

        # Note: MATLAB version has warning check but Julia version doesn't implement it
        # (not critical for functionality)

        # Test array broadcasting
        td_array = [263.15, 268.15, 273.15]
        vp_array = dewpoint_to_vapor_pressure(td_array)

        @test length(vp_array) == 3
        @test all(vp_array .> 0.0)
        @test vp_array[1] < vp_array[2] < vp_array[3]
        @test vp_array[3] ≈ 611.0 atol=1.0
    end

    @testset "Cross-validation: vapor pressure conversions" begin
        # Test internal consistency
        temperature_air = 260.0  # K
        rh_original = 75.0  # %

        # Convert RH -> VP -> RH
        vp = relative_humidity_to_vapor_pressure(temperature_air, rh_original)
        rh_recovered = vapor_pressure_to_relative_humidity(vp, temperature_air)

        @test rh_recovered ≈ rh_original atol=0.01

        # Test with multiple temperatures
        temps = [250.0, 260.0, 270.0, 273.15, 280.0]
        rhs = [50.0, 75.0, 90.0, 100.0, 80.0]

        for (T, RH) in zip(temps, rhs)
            vp = relative_humidity_to_vapor_pressure(T, RH)
            rh_check = vapor_pressure_to_relative_humidity(vp, T)
            @test rh_check ≈ RH atol=0.01
        end

        # Test relationship between dewpoint and vapor pressure
        # At 100% RH, dewpoint = air temperature
        temp = 270.0
        vp_at_100 = relative_humidity_to_vapor_pressure(temp, 100.0)
        vp_from_dewpoint = dewpoint_to_vapor_pressure(temp)

        # These should be very close (same formula, minor coefficient differences)
        @test vp_at_100 ≈ vp_from_dewpoint rtol=0.01
    end
end
