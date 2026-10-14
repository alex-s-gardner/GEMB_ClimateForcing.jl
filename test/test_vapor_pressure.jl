using Test
using GEMB_ClimateForcing
using GEMB_ClimateForcing: saturation_vapor_pressure

# IFS saturation vapor pressure over water (Cy45r1 Part IV Eq. 7.5, Buck 1981 constants), coded
# here independently of the package.
_ifs_esw(T) = 611.21 * exp(17.502 * (T - 273.16) / (T - 32.19))

@testset "Vapor Pressure Conversions" begin
    @testset "saturation_vapor_pressure is the IFS curve over water" begin
        for T in 223.15:2.5:313.15
            @test saturation_vapor_pressure(T) ≈ _ifs_esw(T) rtol = 1e-14
        end
        # The triple point is the curve's reference: exactly the leading constant.
        @test saturation_vapor_pressure(273.16) == 611.21
    end

    @testset "dewpoint against IFS reference values" begin
        # Td = −40, −20, 0, 20 °C, values of the IFS water curve worked by hand.
        for (Td, e) in ((233.15, 18.74), (253.15, 125.27), (273.15, 610.77), (293.15, 2335.84))
            @test dewpoint_to_vapor_pressure(Td) ≈ e rtol = 2e-3
        end
        # Broadcasts and promotes Float32 to Float64.
        v = dewpoint_to_vapor_pressure(Float32[263.15, 268.15, 273.15])
        @test v isa Vector{Float64}
        @test v[1] < v[2] < v[3]
    end

    @testset "relative humidity is defined against the dewpoint's own curve" begin
        # Dewpoint equal to the air temperature is saturation: RH = 1 with no clamp involved.
        for T in 233.15:1.0:303.15
            @test dewpoint_to_vapor_pressure(T) / saturation_vapor_pressure(T) ≈ 1 atol = 1e-3
            @test vapor_pressure_to_relative_humidity(dewpoint_to_vapor_pressure(T), T) ≈ 100 atol = 1e-9
        end
    end

    @testset "vapor_pressure_to_relative_humidity" begin
        @test vapor_pressure_to_relative_humidity(313.9, 265.3) ≈ 100 * 313.9 / _ifs_esw(265.3)
        @test vapor_pressure_to_relative_humidity(0.0, 273.15) == 0.0
        # Clamped to [0, 100].
        @test vapor_pressure_to_relative_humidity(10000.0, 273.15) == 100.0
        rh = vapor_pressure_to_relative_humidity([313.9, 0.0], [265.3, 273.15])
        @test rh ≈ [100 * 313.9 / _ifs_esw(265.3), 0.0]
    end

    @testset "relative_humidity_to_vapor_pressure" begin
        @test relative_humidity_to_vapor_pressure(265.3, 80.0) ≈ 0.8 * _ifs_esw(265.3)
        @test relative_humidity_to_vapor_pressure(273.15, 0.0) == 0.0
        @test relative_humidity_to_vapor_pressure([265.3, 280.0], [80.0, 50.0]) ≈
              [0.8 * _ifs_esw(265.3), 0.5 * _ifs_esw(280.0)]
        # Round trip below saturation.
        for (T, RH) in zip([250.0, 260.0, 270.0, 273.15, 280.0], [50.0, 75.0, 90.0, 100.0, 80.0])
            @test vapor_pressure_to_relative_humidity(relative_humidity_to_vapor_pressure(T, RH), T) ≈ RH
        end
        # At 100 % the vapor pressure is the dewpoint conversion at the air temperature.
        @test relative_humidity_to_vapor_pressure(270.0, 100.0) ≈ dewpoint_to_vapor_pressure(270.0)
    end
end
