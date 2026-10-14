#!/usr/bin/env julia

"""
Solar geometry and the shortwave layers derived from it, against published ephemeris values,
closed-form insolation, direct quadrature, and the Erbs et al. (1982) and Lacis & Hansen (1974)
relations. Offline.
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData
using Statistics

using GEMB_ClimateForcing: _partition_shortwave, CLOUD_OPTICAL_THICKNESS_MAX,
    validate_climate_forcing_units, _cell_forcing_at_interval

@testset "Solar geometry" begin
    @testset "declination, equation of time and distance against the ephemeris" begin
        # Solstices and equinoxes (UTC), obliquity 23.44°.
        @test solar_declination(DateTime(2020, 6, 20, 21, 44)) ≈ 23.44 atol = 0.05
        @test solar_declination(DateTime(2020, 12, 21, 10, 2)) ≈ -23.44 atol = 0.05
        @test solar_declination(DateTime(2020, 3, 20, 3, 50)) ≈ 0.0 atol = 0.01
        @test solar_declination(DateTime(2020, 9, 22, 13, 31)) ≈ 0.0 atol = 0.01
        @test solar_declination(DateTime(1979, 3, 21, 5, 22)) ≈ 0.0 atol = 0.01
        # Equation of time extremes: about −14.2 min on 11 February, +16.4 min on 3 November.
        @test equation_of_time(DateTime(2020, 2, 11, 12)) ≈ -14.2 atol = 0.5
        @test equation_of_time(DateTime(2020, 11, 3, 12)) ≈ 16.4 atol = 0.5
        # Perihelion 0.98329 au (5 January 2020), aphelion 1.01671 au (4 July 2020).
        @test earth_sun_distance_factor(DateTime(2020, 1, 5, 8)) ≈ 1 / 0.98329^2 atol = 1e-3
        @test earth_sun_distance_factor(DateTime(2020, 7, 4, 12)) ≈ 1 / 1.01671^2 atol = 1e-3
    end

    @testset "solar noon: cos θ = cos(φ − δ)" begin
        for (φ, λ, day) in ((46.5, 8.0, Date(2020, 6, 21)), (-33.9, 151.2, Date(2020, 1, 15)),
                            (61.0, -147.0, Date(2020, 9, 1)), (0.0, -78.5, Date(2020, 3, 20)))
            # Solar noon in UTC: 12:00 less 4 min per degree east and the equation of time,
            # iterated once so the equation of time is taken at noon itself.
            noon = DateTime(day) + Hour(12)
            for _ in 1:2
                noon = DateTime(day) + Hour(12) -
                       Millisecond(round(Int, (4λ + equation_of_time(noon)) * 60_000))
            end
            δ = solar_declination(noon)
            @test cos_solar_zenith(φ, λ, noon) ≈ cosd(φ - δ) atol = 1e-9
        end
    end

    @testset "hourly top-of-atmosphere sums to the daily insolation" begin
        # Q = (S₀/π)(d₀/d)²(h₀ sin φ sin δ + cos φ cos δ sin h₀), cos h₀ = −tan φ tan δ, for the day
        # centered on local solar noon; 24 hourly intervals cover one full rotation.
        function daily_mean(φ, t)
            δ = deg2rad(solar_declination(t))
            c = -tand(φ) * tan(δ)
            h₀ = c >= 1 ? 0.0 : c <= -1 ? Float64(π) : acos(c)
            return SOLAR_CONSTANT / π * earth_sun_distance_factor(t) *
                   (h₀ * sind(φ) * sin(δ) + cosd(φ) * cos(δ) * sin(h₀))
        end
        for (φ, day) in ((46.5, Date(2020, 6, 20)), (0.0, Date(2020, 3, 20)),
                         (-45.0, Date(2020, 12, 21)), (80.0, Date(2020, 6, 20)),
                         (70.0, Date(2020, 12, 21)), (30.0, Date(2020, 11, 5)))
            λ = 0.0
            times = DateTime(day) .+ Hour.(1:24)                    # intervals (00, 01] … (23, 24]
            toa = solar_geometry(φ, λ, times).shortwave_toa
            Q = daily_mean(φ, DateTime(day) + Hour(12))
            if Q == 0
                @test all(iszero, toa)                              # polar night
            else
                @test mean(toa) ≈ Q rtol = 1e-3
            end
        end
    end

    @testset "interval averages match direct quadrature" begin
        # Includes the hours containing sunrise and sunset, where the closed form has to clip the
        # integral at the horizon.
        φ, λ = 61.0, -147.0
        times = DateTime(2020, 9, 1) .+ Hour.(1:48)
        g = solar_geometry(φ, λ, times)
        for (i, t) in pairs(times)
            τs = [t - Hour(1) + Millisecond(180j) for j in 0:20_000]
            μ = [max(cos_solar_zenith(φ, λ, s), 0.0) for s in τs]
            w = fill(1.0, length(μ)); w[1] = w[end] = 0.5          # trapezoid weights
            I₁, I₂ = sum(w .* μ), sum(w .* μ .^ 2)
            toa = SOLAR_CONSTANT * earth_sun_distance_factor(t - Minute(30)) * I₁ / sum(w)
            # The closed form holds declination and the equation of time at the midpoint, which
            # moves a horizon crossing by seconds: a few hundredths of a W/m² in a sunrise hour.
            @test g.shortwave_toa[i] ≈ toa atol = 0.05
            if I₁ > 0
                @test cosd(g.solar_zenith_angle[i]) ≈ I₂ / I₁ atol = 1e-4
            else
                @test g.solar_zenith_angle[i] == 90.0
                @test g.shortwave_toa[i] == 0.0
            end
        end
        @test any(==(90.0), g.solar_zenith_angle)
        @test all(z -> 0 <= z <= 90, g.solar_zenith_angle)
        @test_throws ArgumentError solar_geometry(φ, λ, times; step = Hour(25))
    end

    @testset "Erbs, Klein & Duffie (1982) diffuse fraction" begin
        @test erbs_diffuse_fraction(0.0) == 1.0
        @test erbs_diffuse_fraction(0.22) ≈ 1 - 0.09 * 0.22
        @test erbs_diffuse_fraction(0.5) ≈ 0.9511 - 0.1604 * 0.5 + 4.388 * 0.25 - 16.638 * 0.125 +
                                           12.336 * 0.0625
        @test erbs_diffuse_fraction(0.81) == 0.165
        @test erbs_diffuse_fraction(1.0) == 0.165
        # The published pieces meet to within 1e-3 at both breakpoints.
        @test erbs_diffuse_fraction(nextfloat(0.22)) ≈ erbs_diffuse_fraction(0.22) atol = 1e-3
        @test erbs_diffuse_fraction(0.80) ≈ 0.165 atol = 1e-3
        @test_throws ArgumentError erbs_diffuse_fraction(-0.1)
    end

    @testset "cloud optical thickness inverts the cloud transmissivity" begin
        @test cloud_optical_thickness(1.0) == 0.0
        @test cloud_transmissivity(0.0) == 1.0
        for τ in (0.1, 1.0, 5.0, 20.0, 99.0)
            @test cloud_optical_thickness(cloud_transmissivity(τ)) ≈ τ rtol = 1e-12
        end
        @test cloud_transmissivity(10.0) ≈ 1 / 2.3
        @test cloud_optical_thickness(0.01) == CLOUD_OPTICAL_THICKNESS_MAX
        @test_throws ArgumentError cloud_optical_thickness(0.0)
        @test_throws ArgumentError cloud_optical_thickness(1.1)
    end

    @testset "partition of one step" begin
        z = 2000.0
        k_clear = 0.75 + 2e-5 * z
        # Night: all shortwave diffuse, τ = 0.
        @test _partition_shortwave(0.0, 0.0, z) == (0.0, 0.0)
        @test _partition_shortwave(0.4, 0.0, z) == (0.4, 0.0)
        # Clear: the clear-sky clearness or more gives τ = 0.
        d, τ = _partition_shortwave(k_clear * 1000, 1000.0, z)
        @test τ == 0.0
        @test d ≈ erbs_diffuse_fraction(k_clear) * k_clear * 1000
        # Shortwave above the geometric top of atmosphere is treated as k_t = 1.
        @test _partition_shortwave(120.0, 100.0, z) == (0.165 * 120.0, 0.0)
        # Overcast: τ from the transmissivity k_t / k_clear.
        d, τ = _partition_shortwave(200.0, 1000.0, z)
        @test τ ≈ cloud_optical_thickness(0.2 / k_clear)
        @test d ≈ erbs_diffuse_fraction(0.2) * 200.0
        # Daylight with no shortwave at all is the thickest cloud the inversion reports.
        @test _partition_shortwave(0.0, 500.0, z) == (0.0, CLOUD_OPTICAL_THICKNESS_MAX)
        @test all(isnan, _partition_shortwave(NaN, 500.0, z))
    end

    @testset "shortwave_forcing_layers" begin
        times = DateTime(2020, 7, 1) .+ Hour.(1:48)
        g = solar_geometry(46.5, 8.0, times)
        sw = 0.6 .* g.shortwave_toa
        L = shortwave_forcing_layers(sw, 46.5, 8.0, times, 2500.0)
        @test L.solar_zenith_angle == g.solar_zenith_angle
        @test all(0 .<= L.shortwave_downward_diffuse .<= sw)
        @test all(0 .<= L.cloud_optical_thickness .<= CLOUD_OPTICAL_THICKNESS_MAX)
        @test L.shortwave_downward_diffuse == first.(_partition_shortwave.(sw, g.shortwave_toa, 2500.0))
    end

    @testset "the layers ride through every adjustment unchanged" begin
        times = DateTime(2020, 7, 1) .+ Hour.(1:48)
        L = shortwave_forcing_layers(0.6 .* solar_geometry(46.5, 8.0, times).shortwave_toa,
                                     46.5, 8.0, times, 2500.0)
        n = length(times)
        ti = Ti(times)
        T = [273.15 + 5 * sinpi(2i / 24) for i in 1:n]
        sw = 0.6 .* solar_geometry(46.5, 8.0, times).shortwave_toa
        s = DimStack((temperature_air = DimArray(T, ti), pressure_air = DimArray(fill(75000.0, n), ti),
                      vapor_pressure = DimArray(fill(500.0, n), ti),
                      wind_speed = DimArray(fill(3.0, n), ti), precipitation = DimArray(fill(0.1, n), ti),
                      shortwave_downward = DimArray(sw, ti),
                      longwave_downward = DimArray(fill(290.0, n), ti),
                      solar_zenith_angle = DimArray(L.solar_zenith_angle, ti),
                      shortwave_downward_diffuse = DimArray(L.shortwave_downward_diffuse, ti),
                      cloud_optical_thickness = DimArray(L.cloud_optical_thickness, ti));
                     metadata = Dict{String,Any}("elevation" => 2500.0))
        @test validate_climate_forcing_units(s)
        k = fill(0.8, n)
        for adj in (climate_adjust_for_elevation(s, 400.0), temperature_adjust(s, 2.0),
                    precipitation_adjust(s, 1.3), climate_adjust_for_glacier(s, 0.8),
                    _cell_forcing_at_interval(s, 400.0, k, 5.0, 0.5, k))
            for v in (:shortwave_downward, :solar_zenith_angle, :shortwave_downward_diffuse,
                      :cloud_optical_thickness)
                @test parent(adj[v]) == parent(s[v])
            end
        end

        # Out-of-range derived layers fail validation.
        bad(name, v) = rebuild(s; data = merge(map(parent, NamedTuple(s)), (; name => v)))
        @test_throws "solar_zenith_angle" validate_climate_forcing_units(
            bad(:solar_zenith_angle, fill(95.0, n)))
        @test_throws "exceeds shortwave_downward" validate_climate_forcing_units(
            bad(:shortwave_downward_diffuse, sw .+ 1))
        @test_throws "cloud_optical_thickness" validate_climate_forcing_units(
            bad(:cloud_optical_thickness, fill(-1.0, n)))
    end

    @testset "simulated forcing carries the layers" begin
        ds = simulate_climate_forcing("test_1", 3)
        for v in (:solar_zenith_angle, :shortwave_downward_diffuse, :cloud_optical_thickness)
            @test haskey(ds, v)
        end
        sw = parent(ds[:shortwave_downward])
        @test all(parent(ds[:solar_zenith_angle])[sw .== 0] .== 90.0)
        @test all(0 .<= parent(ds[:shortwave_downward_diffuse]) .<= sw)
    end
end
