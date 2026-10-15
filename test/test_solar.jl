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

using GEMB_ClimateForcing: solar_declination, equation_of_time, earth_sun_distance_factor,
    cos_solar_zenith, erbs_diffuse_fraction, cloud_transmissivity, cloud_optical_thickness,
    _clear_sky_clearness, _inverted_optical_thickness, _partition_shortwave_series,
    CLOUD_OPTICAL_THICKNESS_MAX, CLOUD_OPTICAL_THICKNESS_ZENITH_MAX,
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
        for τ in (0.1, 1.0, 5.0, 20.0, 29.0)
            @test cloud_optical_thickness(cloud_transmissivity(τ)) ≈ τ rtol = 1e-12
        end
        @test cloud_transmissivity(10.0) ≈ 1 / 2.3
        # The cap is the top of the Gardner & Sharp (2010) fitted range, Table B1: τ = 0–30.
        @test CLOUD_OPTICAL_THICKNESS_MAX == 30.0
        # t = 1/(1 + 0.13·40) inverts to 40 and is held at the cap; t = 0.01 to 761.5.
        @test cloud_optical_thickness(1 / (1 + 0.13 * 40)) == 30.0
        @test cloud_optical_thickness(0.01) == 30.0
        @test cloud_optical_thickness(1 / (1 + 0.13 * 30)) ≈ 30.0 rtol = 1e-12
        @test_throws ArgumentError cloud_optical_thickness(0.0)
        @test_throws ArgumentError cloud_optical_thickness(1.1)
    end

    # FAO-56 Annex 3 eqs. 3-17 to 3-20 with K_t = 1, written out here independently of the
    # package: sin φ the sine of the sun's elevation, P and e_a in kPa.
    function fao_clearness(sin_φ, P, e_a)
        W = 0.14 * e_a * P + 2.1                                             # (3-19), mm
        K_B = 0.98 * exp(-0.00146 * P / sin_φ - 0.091 * (W / sin_φ)^0.25)   # (3-18)
        K_D = K_B >= 0.15 ? 0.35 - 0.33 * K_B : 0.18 + 0.82 * K_B           # (3-20)
        return K_B + K_D                                                     # (3-17)
    end

    @testset "FAO-56 Annex 3 clear-sky clearness" begin
        # Reference values from eqs. 3-17 to 3-20 evaluated in 256-bit BigFloat with the
        # coefficients parsed as decimals, at (sin φ, P kPa, e_a kPa):
        #   (0.5, 101.3, 1.0):  K_B = 0.58662, K_D = 0.35 − 0.33 K_B = 0.15642
        #   (0.9,  77.0, 0.4):  K_B = 0.74544, K_D = 0.35 − 0.33 K_B = 0.10401
        #   (0.08, 101.3, 2.0): K_B = 0.10322 < 0.15, K_D = 0.18 + 0.82 K_B = 0.26464
        # The package takes pressures in Pa.
        @test _clear_sky_clearness(0.5, 101.3e3, 1.0e3) ≈ 0.7430339046365992 rtol = 1e-12
        @test _clear_sky_clearness(0.9, 77.0e3, 0.4e3) ≈ 0.8494422160370874 rtol = 1e-12
        @test _clear_sky_clearness(0.08, 101.3e3, 2.0e3) ≈ 0.36785562641330877 rtol = 1e-12
        @test_throws "cosine of the solar zenith angle" _clear_sky_clearness(0.0, 1e5, 1e3)
        @test_throws "cosine of the solar zenith angle" _clear_sky_clearness(1.1, 1e5, 1e3)
        @test_throws "surface pressure must be positive" _clear_sky_clearness(0.5, 0.0, 1e3)
        @test_throws "surface pressure must be positive" _clear_sky_clearness(0.5, NaN, 1e3)
        @test_throws "vapor pressure must be non-negative" _clear_sky_clearness(0.5, 1e5, -1.0)
    end

    @testset "inversion of one sunlit step" begin
        P, e = 75.0e3, 500.0                                # Pa
        toa = 1000.0
        for θ in (30.0, 60.0, 75.0)
            k_clear = fao_clearness(cosd(θ), P / 1000, e / 1000)
            # Shortwave equal to the clear sky inverts to a clear sky.
            @test _inverted_optical_thickness(k_clear * toa, toa, cosd(θ), P, e) ≈ 0 atol = 1e-12
            # Transmissivity 1/(1 + 0.13·10) inverts to τ = 10.
            @test _inverted_optical_thickness(k_clear * toa / (1 + 0.13 * 10), toa, cosd(θ), P, e) ≈
                  10.0 rtol = 1e-12
            # Transmissivity 1/(1 + 0.13·50) would invert to 50 and is held at the cap.
            @test _inverted_optical_thickness(k_clear * toa / (1 + 0.13 * 50), toa, cosd(θ), P, e) ==
                  CLOUD_OPTICAL_THICKNESS_MAX
        end
        # Shortwave above the clear sky, or above the top of atmosphere, is clear.
        @test _inverted_optical_thickness(900.0, toa, cosd(30.0), P, e) == 0.0
        @test _inverted_optical_thickness(1200.0, toa, cosd(30.0), P, e) == 0.0
        # Daylight with no shortwave at all is the thickest cloud the inversion reports.
        @test _inverted_optical_thickness(0.0, toa, cosd(30.0), P, e) == CLOUD_OPTICAL_THICKNESS_MAX
        @test isnan(_inverted_optical_thickness(NaN, toa, cosd(30.0), P, e))
        @test_throws "only on a sunlit step" _inverted_optical_thickness(0.0, 0.0, 0.5, P, e)
    end

    @testset "partition of a series: diffuse and the low-sun rule" begin
        P, e = 70.0e3, 300.0                                 # Pa
        # One synthetic series: night, a day, night, a twilight-only period, night. `target` is
        # the τ the shortwave of each inverted step is built to give; low-sun steps get no
        # shortwave at all, which an inversion would report as the cap.
        sza    = [90.0, 88.0, 82.0, 80.0, 60.0, 40.0, 79.0, 84.0, 89.0, 90.0, 86.0, 85.0, 87.0, 90.0]
        target = [ NaN,  NaN,  NaN,  2.0,  3.0,  4.0,  6.0,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN]
        toa = [z < 90 ? SOLAR_CONSTANT * cosd(z) : 0.0 for z in sza]
        sw = [isnan(τ) ? 0.0 : fao_clearness(cosd(z), P / 1000, e / 1000) * a / (1 + 0.13τ)
              for (z, a, τ) in zip(sza, toa, target)]
        n = length(sza)
        diffuse, τ = _partition_shortwave_series(sw, toa, sza, fill(P, n), fill(e, n))
        # Dawn (88°, 82°) takes the first inverted step (80°, exactly at the threshold); dusk
        # (84°, 89°) the last (79°); night and the twilight-only period are 0.
        expected = [0.0, 2.0, 2.0, 2.0, 3.0, 4.0, 6.0, 6.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        @test CLOUD_OPTICAL_THICKNESS_ZENITH_MAX == 80.0
        @test all(i -> isapprox(τ[i], expected[i]; rtol = 1e-12, atol = 0), eachindex(τ, expected))
        # Diffuse is f(k_t)·SW at every sunlit step, the Erbs fraction of k_t = SW/TOA; at night
        # all the shortwave is diffuse.
        for i in eachindex(sw, toa)
            @test diffuse[i] == (toa[i] > 0 ? erbs_diffuse_fraction(sw[i] / toa[i]) * sw[i] : sw[i])
        end

        # Polar day: one daylight period spanning two "days". Low-sun steps take the latest
        # earlier inversion, the following one only before the first.
        sza = [85.0, 70.0, 85.0, 86.0, 70.0, 85.0]
        target = [NaN, 1.0, NaN, NaN, 5.0, NaN]
        toa = SOLAR_CONSTANT .* cosd.(sza)
        sw = [isnan(τ) ? 1.0 : fao_clearness(cosd(z), P / 1000, e / 1000) * a / (1 + 0.13τ)
              for (z, a, τ) in zip(sza, toa, target)]
        _, τ = _partition_shortwave_series(sw, toa, sza, fill(P, 6), fill(e, 6))
        @test τ ≈ [1.0, 1.0, 1.0, 1.0, 5.0, 5.0] rtol = 1e-12

        # Night shortwave, normally zero, is all diffuse with τ = 0; NaN shortwave is NaN in both.
        diffuse, τ = _partition_shortwave_series([0.4, NaN, NaN], [0.0, 500.0, 50.0],
                                                 [90.0, 60.0, 85.0], fill(P, 3), fill(e, 3))
        @test diffuse[1] == 0.4 && τ[1] == 0.0
        @test all(isnan, diffuse[2:3]) && all(isnan, τ[2:3])
        @test_throws DimensionMismatch _partition_shortwave_series(zeros(3), zeros(3), zeros(3),
                                                                   zeros(2), zeros(3))
    end

    @testset "shortwave_forcing_layers" begin
        times = DateTime(2020, 7, 1) .+ Hour.(1:48)
        g = solar_geometry(46.5, 8.0, times)
        sw = 0.6 .* g.shortwave_toa
        P, e = fill(75000.0, 48), fill(500.0, 48)
        L = shortwave_forcing_layers(sw, P, e, 46.5, 8.0, times)
        @test L.solar_zenith_angle == g.solar_zenith_angle
        @test all(0 .<= L.shortwave_downward_diffuse .<= sw)
        @test all(0 .<= L.cloud_optical_thickness .<= CLOUD_OPTICAL_THICKNESS_MAX)
        @test (L.shortwave_downward_diffuse, L.cloud_optical_thickness) ==
              _partition_shortwave_series(sw, g.shortwave_toa, g.solar_zenith_angle, P, e)
        # An inverted hour's τ is 0.6 against the FAO clear sky at its own zenith angle.
        i = argmin(g.solar_zenith_angle)
        k_clear = fao_clearness(cosd(g.solar_zenith_angle[i]), 75.0, 0.5)
        @test L.cloud_optical_thickness[i] ≈ (k_clear / 0.6 - 1) / 0.13 rtol = 1e-12
        # The same series as views into longer vectors gives the same layers.
        pad(v) = view(vcat([NaN], v, [NaN]), 2:49)
        Lv = shortwave_forcing_layers(pad(sw), pad(P), pad(e), 46.5, 8.0, times)
        @test Lv.cloud_optical_thickness == L.cloud_optical_thickness
        @test Lv.shortwave_downward_diffuse == L.shortwave_downward_diffuse
        @test_throws DimensionMismatch shortwave_forcing_layers(sw, P[1:47], e, 46.5, 8.0, times)
    end

    @testset "the layers ride through every adjustment unchanged" begin
        times = DateTime(2020, 7, 1) .+ Hour.(1:48)
        n = length(times)
        L = shortwave_forcing_layers(0.6 .* solar_geometry(46.5, 8.0, times).shortwave_toa,
                                     fill(75000.0, n), fill(500.0, n), 46.5, 8.0, times)
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
        @test_throws "cloud_optical_thickness: expected [0, 30.0]" validate_climate_forcing_units(
            bad(:cloud_optical_thickness, fill(30.5, n)))
        @test validate_climate_forcing_units(bad(:cloud_optical_thickness, fill(30.0, n)))
    end

    @testset "simulated forcing carries the layers" begin
        ds = simulate_climate_forcing("test_1", 3)
        for v in (:solar_zenith_angle, :shortwave_downward_diffuse, :cloud_optical_thickness)
            @test haskey(ds, v)
        end
        sw = parent(ds[:shortwave_downward])
        @test all(parent(ds[:solar_zenith_angle])[sw .== 0] .== 90.0)
        @test all(0 .<= parent(ds[:shortwave_downward_diffuse]) .<= sw)
        τ = parent(ds[:cloud_optical_thickness])
        @test all(0 .<= τ .<= CLOUD_OPTICAL_THICKNESS_MAX)
        @test all(iszero, τ[sw .== 0])
    end
end
