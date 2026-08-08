#!/usr/bin/env julia

"""
Analytic (offline) tests for elevation adjustment of climate forcing.
No network access required.
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData

using GEMB_ClimateForcing: saturation_vapor_pressure,
    konzelmann_clear_sky_emissivity, empirical_lapse_rate,
    GREENLAND_LAPSE_RATE, ARCTIC_LAPSE_RATE, ANTARCTICA_LAPSE_RATE,
    _resolve_lapse_rate, _G0, _R_DRY, _SIGMA_SB

# Build a physically reasonable stack spanning n_months so month-dependent
# monthly lapse tables are exercised.
function make_forcing_stack(; n_months::Int=12)
    times = [DateTime(2020, m, 15, 12) for m in 1:n_months]
    time_dim = Ti(times)
    npts = length(times)
    T  = fill(260.0, npts)              # K (cold, so over-ice curve is used)
    P  = fill(80000.0, npts)            # Pa (high-elevation surface pressure)
    e  = fill(150.0, npts)              # Pa (low vapor pressure)
    ws = fill(6.0, npts)                # m/s
    pr = fill(0.2, npts)                # kg/m²/hr
    sw = fill(150.0, npts)              # W/m²
    lw = fill(220.0, npts)              # W/m²
    return DimStack((
        temperature_air = DimArray(T, (time_dim,); metadata=Dict("units" => "K")),
        pressure_air = DimArray(P, (time_dim,); metadata=Dict("units" => "Pa")),
        vapor_pressure = DimArray(e, (time_dim,); metadata=Dict("units" => "Pa")),
        wind_speed = DimArray(ws, (time_dim,); metadata=Dict("units" => "m/s")),
        precipitation = DimArray(pr, (time_dim,); metadata=Dict("units" => "kg/m²/hr")),
        shortwave_downward = DimArray(sw, (time_dim,); metadata=Dict("units" => "W/m²")),
        longwave_downward = DimArray(lw, (time_dim,); metadata=Dict("units" => "W/m²")),
    ); metadata=Dict("elevation" => 1500.0, "dataset" => "test"))
end

@testset "Elevation Adjustment" begin

    @testset "Saturation vapor pressure" begin
        # Over water at 0 °C ≈ 611.21 Pa; over ice ≈ 611.15 Pa (Buck 1981).
        @test saturation_vapor_pressure(273.15; over_ice=false) ≈ 611.21 rtol=1e-3
        @test saturation_vapor_pressure(273.15; over_ice=true) ≈ 611.15 rtol=1e-3
        # Monotonic increase with temperature.
        @test saturation_vapor_pressure(283.15) > saturation_vapor_pressure(273.15)
        # Below freezing, ice curve is below water curve at the same T.
        @test saturation_vapor_pressure(253.15; over_ice=true) <
              saturation_vapor_pressure(253.15; over_ice=false)
        # Auto-selection: default uses ice below 0 °C.
        @test saturation_vapor_pressure(253.15) ==
              saturation_vapor_pressure(253.15; over_ice=true)
    end

    @testset "Konzelmann clear-sky emissivity" begin
        ε = konzelmann_clear_sky_emissivity(200.0, 260.0)
        @test 0.0 < ε < 1.0
        @test konzelmann_clear_sky_emissivity(400.0, 260.0) >
              konzelmann_clear_sky_emissivity(100.0, 260.0)
    end

    @testset "Named lapse-rate constants (K/km, positive)" begin
        for Γ in (GREENLAND_LAPSE_RATE, ARCTIC_LAPSE_RATE, ANTARCTICA_LAPSE_RATE)
            @test length(Γ) == 12
            @test all(Γ .> 0)                 # positive = cooling with height
        end
        # Antarctic interior is the steepest; Greenland summer is shallow.
        @test minimum(ANTARCTICA_LAPSE_RATE) > maximum(GREENLAND_LAPSE_RATE)
        @test GREENLAND_LAPSE_RATE[7] < GREENLAND_LAPSE_RATE[2]  # Jul < Feb
    end

    @testset "Empirical (locally-fitted) lapse rate" begin
        # T = -0.006·z + 290  ⇒  lapse rate = +6 K/km.
        z = [0.0, 100.0, 200.0, 300.0, 400.0, 500.0]
        T = @. -0.006 * z + 290.0
        @test empirical_lapse_rate(T, z) ≈ 6.0
        # Robust to noise.
        Tn = T .+ [0.1, -0.05, 0.08, -0.1, 0.03, -0.02]
        @test empirical_lapse_rate(Tn, z) ≈ 6.0 atol=0.5
        # Degenerate inputs are rejected.
        @test_throws ArgumentError empirical_lapse_rate([1.0], [1.0])
        @test_throws ArgumentError empirical_lapse_rate([1.0, 2.0], [5.0, 5.0])
        @test_throws ArgumentError empirical_lapse_rate([1.0, 2.0], [1.0])
        # Composes with the adjuster (positive rate cools going up).
        stack = make_forcing_stack()
        Γ = empirical_lapse_rate(T, z)
        adj = climate_adjust_for_elevation(stack, 1000.0; lapse_rate=Γ)
        ΔT = parent(adj[:temperature_air]) .- parent(stack[:temperature_air])
        @test all(ΔT .≈ -(Γ / 1000.0) * 1000.0)   # = -6 K
    end

    @testset "Lapse-rate argument forms" begin
        months = [1, 4, 7, 10]
        # Scalar broadcasts.
        @test _resolve_lapse_rate(6.5, months, 4) == fill(6.5, 4)
        # Length-12 treated as monthly, indexed by month.
        monthly = collect(1.0:12.0)
        @test _resolve_lapse_rate(monthly, months, 4) == [1.0, 4.0, 7.0, 10.0]
        # Length matching the record is per-time-step.
        perstep = [5.0, 6.0, 7.0, 8.0]
        @test _resolve_lapse_rate(perstep, months, 4) == perstep
        # Mismatched length errors.
        @test_throws ArgumentError _resolve_lapse_rate([1.0, 2.0, 3.0], months, 4)
    end

    @testset "Lapse-rate physical bounds check" begin
        stack = make_forcing_stack()
        # Documented extremes are accepted: Antarctic interior (14.3) and strong
        # inversions (negative).
        @test climate_adjust_for_elevation(stack, 100.0; lapse_rate=14.3) isa DimStack
        @test climate_adjust_for_elevation(stack, 100.0; lapse_rate=-10.0) isa DimStack
        @test climate_adjust_for_elevation(stack, 100.0; lapse_rate=ANTARCTICA_LAPSE_RATE) isa DimStack
        # Out-of-range scalar (unit slip: K/m ~ 0.0065 is fine numerically, but a
        # wildly large value or wrong sign is rejected).
        @test_throws ArgumentError climate_adjust_for_elevation(stack, 100.0; lapse_rate=50.0)
        @test_throws ArgumentError climate_adjust_for_elevation(stack, 100.0; lapse_rate=-40.0)
        # Out-of-range value inside an otherwise-valid monthly vector is caught.
        bad_monthly = fill(6.5, 12); bad_monthly[6] = 99.0
        @test_throws ArgumentError climate_adjust_for_elevation(stack, 100.0; lapse_rate=bad_monthly)
    end

    @testset "Zero delta is identity" begin
        stack = make_forcing_stack()
        adj = climate_adjust_for_elevation(stack, 0.0)
        @test parent(adj[:temperature_air]) ≈ parent(stack[:temperature_air])
        @test parent(adj[:pressure_air]) ≈ parent(stack[:pressure_air])
        @test parent(adj[:vapor_pressure]) ≈ parent(stack[:vapor_pressure])
        @test parent(adj[:longwave_downward]) ≈ parent(stack[:longwave_downward])
        @test parent(adj[:precipitation]) ≈ parent(stack[:precipitation])
        @test parent(adj[:wind_speed]) == parent(stack[:wind_speed])
        @test parent(adj[:shortwave_downward]) == parent(stack[:shortwave_downward])
    end

    @testset "Temperature lapse sign and magnitude" begin
        stack = make_forcing_stack()
        Δz = 500.0
        adj = climate_adjust_for_elevation(stack, Δz; lapse_rate=6.5)
        # 6.5 K/km × 0.5 km = 3.25 K cooling everywhere.
        @test all(parent(adj[:temperature_air]) .≈ (parent(stack[:temperature_air]) .- 3.25))
        # Going up cools; going down warms.
        adj_down = climate_adjust_for_elevation(stack, -Δz; lapse_rate=6.5)
        @test all(parent(adj_down[:temperature_air]) .> parent(stack[:temperature_air]))
    end

    @testset "Monthly lapse varies by month" begin
        stack = make_forcing_stack()
        adj = climate_adjust_for_elevation(stack, 1000.0; lapse_rate=GREENLAND_LAPSE_RATE)
        ΔT = parent(adj[:temperature_air]) .- parent(stack[:temperature_air])
        # Greenland July (index 7) cooling is weaker than February (index 2).
        @test abs(ΔT[7]) < abs(ΔT[2])
        # Matches the table exactly: ΔT = -(Γ/1000)·1000 = -Γ.
        @test ΔT ≈ -GREENLAND_LAPSE_RATE
    end

    @testset "Per-time-step lapse rate" begin
        # Use a record whose length != 12 (6 steps) so a length-6 vector is
        # unambiguously per-time-step, not monthly.
        stack6 = make_forcing_stack(n_months=6)
        perstep = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]   # K/km, length 6 == record
        adj = climate_adjust_for_elevation(stack6, 1000.0; lapse_rate=perstep)
        ΔT = parent(adj[:temperature_air]) .- parent(stack6[:temperature_air])
        @test ΔT ≈ -perstep
    end

    @testset "Pressure barometric law" begin
        stack = make_forcing_stack()
        Δz = 300.0
        adj = climate_adjust_for_elevation(stack, Δz; lapse_rate=6.5)
        P  = parent(stack[:pressure_air])
        T  = parent(stack[:temperature_air])
        T′ = parent(adj[:temperature_air])
        expected = @. P * exp(-_G0 * Δz / (_R_DRY * 0.5 * (T + T′)))
        @test parent(adj[:pressure_air]) ≈ expected
        @test all(parent(adj[:pressure_air]) .< P)   # higher target ⇒ lower P
    end

    @testset "Constant relative humidity preserved" begin
        stack = make_forcing_stack()
        adj = climate_adjust_for_elevation(stack, 400.0; lapse_rate=GREENLAND_LAPSE_RATE)
        T  = parent(stack[:temperature_air]); e  = parent(stack[:vapor_pressure])
        T′ = parent(adj[:temperature_air]);   e′ = parent(adj[:vapor_pressure])
        RH  = e ./ saturation_vapor_pressure.(T)
        RH′ = e′ ./ saturation_vapor_pressure.(T′)
        @test RH′ ≈ RH
        @test all(e′ .< e)                     # cooling aloft ⇒ lower vapor pressure
    end

    @testset "Longwave recompute preserves cloud increment" begin
        stack = make_forcing_stack()
        adj0 = climate_adjust_for_elevation(stack, 0.0)
        @test parent(adj0[:longwave_downward]) ≈ parent(stack[:longwave_downward])
        adj = climate_adjust_for_elevation(stack, 500.0; lapse_rate=GREENLAND_LAPSE_RATE)
        @test all(parent(adj[:longwave_downward]) .< parent(stack[:longwave_downward]))
    end

    @testset "Precipitation scaling method" begin
        stack = make_forcing_stack()
        # Default: precipitation unchanged.
        adj_none = climate_adjust_for_elevation(stack, 500.0)
        @test parent(adj_none[:precipitation]) ≈ parent(stack[:precipitation])
        # Clausius–Clapeyron: cooling aloft ⇒ less precipitation.
        adj_cc = climate_adjust_for_elevation(stack, 500.0;
                                              lapse_rate=6.5,
                                              precip_scaling_method=:clausius_clapeyron)
        @test all(parent(adj_cc[:precipitation]) .< parent(stack[:precipitation]))
        # Going down (warmer) ⇒ more precipitation.
        adj_dn = climate_adjust_for_elevation(stack, -500.0;
                                              lapse_rate=6.5,
                                              precip_scaling_method=:clausius_clapeyron)
        @test all(parent(adj_dn[:precipitation]) .> parent(stack[:precipitation]))
        # Δz = 0 identity even with scaling on.
        adj0 = climate_adjust_for_elevation(stack, 0.0;
                                            precip_scaling_method=:clausius_clapeyron)
        @test parent(adj0[:precipitation]) ≈ parent(stack[:precipitation])
        # Unknown method rejected.
        @test_throws ArgumentError climate_adjust_for_elevation(
            stack, 100.0; precip_scaling_method=:bogus)
    end

    @testset "Metadata bookkeeping" begin
        stack = make_forcing_stack()
        adj = climate_adjust_for_elevation(stack, 1000.0; lapse_rate=5.0)
        ΔT = parent(adj[:temperature_air]) .- parent(stack[:temperature_air])
        @test all(ΔT .≈ -5.0)
        meta = metadata(adj)
        @test meta["delta_elevation"] == 1000.0
        @test meta["elevation"] == 2500.0            # 1500 + 1000
        @test meta["elevation_reanalysis"] == 1500.0
    end

end

println("\n✓ All elevation adjustment tests passed!")
