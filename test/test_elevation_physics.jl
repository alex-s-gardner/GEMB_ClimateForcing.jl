#!/usr/bin/env julia

"""
Analytical physics tests for elevation adjustment.

Each testset validates one physics method against an *independently* coded
closed-form reference (published equations/constants written out fresh here, not
reusing the module's implementation), so a coefficient or sign error in the
source is caught rather than mirrored. No network access required.

References cross-checked:
- Glover (1999), J. Climate 12, 551–563 (Eqs. 15–20): adjustment scheme.
- Buck (1981), J. Appl. Meteorol. 20, 1527–1532: saturation vapor pressure.
- Konzelmann et al. (1994), Glob. Planet. Change 9, 143–164: clear-sky emissivity.
- Standard atmosphere / hydrostatic equation: barometric pressure.
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData

using GEMB_ClimateForcing: saturation_vapor_pressure,
    konzelmann_clear_sky_emissivity, empirical_lapse_rate,
    GREENLAND_LAPSE_RATE, ARCTIC_LAPSE_RATE, ANTARCTICA_LAPSE_RATE

# ---------------------------------------------------------------------------
# Independent reference implementations (published constants, coded fresh).
# ---------------------------------------------------------------------------
const G_REF = 9.80665       # m/s²  (standard gravity)
const RD_REF = 287.05       # J/(kg·K)
const SIGMA_REF = 5.67e-8   # W/m²/K⁴

# Buck (1981) saturation vapor pressure, Pa (Tc in °C).
ref_esw(Tc) = 611.21 * exp(17.502 * Tc / (240.97 + Tc))   # over water
ref_esi(Tc) = 611.15 * exp(22.452 * Tc / (272.55 + Tc))   # over ice

# Konzelmann et al. (1994) clear-sky emissivity (e in Pa, T in K).
ref_konz(e, T) = 0.23 + 0.484 * (e / T)^(1 / 8)

# Build a single-step stack at prescribed state (one timestep in given month).
function stack_at(; T=260.0, P=80000.0, e=150.0, ws=6.0, pr=0.2, sw=150.0, lw=220.0,
                    month::Int=1, elevation=1500.0)
    time_dim = Ti([DateTime(2020, month, 15, 12)])
    return DimStack((
        temperature_air = DimArray([T], (time_dim,); metadata=Dict("units" => "K")),
        pressure_air = DimArray([P], (time_dim,); metadata=Dict("units" => "Pa")),
        vapor_pressure = DimArray([e], (time_dim,); metadata=Dict("units" => "Pa")),
        wind_speed = DimArray([ws], (time_dim,); metadata=Dict("units" => "m/s")),
        precipitation = DimArray([pr], (time_dim,); metadata=Dict("units" => "kg/m²/hr")),
        shortwave_downward = DimArray([sw], (time_dim,); metadata=Dict("units" => "W/m²")),
        longwave_downward = DimArray([lw], (time_dim,); metadata=Dict("units" => "W/m²")),
    ); metadata=Dict("elevation" => elevation))
end

# Convenience: pull the single adjusted scalar for a variable.
val(stack, var) = parent(stack[var])[1]

@testset "Elevation Physics (analytical)" begin

    # -----------------------------------------------------------------------
    # 1. TEMPERATURE — T′ = T − (Γ/1000)·Δz  (Glover 1999 Eq. 15).
    # -----------------------------------------------------------------------
    @testset "Temperature lapse (Glover Eq. 15)" begin
        T0 = 260.0
        for (Γ, Δz) in [(6.5, 500.0), (4.6, 1000.0), (14.3, -250.0), (0.0, 800.0)]
            adj = climate_adjust_for_elevation(stack_at(T=T0), Δz; lapse_rate=Γ)
            @test val(adj, :temperature_air) ≈ T0 - (Γ / 1000.0) * Δz
        end
        # Sign: ascending (Δz>0) cools, descending warms.
        @test val(climate_adjust_for_elevation(stack_at(T=T0), 300.0; lapse_rate=6.5), :temperature_air) < T0
        @test val(climate_adjust_for_elevation(stack_at(T=T0), -300.0; lapse_rate=6.5), :temperature_air) > T0
    end

    # -----------------------------------------------------------------------
    # 2. PRESSURE — hydrostatic P′ = P·exp(−g·Δz/(R_d·T̄)), T̄ = ½(T+T′).
    #    With Γ=0, T̄ = T exactly ⇒ closed-form isothermal barometric law.
    # -----------------------------------------------------------------------
    @testset "Pressure barometric (isothermal, Γ=0)" begin
        T0 = 273.15; P0 = 90000.0; Δz = 1000.0
        adj = climate_adjust_for_elevation(stack_at(T=T0, P=P0), Δz; lapse_rate=0.0)
        H = RD_REF * T0 / G_REF                       # scale height
        @test val(adj, :pressure_air) ≈ P0 * exp(-Δz / H)
        # Independent literal: at T=273.15, dz=1000 ⇒ ratio ≈ 0.88243.
        @test val(adj, :pressure_air) / P0 ≈ 0.88243 atol = 1e-4
    end

    @testset "Pressure barometric (with lapse, mean-layer T)" begin
        T0 = 260.0; P0 = 80000.0; Δz = 400.0; Γ = 6.5
        adj = climate_adjust_for_elevation(stack_at(T=T0, P=P0), Δz; lapse_rate=Γ)
        T′ = T0 - (Γ / 1000.0) * Δz
        Tbar = 0.5 * (T0 + T′)
        @test val(adj, :pressure_air) ≈ P0 * exp(-G_REF * Δz / (RD_REF * Tbar))
        # First-order equivalence to Glover Eq. 16 linear form p − gρΔz for small Δz.
        ρ = P0 / (RD_REF * Tbar)
        @test val(adj, :pressure_air) ≈ (P0 - G_REF * ρ * Δz) rtol = 5e-3
    end

    # -----------------------------------------------------------------------
    # 3. SATURATION VAPOR PRESSURE — Buck (1981) exact reference points.
    # -----------------------------------------------------------------------
    @testset "Saturation vapor pressure (Buck 1981)" begin
        # Literal reference values (Pa) computed from Buck coefficients.
        @test saturation_vapor_pressure(273.15; over_ice=false) ≈ 611.21 atol = 1e-2
        @test saturation_vapor_pressure(293.15; over_ice=false) ≈ 2337.28 atol = 1e-2
        @test saturation_vapor_pressure(303.15; over_ice=false) ≈ 4243.51 atol = 1e-1
        @test saturation_vapor_pressure(263.15; over_ice=true) ≈ 259.872 atol = 1e-2
        @test saturation_vapor_pressure(253.15; over_ice=true) ≈ 103.267 atol = 1e-2
        # Matches the independent reference formulas across a range.
        for Tc in -40.0:5.0:40.0
            @test saturation_vapor_pressure(273.15 + Tc; over_ice=false) ≈ ref_esw(Tc)
        end
        for Tc in -50.0:5.0:-1.0
            @test saturation_vapor_pressure(273.15 + Tc; over_ice=true) ≈ ref_esi(Tc)
        end
        # Clausius–Clapeyron sanity: water es rises ~7.5%/K near 0 °C.
        @test saturation_vapor_pressure(274.15; over_ice=false) /
              saturation_vapor_pressure(273.15; over_ice=false) ≈ 1.075 atol = 2e-3
        # Ice curve below water curve below freezing.
        @test saturation_vapor_pressure(253.15; over_ice=true) <
              saturation_vapor_pressure(253.15; over_ice=false)
    end

    # -----------------------------------------------------------------------
    # 4. CONSTANT RELATIVE HUMIDITY — RH preserved; e′/e = es(T′)/es(T).
    # -----------------------------------------------------------------------
    @testset "Constant-RH humidity (Glover Eq. 20)" begin
        T0 = 260.0; e0 = 150.0; Δz = 400.0; Γ = 6.5
        adj = climate_adjust_for_elevation(stack_at(T=T0, e=e0), Δz; lapse_rate=Γ)
        T′ = T0 - (Γ / 1000.0) * Δz
        # RH preserved exactly (both use over-ice curve here, T<0 °C).
        RH0 = e0 / ref_esi(T0 - 273.15)
        e′ = val(adj, :vapor_pressure)
        RH′ = e′ / ref_esi(T′ - 273.15)
        @test RH′ ≈ RH0
        # e′ equals RH·es(T′) with the published ice coefficients.
        @test e′ ≈ RH0 * ref_esi(T′ - 273.15)
        # Ratio form matches Glover Eq. 20 (q ∝ es proxy).
        @test e′ / e0 ≈ ref_esi(T′ - 273.15) / ref_esi(T0 - 273.15)
        # Supersaturation guard: a supersaturated input (e > es(T), RH>1) is
        # clamped to RH=1 so the adjusted vapor pressure never exceeds es(T′).
        e_super = 1.2 * ref_esi(260.0 - 273.15)         # RH ≈ 1.2 at the source
        wet = climate_adjust_for_elevation(stack_at(T=260.0, e=e_super), 200.0; lapse_rate=6.5)
        T′w = 260.0 - 6.5 / 1000.0 * 200.0
        @test val(wet, :vapor_pressure) ≈ saturation_vapor_pressure(T′w)  # clamped to RH=1
    end

    # -----------------------------------------------------------------------
    # 5. LONGWAVE — Konzelmann clear-sky emissivity + preserved cloud increment.
    # -----------------------------------------------------------------------
    @testset "Konzelmann clear-sky emissivity (1994)" begin
        @test konzelmann_clear_sky_emissivity(200.0, 260.0) ≈ 0.69838 atol = 1e-5
        @test konzelmann_clear_sky_emissivity(200.0, 260.0) ≈ ref_konz(200.0, 260.0)
        # Emissivity increases with vapor pressure (more emitters), 0<ε<1.
        @test 0 < konzelmann_clear_sky_emissivity(50.0, 260.0) < 1
        @test konzelmann_clear_sky_emissivity(500.0, 260.0) >
              konzelmann_clear_sky_emissivity(50.0, 260.0)
    end

    @testset "Longwave recompute (clear-sky exact, Δε=0)" begin
        # Construct a stack whose LW is EXACTLY clear-sky (Δε = 0), so the
        # recompute reduces to ε_cs(e′,T′)·σ·T′⁴ with no cloud residual.
        T0 = 260.0; e0 = 150.0
        LW_clear = ref_konz(e0, T0) * SIGMA_REF * T0^4
        Δz = 500.0; Γ = 6.5
        adj = climate_adjust_for_elevation(stack_at(T=T0, e=e0, lw=LW_clear), Δz; lapse_rate=Γ)
        T′ = T0 - (Γ / 1000.0) * Δz
        RH0 = e0 / ref_esi(T0 - 273.15)
        e′ = RH0 * ref_esi(T′ - 273.15)
        LW_expected = ref_konz(e′, T′) * SIGMA_REF * T′^4
        @test val(adj, :longwave_downward) ≈ LW_expected
        # Cooling aloft reduces clear-sky LW.
        @test val(adj, :longwave_downward) < LW_clear
    end

    @testset "Longwave preserves cloud increment Δε" begin
        # All-sky LW above clear-sky (cloudy). The cloud emissivity increment
        # Δε = LW/(σT⁴) − ε_cs must be carried unchanged to the new level.
        T0 = 260.0; e0 = 150.0
        LW_all = 300.0                                  # > clear-sky ⇒ Δε > 0
        Δε = LW_all / (SIGMA_REF * T0^4) - ref_konz(e0, T0)
        @test Δε > 0
        Δz = 500.0; Γ = 6.5
        adj = climate_adjust_for_elevation(stack_at(T=T0, e=e0, lw=LW_all), Δz; lapse_rate=Γ)
        T′ = T0 - (Γ / 1000.0) * Δz
        RH0 = e0 / ref_esi(T0 - 273.15)
        e′ = RH0 * ref_esi(T′ - 273.15)
        LW_expected = (ref_konz(e′, T′) + Δε) * SIGMA_REF * T′^4
        @test val(adj, :longwave_downward) ≈ LW_expected
        # Recovered increment matches at the new level.
        Δε′ = val(adj, :longwave_downward) / (SIGMA_REF * T′^4) - ref_konz(e′, T′)
        @test Δε′ ≈ Δε
    end

    @testset "Longwave identity at Δz=0" begin
        adj = climate_adjust_for_elevation(stack_at(lw=245.0), 0.0)
        @test val(adj, :longwave_downward) ≈ 245.0
    end

    # -----------------------------------------------------------------------
    # 6. PRECIPITATION — Clausius–Clapeyron scaling (Glover Eq. 19).
    # -----------------------------------------------------------------------
    @testset "Clausius–Clapeyron precipitation (Glover Eq. 19)" begin
        T0 = 260.0; pr0 = 0.2; Δz = 500.0; Γ = 6.5
        adj = climate_adjust_for_elevation(stack_at(T=T0, pr=pr0), Δz;
                                           lapse_rate=Γ, precip_scaling_method=:clausius_clapeyron)
        T′ = T0 - (Γ / 1000.0) * Δz
        # P′/P = es(T′)/es(T) with published ice coefficients.
        @test val(adj, :precipitation) ≈ pr0 * ref_esi(T′ - 273.15) / ref_esi(T0 - 273.15)
        @test val(adj, :precipitation) < pr0            # cooler aloft ⇒ less precip
        # Default (nothing) leaves precip untouched.
        adj0 = climate_adjust_for_elevation(stack_at(T=T0, pr=pr0), Δz; lapse_rate=Γ)
        @test val(adj0, :precipitation) ≈ pr0
        # Δz = 0 identity even with scaling on.
        adjz = climate_adjust_for_elevation(stack_at(T=T0, pr=pr0), 0.0;
                                            precip_scaling_method=:clausius_clapeyron)
        @test val(adjz, :precipitation) ≈ pr0
    end

    # -----------------------------------------------------------------------
    # 7. EMPIRICAL LAPSE RATE — OLS slope, returned K/km positive = cooling.
    # -----------------------------------------------------------------------
    @testset "Empirical lapse rate (OLS, RACMO-style)" begin
        # Exact linear field T = −0.0065·z + 285  ⇒  lapse = +6.5 K/km.
        z = collect(0.0:100.0:1000.0)
        T = @. -0.0065 * z + 285.0
        @test empirical_lapse_rate(T, z) ≈ 6.5
        # Inversion (T increases with height) ⇒ negative lapse rate.
        Tinv = @. 0.003 * z + 250.0
        @test empirical_lapse_rate(Tinv, z) ≈ -3.0
        # Matches the analytic OLS slope formula independently.
        zr = [10.0, 40.0, 55.0, 90.0, 120.0]
        Tr = [270.0, 268.5, 268.0, 266.9, 265.5]
        z̄ = sum(zr) / length(zr); T̄ = sum(Tr) / length(Tr)
        slope = sum((zr .- z̄) .* (Tr .- T̄)) / sum((zr .- z̄) .^ 2)
        @test empirical_lapse_rate(Tr, zr) ≈ -1000.0 * slope
    end

    # -----------------------------------------------------------------------
    # 8. UNCHANGED VARIABLES — shortwave and wind pass through untouched.
    # -----------------------------------------------------------------------
    @testset "Shortwave and wind unchanged" begin
        adj = climate_adjust_for_elevation(stack_at(sw=175.0, ws=8.0), 750.0;
                                           lapse_rate=6.5, precip_scaling_method=:clausius_clapeyron)
        @test val(adj, :shortwave_downward) == 175.0
        @test val(adj, :wind_speed) == 8.0
    end

    # -----------------------------------------------------------------------
    # 9. NAMED CONSTANTS — magnitudes/ordering consistent with the literature.
    # -----------------------------------------------------------------------
    @testset "Literature lapse-rate constants" begin
        # Greenland (Fausto 2009 Table 4): summer shallow, winter steep;
        # annual mean 6.8 K/km. Exact table values (negated to positive=cooling).
        @test GREENLAND_LAPSE_RATE == [7.9, 8.9, 7.9, 7.3, 5.9, 4.7, 4.6, 5.7, 6.9, 7.3, 6.5, 7.6]
        @test GREENLAND_LAPSE_RATE[7] ≈ 4.6            # July, shallowest
        @test GREENLAND_LAPSE_RATE[2] ≈ 8.9            # February, steepest
        @test sum(GREENLAND_LAPSE_RATE) / 12 ≈ 6.8 atol = 0.05
        # Arctic (Gardner 2009): summer mean ~4.9, winter ~3.2 K/km.
        @test ARCTIC_LAPSE_RATE[7] ≈ 4.9
        @test ARCTIC_LAPSE_RATE[1] ≈ 3.2
        # Antarctica interior (Fortuin & Oerlemans 1990): steep ~14.3 K/km.
        @test all(ANTARCTICA_LAPSE_RATE .≈ 14.3)
        # Ordering: Antarctic interior > Greenland > Arctic (summer melt regimes).
        @test minimum(ANTARCTICA_LAPSE_RATE) > maximum(GREENLAND_LAPSE_RATE)
        @test sum(GREENLAND_LAPSE_RATE) > sum(ARCTIC_LAPSE_RATE)
    end

end

println("\n✓ All analytical physics tests passed!")
