#!/usr/bin/env julia

"""
Analytic (offline) tests for the per-region downscaling-parameter fits
(`derive_decoupling_factor`, `derive_lapse_rate`, `decoupling_factor_at_elevation`).

The fits take an accumulator of per-timestep cross-cell regression sums rather than forcing
stacks, so these drive them directly from a hand-built `acc` — no table, no loader, no network.
`_acc` below is the consumer-side statement of that contract: if a producer's sums stop matching
it, the recovery tests here fail rather than a downstream aggregation quietly returning `NaN`.
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData
using Statistics

using GEMB_ClimateForcing: _MIN_AMBIENT_EXCESS, _MIN_CELLS_DEFAULT, _LAPSE_RATE_LIMITS,
    _DECOUPLING_FACTOR_LIMITS, _DECOUPLING_REFERENCE_TEMPERATURE, _solve_sym4, _r2,
    _make_applicable, _decouple_per_timestep, saturation_vapor_pressure,
    _downscaling_sums, _accumulate_downscaling_cell!

const _DP_TIME = collect(DateTime(2000, 1, 1):Hour(1):DateTime(2000, 1, 1, 5))

# The cross-cell sums the two fits read, accumulated exactly as a streaming producer does: the
# design vector is `x = [1, z, glm, glm*z]`, the responses are the warm excess `E` and the
# temperature `T`, and the decoupling fit's sums run over the cells above melting only. `cells` is
# a vector of `(; z, glm)`; `temperature(cell)` returns that cell's temperature vector over `time`.
function _acc(cells, temperature; time = _DP_TIME)
    n = length(time)
    S = (; s1 = zeros(Int, n), z = zeros(n), zz = zeros(n), T = zeros(n), zT = zeros(n),
         g = zeros(n), w = zeros(n), zw = zeros(n),
         s1_warm = zeros(Int, n), z_warm = zeros(n), g_warm = zeros(n), w_warm = zeros(n),
         zz_warm = zeros(n), zg_warm = zeros(n), zw_warm = zeros(n),
         gg_warm = zeros(n), gw_warm = zeros(n), ww_warm = zeros(n),
         E = zeros(n), zE = zeros(n), gE = zeros(n), wE = zeros(n), EE = zeros(n))
    for c in cells
        z, glm, gz = Float64(c.z), Float64(c.glm), Float64(c.glm) * Float64(c.z)
        T = temperature(c)
        for t in 1:n
            Tt = T[t]
            S.s1[t] += 1
            S.z[t] += z; S.zz[t] += z * z; S.T[t] += Tt; S.zT[t] += z * Tt
            S.g[t] += glm; S.w[t] += gz; S.zw[t] += z * gz
            Et = Tt - _DECOUPLING_REFERENCE_TEMPERATURE
            Et > 0 || continue
            S.s1_warm[t] += 1
            S.z_warm[t]  += z;         S.g_warm[t]  += glm;        S.w_warm[t]  += gz
            S.zz_warm[t] += z * z;     S.zg_warm[t] += z * glm;    S.zw_warm[t] += z * gz
            S.gg_warm[t] += glm * glm; S.gw_warm[t] += glm * gz;   S.ww_warm[t] += gz * gz
            S.E[t]  += Et;        S.zE[t] += z * Et
            S.gE[t] += glm * Et;  S.wE[t] += gz * Et
            S.EE[t] += Et * Et
        end
    end
    # Unweighted mean elevation: these fixtures give every cell the same glacier area, so the
    # area-weighted reference a real producer computes reduces to this.
    return (; time, n, sums = S,
            mean_elevation = Statistics.mean(Float64(c.z) for c in cells))
end

# Four cells spanning 1000-2500 m with a range of glacier fractions — the same fixture geometry the
# downstream regional derivation is tested on.
const _DP_CELLS = [(z = 1000.0, glm = 0.0),
                   (z = 1500.0, glm = 0.5),
                   (z = 2000.0, glm = 1.0),
                   (z = 2500.0, glm = 0.25)]

const _DP_K_TRUE = 0.7
const _DP_GAMMA_TRUE = 5.0      # K/km, ambient
# Ambient temperature at sea level, high enough that every cell stays above melting: the estimator
# assumes the ambient excess is linear in elevation, and it is only linear while no cell clips at
# the melting point.
const _DP_T0 = 292.0

_dp_ambient(c) = _DP_T0 - _DP_GAMMA_TRUE * c.z / 1000.0

# What the reanalysis reports: a glm=1 cell it runs as ice is fully decoupled, a glm=0 cell it runs
# as bare ground is left ambient, and the correction it has already applied scales with glm.
function _dp_observed(c)
    a = _dp_ambient(c)
    k_applied = 1 - c.glm * (1 - _DP_K_TRUE)
    return fill(a + (k_applied - 1) * max(a - _DECOUPLING_REFERENCE_TEMPERATURE, 0.0),
                length(_DP_TIME))
end

@testset "Downscaling Parameters" begin
    acc = _acc(_DP_CELLS, _dp_observed)
    decoupling = derive_decoupling_factor(acc; min_cells = 4)
    lapse = derive_lapse_rate(acc, decoupling; min_cells = 4)

    @testset "decoupling factor" begin
        # The fixture is exactly the model the estimator inverts, so `k` comes back to roundoff.
        @test all(k -> isapprox(k, _DP_K_TRUE, atol = 1e-6), decoupling.decoupling_factor)
        @test DimensionalData.metadata(decoupling)["n_fitted"] == length(_DP_TIME)
        @test DimensionalData.metadata(decoupling)["reference_elevation"] == acc.mean_elevation
        @test DimensionalData.metadata(decoupling)["min_cells"] == 4
        @test all(==(4), decoupling.n_cells)
        @test collect(dims(decoupling, Ti)) == _DP_TIME

        # A four-term surface fit through four points that lie exactly on it: no residual.
        @test all(r -> isapprox(r, 1.0, atol = 1e-8), decoupling.r2)

        # `ambient_excess` is the fitted ice-free excess at the reference elevation, which is the
        # denominator of `k` and the quantity the measurement guard is applied to.
        z̄ = acc.mean_elevation
        @test all(a -> isapprox(a, _dp_ambient((; z = z̄)) - _DECOUPLING_REFERENCE_TEMPERATURE,
                                atol = 1e-6),
                  decoupling.ambient_excess)

        # All four coefficients are reported, so `k` is re-evaluable off-reference.
        for layer in (:coef_alpha, :coef_beta, :coef_gamma, :coef_delta)
            @test all(isfinite, decoupling[layer])
        end
    end

    @testset "lapse rate" begin
        # The ON-GLACIER lapse rate, which is not the ambient one. Above melting the decoupling maps
        # T -> 273.15 + k(T - 273.15), so it compresses the profile and the slope the fit reports is
        # k*Γ. Positive for cooling with height.
        @test all(g -> isapprox(g, _DP_K_TRUE * _DP_GAMMA_TRUE, atol = 1e-6), lapse.lapse_rate)
        @test DimensionalData.metadata(lapse)["n_fitted"] == length(_DP_TIME)
        @test all(==(4), lapse.n_cells)
        @test collect(dims(lapse, Ti)) == _DP_TIME
        # Elevation spread is the diagnostic that governs whether a slope is determined at all.
        @test all(s -> isapprox(s, Statistics.std([c.z for c in _DP_CELLS]; corrected = false),
                                atol = 1e-6),
                  lapse.elevation_spread)
    end

    @testset "NaN means the forcing carried no information" begin
        # Fewer cells than the four the four-term fit needs. `min_cells` is raised above the cell
        # count rather than the cells removed, so `n_cells` still reports what was available.
        gated = derive_decoupling_factor(acc; min_cells = 5)
        @test all(isnan, gated.decoupling_factor)
        @test all(==(4), gated.n_cells)
        @test DimensionalData.metadata(gated)["n_fitted"] == 0
        # The lapse fit is gated independently, and by the same keyword.
        @test all(isnan, derive_lapse_rate(acc, gated; min_cells = 5).lapse_rate)

        # `min_cells` below 4 does not rescue a three-cell region: four parameters need four cells,
        # and lowering the gate cannot make a rank-deficient system solvable.
        thin = _acc(_DP_CELLS[1:3], _dp_observed)
        thin_k = derive_decoupling_factor(thin; min_cells = 2)
        @test all(isnan, thin_k.decoupling_factor)
        # The slope is still reported, because it only needs two parameters. With `k` unmeasured the
        # fit falls back to `k = 1`, i.e. it regresses the raw temperatures as given without undoing
        # any decoupling — so what comes back is the slope of the reanalysis' own profile, steeper
        # than the ambient Γ here because the glacier fraction rises with elevation in this fixture.
        thin_slope = let z = [c.z for c in _DP_CELLS[1:3]],
                         T = [first(_dp_observed(c)) for c in _DP_CELLS[1:3]]
            -1000 * Statistics.cov(z, T; corrected = false) /
                    Statistics.var(z; corrected = false)
        end
        @test thin_slope > _DP_GAMMA_TRUE
        @test all(g -> isapprox(g, thin_slope, atol = 1e-6),
                  derive_lapse_rate(thin, thin_k; min_cells = 2).lapse_rate)

        # Every cell below melting: no warm excess anywhere, so there is no ratio to form. The
        # ambient guard, not a divide-by-zero, is what catches this.
        frozen = _acc(_DP_CELLS, c -> fill(250.0 - _DP_GAMMA_TRUE * c.z / 1000.0,
                                           length(_DP_TIME)))
        frozen_k = derive_decoupling_factor(frozen; min_cells = 4)
        @test all(isnan, frozen_k.decoupling_factor)
        @test all(isnan, frozen_k.ambient_excess)
        # The slope is still fitted there — with `k = 1` locally, i.e. on the raw temperatures,
        # which is the ambient lapse rate and an honest answer for an unmeasured decoupling.
        @test all(g -> isapprox(g, _DP_GAMMA_TRUE, atol = 1e-6),
                  derive_lapse_rate(frozen, frozen_k; min_cells = 4).lapse_rate)

        # Barely above melting: the fit succeeds but the ambient excess is below the measurement
        # threshold, so `k` is withheld rather than computed from two ~0.1 K numbers.
        marginal = _acc(_DP_CELLS,
                        c -> fill(273.35 - 0.05 * (c.z - 1000.0) / 1000.0, length(_DP_TIME)))
        marginal_k = derive_decoupling_factor(marginal; min_cells = 4)
        @test all(isnan, marginal_k.decoupling_factor)
        @test all(a -> !isfinite(a) || a >= _MIN_AMBIENT_EXCESS, marginal_k.ambient_excess)

        # No elevation spread: every cell at the same reanalysis surface. `glm*z` is then collinear
        # with `glm`, so the decoupling fit is rank-deficient and the slope is undefined.
        flat = _acc([(z = 1500.0, glm = g) for g in (0.0, 0.25, 0.5, 1.0)], _dp_observed)
        flat_k = derive_decoupling_factor(flat; min_cells = 4)
        @test all(isnan, flat_k.decoupling_factor)
        @test all(isnan, derive_lapse_rate(flat, flat_k; min_cells = 4).lapse_rate)
        @test all(isnan, derive_lapse_rate(flat, flat_k; min_cells = 4).elevation_spread)

        # No glm spread: nothing varies the reanalysis's own degree of decoupling, so `k` is not
        # identified even though the elevation slope still is.
        same_glm = _acc([(z = z, glm = 0.4) for z in (1000.0, 1500.0, 2000.0, 2500.0)],
                        _dp_observed)
        same_k = derive_decoupling_factor(same_glm; min_cells = 4)
        @test all(isnan, same_k.decoupling_factor)
        @test all(isfinite, derive_lapse_rate(same_glm, same_k; min_cells = 4).lapse_rate)

        # An accumulator that saw no usable cell at all reports `NaN` rather than throwing.
        empty_acc = (; time = _DP_TIME, n = length(_DP_TIME), sums = nothing,
                     mean_elevation = NaN)
        empty_k = derive_decoupling_factor(empty_acc; min_cells = 4)
        @test all(isnan, empty_k.decoupling_factor)
        @test all(==(0), empty_k.n_cells)
        @test all(isnan, derive_lapse_rate(empty_acc, empty_k; min_cells = 4).lapse_rate)
    end

    @testset "the fits report raw, and a bad k does not reach the slope" begin
        # Warm excess *rising* with glm — the opposite of decoupling. The fit reports what it sees,
        # so `k` lands above 1, outside the domain `climate_adjust_for_glacier` accepts.
        function inverted(c)
            a = _dp_ambient(c)
            return fill(a + c.glm * 0.5 * max(a - _DECOUPLING_REFERENCE_TEMPERATURE, 0.0),
                        length(_DP_TIME))
        end
        inv_acc = _acc(_DP_CELLS, inverted)
        inv_k = derive_decoupling_factor(inv_acc; min_cells = 4)
        @test all(k -> k > _DECOUPLING_FACTOR_LIMITS[2], inv_k.decoupling_factor)
        @test all(isfinite, inv_k.decoupling_factor)

        # That `k` is *not* applied by the lapse fit: outside `(0, 1]` it falls back to `k = 1`, the
        # bit-exact no-op, so the slope is the ambient one rather than an arbitrary rescaling. Same
        # fixture through the ambient temperatures gives the same answer.
        inv_lapse = derive_lapse_rate(inv_acc, inv_k; min_cells = 4)
        raw = derive_lapse_rate(inv_acc, derive_decoupling_factor(_acc(_DP_CELLS, inverted);
                                                                  min_cells = 99);
                                min_cells = 4)
        @test all(isapprox.(inv_lapse.lapse_rate, raw.lapse_rate; atol = 1e-9))
        @test all(g -> _LAPSE_RATE_LIMITS[1] <= g <= _LAPSE_RATE_LIMITS[2], inv_lapse.lapse_rate)
    end

    @testset "a wide spread is reported, not discarded" begin
        # The regression test for the behaviour this replaced: an interquartile threshold used to
        # reject a region's whole fitted series and substitute a constant. Lapse rates and `k` both
        # vary strongly on diurnal and seasonal timescales, so a wide spread is usually physics, and
        # nothing here can tell it from fit error. So every per-timestep fit survives.
        function scatter(c)
            v = _dp_observed(c)
            for t in eachindex(v)
                v[t] += c.glm * 4.0 * (isodd(t) ? 1 : -1)
            end
            return v
        end
        spread = derive_decoupling_factor(_acc(_DP_CELLS, scatter); min_cells = 4)
        k = collect(filter(isfinite, spread.decoupling_factor))
        @test length(k) == length(_DP_TIME)
        # Not collapsed to a constant — the whole point.
        @test length(unique(round.(k; digits = 6))) > 1
        @test Statistics.quantile(k, 0.75) - Statistics.quantile(k, 0.25) > 0.5
        # And still usable: the caller's own median over the raw series.
        @test isfinite(Statistics.median(k))
    end

    @testset "decoupling_factor_at_elevation" begin
        z̄ = DimensionalData.metadata(decoupling)["reference_elevation"]
        # At the reference elevation the re-evaluation must reproduce the original fit exactly —
        # same coefficients, same formula, same guard.
        at_ref = decoupling_factor_at_elevation(decoupling, z̄)
        @test all(t -> isapprox(at_ref[t], decoupling.decoupling_factor[t], atol = 1e-9),
                  eachindex(at_ref))
        @test at_ref isa DimArray
        @test collect(dims(at_ref, Ti)) == _DP_TIME
        @test DimensionalData.metadata(at_ref)["elevation"] == z̄

        # The fixture's ambient excess is linear in elevation by construction, and `k_true` is
        # constant, so `k(z)` is `k_true` at every elevation the ambient excess still supports.
        for z in (1200.0, 1800.0, 2400.0)
            at = decoupling_factor_at_elevation(decoupling, z)
            @test all(k -> isapprox(k, _DP_K_TRUE, atol = 1e-6), at)
            @test DimensionalData.metadata(at)["n_fitted"] == length(_DP_TIME)
            @test DimensionalData.metadata(at)["n_held"] == 0
        end

        # The guard that motivates the function: high enough up, the fitted ambient excess falls
        # below `_MIN_AMBIENT_EXCESS` and `k` stops being a ratio of anything. The fixture's ambient
        # reaches melting at 292.0/5.0 = 3840 m, so at 5000 m every timestep is held — evaluated at
        # the highest elevation its own fit still supports rather than reverted to the identity,
        # which would inject a false gradient at the interval boundary it flips on.
        far = decoupling_factor_at_elevation(decoupling, 5000.0)
        @test DimensionalData.metadata(far)["n_held"] == length(_DP_TIME)
        @test DimensionalData.metadata(far)["n_fitted"] == length(_DP_TIME)
        @test all(k -> isapprox(k, _DP_K_TRUE, atol = 1e-6), far)
        @test DimensionalData.metadata(far)["in_range"]

        # Outside the hypsometry there is no glacier to carry a factor, so the answer is `NaN` and
        # says so — a held value there would be a fabrication rather than an extrapolation.
        out = decoupling_factor_at_elevation(decoupling, 6000.0;
                                             elevation_range = (1000.0, 3000.0))
        @test all(isnan, out)
        @test !DimensionalData.metadata(out)["in_range"]
        @test DimensionalData.metadata(out)["n_fitted"] == 0
        # In range, the range argument changes nothing.
        @test all(isapprox.(decoupling_factor_at_elevation(decoupling, 1800.0;
                                                           elevation_range = (1000.0, 3000.0)),
                            decoupling_factor_at_elevation(decoupling, 1800.0); atol = 1e-12))

        # An unfitted timestep stays unfitted at every elevation.
        gated = decoupling_factor_at_elevation(derive_decoupling_factor(acc; min_cells = 5), 1800.0)
        @test all(isnan, gated)
    end

    @testset "the reference accumulator matches the contract" begin
        S = _downscaling_sums(length(_DP_TIME))
        T(c) = [_dp_ambient(c) + (isodd(t) ? 0.0 : -30.0) for t in eachindex(_DP_TIME)]
        for c in _DP_CELLS
            _accumulate_downscaling_cell!(S, T(c), c.z, c.glm)
        end
        ref = _acc(_DP_CELLS, T).sums
        @test keys(S) == keys(ref)
        for f in keys(ref)
            @test getfield(S, f) == getfield(ref, f)
        end
        # Cold steps contribute to the all-cell sums only.
        @test S.s1 == fill(4, length(_DP_TIME))
        @test S.s1_warm == [isodd(t) ? 4 : 0 for t in eachindex(_DP_TIME)]
        @test_throws DimensionMismatch _accumulate_downscaling_cell!(S, zeros(3), 1000.0, 0.5)
    end

    @testset "truncated excess: the top of the region is below melting" begin
        # Twenty cells from 500 to 2500 m, glacier fraction rising with elevation plus a fixed
        # perturbation, ambient excess 14 − 0.006 z, so every cell above 2333 m is below melting.
        # Fitting only the cells above melting recovers k exactly; fitting the truncated excess as
        # if it were linear would give 0.74.
        k_true = 0.8
        zs = range(500.0, 2500.0; length = 20)
        cells = [(z = z, glm = clamp(clamp((z - 500) / 2000, 0, 1)^2 + 0.05 * sin(3i), 0, 1))
                 for (i, z) in enumerate(zs)]
        E_amb(z) = 14 - 0.006z
        observed(c) = fill(_DECOUPLING_REFERENCE_TEMPERATURE +
                           (E_amb(c.z) > 0 ? (1 - c.glm * (1 - k_true)) * E_amb(c.z) : E_amb(c.z)),
                           length(_DP_TIME))
        acc53 = _acc(cells, observed)
        d = derive_decoupling_factor(acc53)
        @test all(k -> isapprox(k, k_true; atol = 0.01), d.decoupling_factor)
        @test all(==(count(c -> E_amb(c.z) > 0, cells)), d.n_cells_above_melting)
        @test all(==(20), d.n_cells)
        # Too few cells above melting: NaN, even with plenty of cells in the region.
        @test all(isnan, derive_decoupling_factor(acc53; min_cells = 19).decoupling_factor)
    end

    @testset "k is held at the edge of the ambient excess" begin
        # A fit whose k varies with elevation: k(z) = 1 + (γ + δz)/(α + βz), ambient excess
        # α + βz = 10 − 0.004z reaching the 0.5 K threshold at z* = 2375 m.
        n = length(_DP_TIME)
        ti = Ti(_DP_TIME)
        col(v) = DimArray(fill(v, n), ti)
        α, β, γ, δ = 10.0, -0.004, -1.0, -0.0005
        dec = DimStack((decoupling_factor = col(NaN), coef_alpha = col(α), coef_beta = col(β),
                        coef_gamma = col(γ), coef_delta = col(δ)))
        z_star = (_MIN_AMBIENT_EXCESS - α) / β
        k_star = (_MIN_AMBIENT_EXCESS + γ + δ * z_star) / _MIN_AMBIENT_EXCESS
        below = decoupling_factor_at_elevation(dec, 2000.0)
        @test all(k -> k ≈ 1 + (γ + δ * 2000.0) / (α + β * 2000.0), below)
        @test DimensionalData.metadata(below)["n_held"] == 0
        for z in (2400.0, 3000.0, 6000.0)
            above = decoupling_factor_at_elevation(dec, z)
            @test all(==(k_star), above)
            @test DimensionalData.metadata(above)["n_held"] == n
        end
        # Continuous at the edge: approaching z* from below reaches the held value.
        @test only(unique(decoupling_factor_at_elevation(dec, z_star - 1e-6))) ≈ k_star rtol = 1e-6
    end

    @testset "_decouple_per_timestep" begin
        # Per-timestep k, warm and cold steps, one near-saturated warm step. Temperature follows
        # T′ = T + (k − 1)·max(T − T_ref, 0); longwave is unchanged; vapor pressure is capped at
        # saturation of the cooled air and untouched elsewhere.
        T = [283.15, 290.15, 268.15, 278.15]
        k = [0.8, 0.5, 0.3, 1.0]
        e = [800.0, 0.99 * saturation_vapor_pressure(290.15), 300.0, 700.0]
        ti = Ti(_DP_TIME[1:4])
        layer(v) = DimArray(v, ti)
        s = DimStack((temperature_air = layer(T), pressure_air = layer(fill(80000.0, 4)),
                      vapor_pressure = layer(e), wind_speed = layer(fill(4.0, 4)),
                      precipitation = layer(zeros(4)), shortwave_downward = layer(fill(300.0, 4)),
                      longwave_downward = layer([310.0, 330.0, 250.0, 300.0]));
                     metadata = Dict{String,Any}("elevation" => 1500.0))
        d = _decouple_per_timestep(s, k)
        T′ = @. T + (k - 1) * max(T - _DECOUPLING_REFERENCE_TEMPERATURE, 0.0)
        @test parent(d[:temperature_air]) == T′
        @test parent(d[:longwave_downward]) == parent(s[:longwave_downward])
        @test parent(d[:vapor_pressure]) == [800.0, saturation_vapor_pressure(T′[2]), 300.0, 700.0]
        for v in (:pressure_air, :wind_speed, :precipitation, :shortwave_downward)
            @test parent(d[v]) == parent(s[v])
        end
        # Every k = 1 is the identity, returned as the same stack.
        @test _decouple_per_timestep(s, ones(4)) === s
    end

    @testset "_make_applicable" begin
        # Gaps filled and counted; nothing else touched.
        applied, n_filled, n_clamped =
            _make_applicable([0.5, NaN, 0.8], 1.0, _DECOUPLING_FACTOR_LIMITS;
                             open_lower = true, clamp_to_valid_domain = true)
        @test applied == [0.5, 1.0, 0.8]
        @test n_filled == 1
        @test n_clamped == 0

        # Out-of-domain values clamped and counted. `k`'s domain is half-open at zero, so a value at
        # or below the bound is moved just *inside* it — clamping to 0.0 would leave something the
        # consumer still rejects.
        applied, n_filled, n_clamped =
            _make_applicable([-0.2, 0.0, 3.2, 0.6], 1.0, _DECOUPLING_FACTOR_LIMITS;
                             open_lower = true, clamp_to_valid_domain = true)
        @test n_filled == 0
        @test n_clamped == 3
        @test applied[1] == applied[2] == nextfloat(0.0)
        @test applied[1] > 0.0
        @test applied[3] == 1.0
        @test applied[4] == 0.6

        # The lapse rate's domain is closed, so its lower bound is clamped to itself.
        applied, _, n_clamped =
            _make_applicable([-99.0, 6.5, 99.0], 6.5, _LAPSE_RATE_LIMITS;
                             open_lower = false, clamp_to_valid_domain = true)
        @test applied == [_LAPSE_RATE_LIMITS[1], 6.5, _LAPSE_RATE_LIMITS[2]]
        @test n_clamped == 2

        # Opting out fills but does not clamp, and reports zero clamped rather than lying.
        applied, n_filled, n_clamped =
            _make_applicable([3.2, NaN], 1.0, _DECOUPLING_FACTOR_LIMITS;
                             open_lower = true, clamp_to_valid_domain = false)
        @test applied == [3.2, 1.0]
        @test (n_filled, n_clamped) == (1, 0)

        # Returns a plain Vector even from a DimArray layer, which is what the callers index into —
        # `collect` on a stack layer would preserve the `Ti` axis rather than dropping to a Vector.
        applied, = _make_applicable(decoupling.decoupling_factor, 1.0,
                                    _DECOUPLING_FACTOR_LIMITS;
                                    open_lower = true, clamp_to_valid_domain = true)
        @test applied isa Vector{Float64}
        @test length(applied) == length(_DP_TIME)
    end

    @testset "_solve_sym4 and _r2" begin
        # Identity system: the solution is the right-hand side.
        @test all(isapprox.(_solve_sym4((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0),
                                        (1.0, 2.0, 3.0, 4.0)),
                            (1.0, 2.0, 3.0, 4.0)))
        # Rank-deficient (a zero pivot) returns `nothing` rather than a garbage fit: Cholesky is its
        # own rank test, which is why the fits can rely on it to catch collinear regressors.
        @test _solve_sym4((1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0),
                          (1.0, 2.0, 3.0, 4.0)) === nothing
        # A non-finite or non-positive diagonal scale is rejected up front.
        @test _solve_sym4((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          (1.0, 1.0, 1.0, 1.0)) === nothing

        # `_r2` is 1 for a perfect fit and `NaN` where the response has no variance to explain.
        @test _r2(4, 8.0, 20.0, 20.0) == 1.0
        @test isnan(_r2(0, 0.0, 0.0, 0.0))
        @test isnan(_r2(4, 8.0, 16.0, 16.0))     # SS_tot == 0
        # Clamped into [0, 1] so accumulated roundoff cannot report a nonsensical value.
        @test _r2(4, 8.0, 20.0, 100.0) == 1.0
        @test _r2(4, 8.0, 20.0, 0.0) == 0.0
    end
end
