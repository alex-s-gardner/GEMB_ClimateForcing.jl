# Downscaling parameters: what a region's reanalysis forcing says about its own decoupling and
# lapse rate.
#
# Two fits, each a *relationship across* the grid cells of a region rather than a property of any
# one of them:
#
#   1. the local on-glacier temperature decoupling factor `k`, from how the warm excess varies with
#      the reanalysis glacier mask fraction at matched elevation ([`derive_decoupling_factor`](@ref));
#   2. the local on-glacier temperature lapse rate, from how decoupled temperature varies with
#      elevation ([`derive_lapse_rate`](@ref)).
#
# Both *measure*: they report what the region's forcing supports at the forcing's own time
# resolution, and `NaN` where it supports nothing. They substitute nothing and clamp nothing,
# because both quantities are validated by the functions that consume them
# (`climate_adjust_for_glacier` rejects `k` outside `(0, 1]`, `climate_adjust_for_elevation` rejects
# a lapse rate outside `_LAPSE_RATE_LIMITS`), and because deciding what a gap or a wide spread
# *means* is the caller's business, not this file's. Downstream grouping is a `groupby` over the
# returned `Ti` axis.
#
# These invert the corrections in `glacier_adjustment.jl` and `elevation_adjustment.jl`: `k` is read
# back out of a forcing record that has already had the reanalysis's own decoupling baked into it,
# and the lapse rate out of temperatures brought to a common on-glacier state. That is why they live
# here rather than with whatever assembles a region — the arithmetic is entirely about this
# package's own adjustment conventions.
#
# ## The `acc` argument
#
# Both fits take an accumulator of per-timestep cross-cell regression sums rather than a stack of
# forcing, so a region of hundreds of cells over decades of hourly steps never has to be held in
# memory at once. Its producer lives in the caller (GEMB_GlacierSims builds one by streaming a
# glacier elevation-class table), so the contract is documented here rather than being visible in
# one file:
#
#   acc.n              — number of timesteps
#   acc.time           — the `Ti` axis values, `length(acc.time) == acc.n`
#   acc.mean_elevation — area-weighted mean of the contributing cells' reanalysis elevations (m),
#                        the reference `k` is reported at
#   acc.sums           — the sums below, or `nothing` if no cell had usable forcing
#
# `acc.sums` holds one length-`n` vector per field. With the design vector `x = [1, z, glm, glm*z]`
# per cell — `z` its reanalysis elevation (m), `glm` its glacier-mask fraction — and `w = glm*z`:
#
#   s1                       — count of contributing cells (Int)
#   z, g, w                  — Σz, Σglm, Σw
#   zz, zg, zw, gg, gw, ww   — the upper triangle of XᵀX
#   E, zE, gE, wE            — Xᵀ against the warm excess `E = max(T - T_ref, 0)`
#   EE                       — ΣE², for the decoupling fit's `r2`
#   T, zT                    — ΣT and ΣzT, for the lapse fit
#
# Duck-typed on purpose: a NamedTuple of `Vector{Float64}` (plus `Vector{Int}` for `s1`) is all
# either fit reads, and pinning a struct here would make the accumulator's own layout this
# package's business.

# Melting point (K). The decoupling correction damps only the excess above it, matching
# `climate_adjust_for_glacier`'s `reference_temperature` default.
const _DECOUPLING_REFERENCE_TEMPERATURE = 273.15

# Smallest ambient warm excess (K) that can support a decoupling fit.
#
# This is a *measurement* threshold, not a policy one: `k` is a ratio of two fitted excesses, so a
# timestep where the region is barely above melting divides two ~0.1 K numbers and returns noise. On
# real Alpine forcing those winter timesteps fit wildly and produced a confident-looking seasonal
# `k` of 0.2 that was pure fit error. A timestep below this threshold carries no information about
# `k`, so it reports `NaN` — which is a different claim from "the fit is unreasonable", and the only
# kind of claim this file makes. 0.5 K is well above the fit residual while still admitting every
# melt-season timestep, which is where the decoupling matters at all.
const _MIN_AMBIENT_EXCESS = 0.5

# Fewest usable cells a region needs before either fit is attempted. See `derive_decoupling_factor`
# for the subsampling that picks 8. Defined once so the two fits and their callers cannot drift into
# gating different regions.
const _MIN_CELLS_DEFAULT = 8

# Physically plausible lapse-rate bounds (K/km) — exactly the bounds
# `climate_adjust_for_elevation` validates against and throws outside (`elevation_adjustment.jl`),
# named here so a clamped value is *guaranteed* acceptable there.
#
# The *fit* does not clamp to these — a slope of 200 K/km is real evidence that the region does not
# constrain a lapse rate at that timestep, and clamping to 25 would hide exactly that. They are the
# domain a *consumer* clamps into when it applies a fitted series, since that is where a value has
# to be acceptable.
const _LAPSE_RATE_LIMITS = (_LAPSE_RATE_MIN, _LAPSE_RATE_MAX)

# The domain `climate_adjust_for_glacier` accepts for `k`, half-open at zero: 1 is no damping and
# anything above it warms the surface instead, while 0 would erase the warm excess entirely. Same
# division of labour as `_LAPSE_RATE_LIMITS` — the fit reports outside it, the application clamps
# into it.
const _DECOUPLING_FACTOR_LIMITS = (0.0, 1.0)

# The seven forcing layers `climate_forcing` returns. Named here rather than derived from a stack so
# a caller accumulating layer-wise can allocate before it has seen one.
const _FORCING_VARIABLES = (:temperature_air, :pressure_air, :vapor_pressure, :wind_speed,
                            :precipitation, :shortwave_downward, :longwave_downward)

"""
    derive_decoupling_factor(acc; min_cells = $(_MIN_CELLS_DEFAULT)) -> DimStack

Fit the on-glacier temperature decoupling factor `k` per timestep from cross-cell sums, as
accumulated over a region's grid cells.

`climate_adjust_for_glacier` damps the warm excess above melting, `T′ = T + (k - 1)·max(T - T_ref, 0)`,
so `E′ = k·E`. The estimator inverts that using the one place the reanalysis already varies its own
degree of decoupling: `glm`, the fraction of each cell its land-surface scheme treats as ice.

**Which end of `glm` is decoupled matters, and it is the high end.** A cell with `glm = 1` is one
ERA5-Land already runs as ice, so its near-surface temperature is *already* damped and needs no
correction — which is precisely why `GEMB_GlacierSims.cell_decoupling_factor` weights the correction it
applies by `1 - glm`. A cell with `glm = 0` is run as ice-free ground and its temperature is
ambient, so it needs the full correction. Mixing the two, a cell's observed excess is

    E_i = [1 - glm_i·(1 - k)] · E_ambient(z_i)

falling from `E_ambient` at `glm = 0` to `k·E_ambient` at `glm = 1`. So `k` is read off as the ratio
of the fitted excess at the two ends of `glm`, and the warm excess must *decrease* with `glm` for the
fit to mean anything.

Note the *product*: the `glm` weighting scales an ambient excess that itself varies with elevation.
Taking `E_ambient` locally linear in `z`, `E_ambient ≈ a + b·z`, expands that product into

    E_i = [1 - glm_i(1 - k)]·(a + b·z_i)
        = α + β·z_i + γ·glm_i + δ·(glm_i·z_i)

with `α = a`, `β = b`, `γ = -(1-k)a`, `δ = -(1-k)b`. So the fit is on four terms, and the `glm × z`
interaction is not optional — dropping it biases `k` toward 1, because the `glm` slope then has to
absorb part of the elevation dependence of the very quantity it multiplies. Evaluating at the
region's area-weighted mean elevation `z̄`:

    k = [(α + γ) + (β + δ)·z̄] / (α + β·z̄)

i.e. the fitted excess of a fully-glaciated cell over that of an ice-free one at the same elevation.

Elevation is a control, not a nuisance: `glm` and elevation covary strongly (high cells are more
glaciated), so without it the `glm` slope is partly a topographic signal.

Returns a `DimStack` on the accumulated `Ti` axis:

| layer | meaning |
|:--|:--|
| `decoupling_factor` | `k` at the reference elevation, raw. `NaN` where the region supported no fit |
| `n_cells` | grid cells contributing to that timestep |
| `r2` | coefficient of determination of the four-term fit |
| `ambient_excess` | the fitted ice-free excess `α + β·z̄` (K) — what the fit's precision hinges on |
| `coef_alpha`, `coef_beta`, `coef_gamma`, `coef_delta` | all four coefficients, so `k` is re-evaluable at any elevation — see [`decoupling_factor_at_elevation`](@ref) |

`metadata` carries `reference_elevation` (`z̄`), `min_cells`, `n_fitted`, and `n_timesteps`.

**`NaN` means "the forcing carried no information here", and nothing more.** A timestep reports `NaN`
when it has fewer than `max(min_cells, 4)` cells (four parameters need four cells), no spread in `glm`
or elevation, collinear regressors, or an ambient excess below `$(_MIN_AMBIENT_EXCESS)` K. That last
case is most of a glacier record and is a real identifiability limit rather than a divide-by-zero
guard: `k` is a *ratio* of fitted excesses, so a region sitting near melting fits one out of noise —
on real Alpine forcing that produced a confident-looking winter `k` of 0.2 from a 0.1 K excess.

**Nothing is substituted and nothing is clamped.** A fit of 3.2 is reported as 3.2, even though
`climate_adjust_for_glacier` would reject it: that value is evidence about the fit, and clamping it to
1.0 would silently convert a diagnostic into a plausible-looking measurement. Filling and clamping
happen only where a series is *applied*, in
`GEMB_GlacierSims.elevation_interval_forcing`, which is the one consumer that cannot
accept a `NaN`.

The per-timestep series is legitimately broad — each value is a ratio of two regression outputs, so on
real Alpine forcing the hourly `k` spans -1.3 to +3.2 between the 1st and 99th percentiles. Summarise
it with a **median**, not a mean; the tail is one-sided and drags a mean badly (against a synthetic
truth of 0.930 with 0.6 K of noise, the mean returns 0.875 and the median 0.926). Grouping is a
`groupby` over the `Ti` axis, so the resolution is the caller's choice:

```julia
using Statistics
k = derive_decoupling_factor(acc)
median(filter(isfinite, k.decoupling_factor))                                    # pooled
map(g -> median(filter(isfinite, g.decoupling_factor)), groupby(k, Ti => month)) # by month
```

`ambient_excess` is the `α + β·z̄` profile [`derive_lapse_rate`](@ref) needs in order to undo the
reanalysis's own partial decoupling cell by cell, and it is the honest measure of how well determined
a given timestep's `k` is: the closer it sits to `$(_MIN_AMBIENT_EXCESS)` K, the more of the ratio is
noise.

## Screening `k` downstream, and why it also screens the lapse rate

`k` is the ratio `[(α+γ) + (β+δ)·z̄] / (α + β·z̄)`, whose **denominator is exactly
`ambient_excess`**. `$(_MIN_AMBIENT_EXCESS)` K is only enough to keep the ratio finite, not enough to
make it meaningful, so an out-of-domain `k` is the expected outcome near the threshold rather than a
sign of a broken fit — `r2` is a median 0.97 for those timesteps, statistically indistinguishable from
the well-behaved ones. `r2` cannot be used to screen them; `ambient_excess` can.

Two years of ERA5-Land at Wrangell-St Elias (38 cells, 3427 timesteps with a measurable `k`):

| `_MIN_AMBIENT_EXCESS` | timesteps kept | with `k` outside `(0,1]` | max `k` |
|--:|--:|--:|--:|
| 0.5 (current) | 3427 (100%) | 788 (23.0%) | 7.35 |
| 1.0 | 3174 (93%) | 634 (20.0%) | 4.33 |
| 2.0 | 2770 (81%) | 456 (16.5%) | 2.63 |
| 3.0 | 2346 (68%) | 287 (12.2%) | 2.00 |

Raising the threshold trades coverage for tail control and never eliminates the tail, which is why it
stays at $(_MIN_AMBIENT_EXCESS) K here — the measurement guard — and the rest is a post-processing
decision. **Screen on `k ∈ (0, 1]` directly**, which is what
`climate_adjust_for_glacier` requires anyway, and drop or down-weight `lapse_rate` at the same
timesteps: those `k` values corrupt the lapse fit that consumes them, and the corrupted slopes are not
detectable from the slope series alone. See the warning in [`derive_lapse_rate`](@ref) for the
measurements.

Neither the narrowness nor the spread of the glacier-mask fraction is the cause. `glm` is **static per
cell**, so it is identical at every timestep and cannot explain why 788 of 17521 blow up while the
rest do not. At Wrangell it spans 0.002–1.0 across the 38 cells (sd 0.38, 7 cells below 0.1 and 13
above 0.9) — a fully sampled axis. What is true is that `glm` correlates with elevation
(`r = 0.92` here, VIF 20 for `glm` and 46 for the `glm × z` interaction), so the four coefficients are
individually poorly determined even where the fitted surface is excellent. That is tolerable for `k`
at `z̄`, which depends on the well-determined combinations, and it is what makes `ambient_excess` — not
the coefficients — the diagnostic to screen on.

## Why `min_cells` defaults to 8

Four cells make the four-term fit *determinate*, which is not the same as identified. Because forcing
completeness is screened per cell, the contributing count is constant across time — so `min_cells`
gates whole regions, and the cost of raising it is measured in regions dropped rather than timesteps
pruned.

Subsampling Wrangell's 38 cells to `m`, spread through the elevation ordering, and scoring each
subsample's `k` against the all-38-cell fit:

| `m` | median \\|k − k₃₈\\| | `k` outside `(0,1]` | timesteps fitted |
|--:|--:|--:|--:|
| 4 | 0.253 | 46.3% | 15.1% |
| 5 | 0.538 | 75.1% | 11.6% |
| 6 | 0.200 | 46.7% | 16.6% |
| 8 | 0.141 | 41.9% | 16.2% |
| 12 | 0.105 | 31.9% | 20.1% |
| 20 | 0.059 | 43.4% | 18.6% |
| 30 | 0.041 | 39.7% | 19.9% |

At four to six cells the error on `k` is 0.20–0.25 against a quantity whose whole range is `(0, 1]` —
a quarter of the scale, and one draw at `m = 5` put 75% of its fits outside the valid domain
altogether. The error falls steeply to `m ≈ 8–12` and then flattens; past 20 cells there is little left
to gain. The lapse fit, having two parameters instead of four, is far less demanding: its median error
is already 0.40 K/km at `m = 4` and 0.38 at `m = 8`, against slopes of ~5.5 K/km. So the four-parameter
fit sets the floor and both share it.

Against global coverage in the ERA5-Land elevation-class table (1912 one-degree boxes containing
glacier cells, 697,693 km² total; box cell counts have a median of 15 and a first decile of 1):

| `min_cells` | boxes retained | glacier area retained |
|--:|--:|--:|
| 3 | 84.5% | 99.82% |
| 4 | 79.7% | 99.67% |
| 8 | 66.1% | 98.49% |
| 12 | 56.1% | 96.50% |
| 20 | 43.4% | 91.52% |

8 costs 1.5% of global glacier area while removing the regime where `k` is worse than a coin flip; 12
would buy a little more precision for 3.5%. The boxes it drops are the small ones — a 3-cell box was
never going to identify a four-term surface, and reporting `NaN` for it is the honest outcome. Regions
below the threshold are not silently degraded: every timestep reports `NaN`, and `n_cells` in the
returned stack says why. Raise it toward 12 if precision matters more than coverage; lower it only
with `elevation_spread` and the out-of-domain `k` fraction in hand.
"""
function derive_decoupling_factor(acc; min_cells::Int = _MIN_CELLS_DEFAULT)
    n = acc.n
    z̄ = acc.mean_elevation
    ti = Ti(acc.time)

    # `NaN` is the starting value, not a fallback: a timestep is only overwritten by a fit that
    # survives the measurement guards below, and reporting `NaN` where none did is the whole
    # contract. Insufficient information is the *common* case, not an edge case — on real Alpine
    # forcing only 44% of hourly timesteps carry enough ambient excess to fit at all.
    decoupling_factor = fill(NaN, n)
    r2 = fill(NaN, n)
    n_cells = zeros(Int, n)
    # The fitted ice-free excess at `z̄`, per timestep. `NaN` where nothing was fitted.
    ambient_excess = fill(NaN, n)
    # All four fitted coefficients per timestep, so `k` can be re-evaluated at an elevation other
    # than `z̄` — which is what every elevation interval needs, since `z̄` is a reanalysis surface
    # and the intervals are the glacier. See `decoupling_factor_at_elevation`.
    cα, cβ, cγ, cδ = fill(NaN, n), fill(NaN, n), fill(NaN, n), fill(NaN, n)

    S = acc.sums
    if S !== nothing
        @inbounds for t in 1:n
            n_cells[t] = S.s1[t]
            # Four parameters need four cells before the system is even determined.
            S.s1[t] >= max(min_cells, 4) || continue

            # XᵀX for x = [1, z, glm, glm*z], upper triangle, and Xᵀy for y = E.
            A = (Float64(S.s1[t]), S.z[t], S.g[t],  S.w[t],
                                   S.zz[t], S.zg[t], S.zw[t],
                                            S.gg[t], S.gw[t],
                                                     S.ww[t])
            coef = _solve_sym4(A, (S.E[t], S.zE[t], S.gE[t], S.wE[t]))
            coef === nothing && continue
            α, β, γ, δ = coef

            # The fitted excess at the reference elevation, at each end of `glm`: an ice-free cell
            # (`glm = 0`), which the reanalysis leaves ambient, and a fully-glaciated one
            # (`glm = 1`), which it has already decoupled.
            ambient = α + β * z̄
            decoupled = ambient + (γ + δ * z̄)
            # `k` is a ratio of fitted excesses, so it is only measurable when the denominator is
            # meaningfully above zero. A region at or near melting carries no information about `k`
            # — and, worse, fits one anyway out of noise. Report nothing instead.
            ambient >= _MIN_AMBIENT_EXCESS || continue

            # Raw, and stays raw. A ratio of two regression outputs is legitimately broad; clamping
            # here would hide precisely the timesteps a caller needs to see before trusting a
            # summary. `elevation_interval_forcing` clamps, because it applies.
            decoupling_factor[t] = decoupled / ambient
            r2[t] = _r2(S.s1[t], S.E[t], S.EE[t],
                        α * S.E[t] + β * S.zE[t] + γ * S.gE[t] + δ * S.wE[t])
            ambient_excess[t] = ambient
            cα[t] = α; cβ[t] = β; cγ[t] = γ; cδ[t] = δ
        end
    end

    metadata = Dict{String,Any}("reference_elevation" => z̄,
                                "min_cells" => min_cells,
                                "n_fitted" => count(isfinite, decoupling_factor),
                                "n_timesteps" => n)

    return DimStack((decoupling_factor = DimArray(decoupling_factor, ti),
                     n_cells = DimArray(n_cells, ti),
                     r2 = DimArray(r2, ti),
                     ambient_excess = DimArray(ambient_excess, ti),
                     coef_alpha = DimArray(cα, ti),
                     coef_beta = DimArray(cβ, ti),
                     coef_gamma = DimArray(cγ, ti),
                     coef_delta = DimArray(cδ, ti)); metadata)
end

"""
    decoupling_factor_at_elevation(decoupling, elevation; elevation_range = nothing)
        -> DimArray

Re-evaluate the fitted decoupling factor at `elevation` rather than at the region's area-weighted
mean elevation. Returns a raw `k(z)` series on the same `Ti` axis, `NaN` at every timestep the fit
cannot support one, with `n_held` / `in_range` / `elevation` in its metadata.

[`derive_decoupling_factor`](@ref) reports `k` at `z̄`, the area-weighted mean of the *reanalysis*
surfaces it fit from. But the elevation intervals the forcing is actually needed at are the glacier,
which sits systematically above that — on real Alpine forcing 76% of the glacier area is above every
contributing cell. Since the fit is on four terms including the `glm × z` interaction, `k` is a
function of elevation and can be evaluated anywhere:

    k(z) = [(α + γ) + (β + δ)·z] / (α + β·z)

**The guard that protects `k` at `z̄` does not protect `k(z)`.** The denominator is the ambient warm
excess `E_ambient(z) = α + β·z`, which falls with height and crosses zero somewhere above the
region. `derive_decoupling_factor` checks it only at `z̄`, so a timestep can pass there and still have
no ambient excess left at an interval 1500 m higher, where the ratio explodes and then changes sign.

That ceiling is a real limit, not a formality. Solved per timestep on ERA5-Land, its median lands at
3750 m in the Ötztal (no glacier area above it) but at 5994 m in the Khumbu, where **43% of the
glacier area** sits higher — so in High Mountain Asia most of the ice is above the elevation at
which the reanalysis carries any information about decoupling at all.

**Above the ceiling `k` is held, not defaulted.** Since `E_ambient` is linear in `z`, the set of
elevations where it clears `$(_MIN_AMBIENT_EXCESS)` K is a half-line, and a fitted timestep always
has `z̄` inside it. So an `elevation` outside that half-line is evaluated at its nearest edge
instead: `k` flattens off above the ceiling at the last value the fit supports. Reverting to the
identity there would instead put a discontinuity into the profile — at Khumbu, a jump from 0.61 to
1.0 across a single interval boundary, with 43% of the ice on the far side of it — and that gradient
is exactly what the downstream sweep interprets as a mass-balance signal. `n_held` counts the
timesteps held this way; when it approaches the number of finite values the reported `k` is the
ceiling value rather than a fit at this elevation.

`elevation_range` is the range `k` is meaningful over — pass the region's hypsometry range and an
`elevation` outside it returns all `NaN` with `in_range = false`, since there is no glacier there to
carry a decoupling factor and a held value would be a fabrication rather than an extrapolation. Left
as `nothing`, every elevation is in range. `NaN` is safe here in a way it is not where a series is
applied: nothing outside the hypsometry becomes an elevation interval, so no `NaN` reaches a forcing
stack from this path.
"""
function decoupling_factor_at_elevation(decoupling, elevation::Real; elevation_range = nothing)
    z = Float64(elevation)
    α, β, γ, δ = decoupling.coef_alpha, decoupling.coef_beta,
                 decoupling.coef_gamma, decoupling.coef_delta
    ti = dims(decoupling, Ti)
    n = length(ti)

    # Outside the hypsometry there is no glacier, so there is nothing to report a factor for.
    in_range = elevation_range === nothing ||
               (Float64(elevation_range[1]) <= z <= Float64(elevation_range[2]))
    if !in_range
        return DimArray(fill(NaN, n), ti;
                        metadata = Dict{String,Any}("elevation" => z, "in_range" => false,
                                                    "n_held" => 0, "n_fitted" => 0))
    end

    k = fill(NaN, n)
    n_held = 0
    @inbounds for t in 1:n
        isfinite(α[t]) || continue
        a, b = α[t], β[t]
        # Evaluate at `z` if the ambient excess still clears the threshold there, and otherwise at
        # the nearest elevation where it does. `E_ambient = a + b·z` is linear, so its valid set is a
        # half-line whose edge is `z* = (threshold - a) / b`; a fitted timestep puts `z̄` inside it,
        # so the edge is always on the far side of `z̄` from `z` and clamping toward it cannot
        # overshoot.
        z_eval = z
        if a + b * z < _MIN_AMBIENT_EXCESS
            b == 0 && continue                     # constant and already below: nothing to hold onto
            z_star = (_MIN_AMBIENT_EXCESS - a) / b
            z_eval = b < 0 ? min(z, z_star) : max(z, z_star)
            n_held += 1
        end
        ambient = a + b * z_eval
        ambient >= _MIN_AMBIENT_EXCESS || continue
        k[t] = (ambient + (γ[t] + δ[t] * z_eval)) / ambient
    end

    return DimArray(k, ti;
                    metadata = Dict{String,Any}("elevation" => z, "in_range" => true,
                                                "n_held" => n_held,
                                                "n_fitted" => count(isfinite, k)))
end

"""
    derive_lapse_rate(acc, decoupling; min_cells = $(_MIN_CELLS_DEFAULT)) -> DimStack

Fit the **on-glacier** temperature lapse rate (K/km, positive for cooling with height) per timestep
from cross-cell sums, after decoupling each cell's temperature with `decoupling_factor`.

Decoupling comes first because the lapse rate wanted here is the on-glacier one and not the ambient
one: `k` damps the warm excess, which compresses the near-surface temperature range over melting ice
and so flattens the slope against elevation. Fitting the raw reanalysis temperature would give the
ambient lapse rate, which is not what a glacier surface sees.

The correction applied to each cell is the *effective* factor, weighted by that cell's non-glacier
fraction exactly as `GEMB_GlacierSims.cell_decoupling_factor` weights it:

    T′_i = T_i + (k_eff_i - 1)·E_i,   k_eff_i = 1 - (1 - k)(1 - glm_i)
         = T_i + (k - 1)(1 - glm_i)·E_i

A cell ERA5-Land already treats as fully glaciated (`glm = 1`) therefore gets no correction, and an
ice-free one gets the full correction, so every cell is brought to a common on-glacier state before
the slope is taken. Applying one uniform `k` to all of them instead would over-damp the glaciated
cells — which are systematically the high ones — and so tilt the fitted slope.

`E_i` here is the **ambient** excess at cell `i`'s elevation — what an ice-free cell there would see
— not the cell's observed excess, which is already damped by its own `glm`. That is exactly the
ice-free end of the decoupling fit, `α + β·z`, so the `(α, β)` per-timestep coefficients from
`decoupling` are reused rather than re-derived. The whole correction is then linear in the sums
already accumulated: no second pass over the forcing.

`decoupling` is the `DimStack` [`derive_decoupling_factor`](@ref) returns, read **raw**. A timestep
whose `k` is `NaN` *or* outside `(0, 1]` is corrected with `k = 1` *locally* — the bit-exact no-op,
which is the right treatment because a timestep with no measurable decoupling is one where the
reanalysis warm excess is too small to matter to the slope anyway. Nothing substituted is written back
into either output.

!!! warning "A corrupt `k` corrupts this slope, and the slope alone does not show it"
    The correction above is proportional to `k - 1`, so an out-of-domain `k` scales every cell's
    temperature by an arbitrary amount before the slope is taken — and because the correction carries
    the `(1 - glm_i)` weight, and `glm` rises with elevation in most regions, the damage is *not*
    scattered noise a robust estimator could reject. It is a smooth, near-collinear false gradient
    that any slope estimator fits confidently.

    Measured on two years of ERA5-Land at Wrangell-St Elias (38 cells, 17521 timesteps):

    | `k` at that timestep | fitted lapse rate, median | max |
    |:--|--:|--:|
    | `NaN` (no correction applied) | 5.2 | 8.95 |
    | inside `(0, 1]` | 3.7 | 8.09 |
    | outside `(0, 1]` | 6.88 | **46.09** |

    Every timestep exceeding $(_LAPSE_RATE_LIMITS[2]) K/km came from a `k` outside `(0, 1]`; none
    came from a `k` inside it. `corr(|k - 1|, |Γ - median Γ|) = 0.81`. At the worst timestep
    `k = 7.35` inflated the lowest cell's temperature from 287.5 K to 375.8 K, and the fitted slope
    went from 7.0 K/km on the raw temperatures to 46.1 K/km on the "corrected" ones. Ordinary least
    squares, Theil–Sen, and Huber all return 46.1 ± 0.3 there: **the estimator is not the problem and
    a robust fit does not help.**

    So the correction is applied **only where `k` lands inside `(0, 1]`**. Outside it — as with a `NaN`
    — the timestep is corrected with `k = 1`, the bit-exact no-op, and the slope is taken on the raw
    temperatures. That is the ambient lapse rate rather than the on-glacier one, which is the honest
    answer for a timestep whose decoupling was never measured: the reported slope is then a real
    measurement of something slightly different, rather than a confident measurement of nothing. At
    Wrangell this drops the series maximum from 46.09 to 8.95 K/km and brings every timestep inside
    $(_LAPSE_RATE_LIMITS) K/km, so the clamp in
    `GEMB_GlacierSims.elevation_interval_forcing` no longer engages on the slope at
    all.

    The screen is on `k`, not on the fitted slope — a corrupted slope is not identifiable from the
    slope series. `decoupling_factor` itself is still reported raw, 7.35 and all, since that value is
    evidence about the decoupling fit. To know which timesteps were corrected, test
    `0 < decoupling_factor <= 1` on the companion stack.

Returns a `DimStack` on the accumulated `Ti` axis:

| layer | meaning |
|:--|:--|
| `lapse_rate` | K/km, **positive for cooling with height**, raw. `NaN` where the region supported no fit |
| `n_cells` | grid cells contributing to that timestep |
| `elevation_spread` | standard deviation (m) of the contributing cells' reanalysis elevations |

`metadata` carries `min_cells`, `n_fitted`, and `n_timesteps`.

`elevation_spread` is the diagnostic that matters: it is what actually governs whether a slope is
determined at all, far more than the cell count. A timestep reports `NaN` when it has fewer than
`min_cells` cells or no spread in elevation.

Two parameters rather than four makes this the less demanding of the two fits — subsampling Wrangell,
its median error against the all-cell slope is 0.40 K/km at 4 cells and 0.38 at 8, on slopes of
~5.5 K/km — so `min_cells` is set by the decoupling fit's needs, not this one's. See
[`derive_decoupling_factor`](@ref) for the numbers.

**Nothing is clamped.** A fitted slope of 200 K/km is real evidence that this region does not
constrain a lapse rate at that timestep, and clamping it into $(_LAPSE_RATE_LIMITS) K/km would hide
exactly that. The clamp lives in
`GEMB_GlacierSims.elevation_interval_forcing`, the one place a slope is applied and
so the one place it has to be acceptable to `climate_adjust_for_elevation`.

Summarise with a **median**, not a mean, and expect real spread: a per-timestep slope over a handful
of cells is broad, and much of that breadth is physical — nocturnal valley inversions and steep
midday cooling are both in the record. At Wrangell the fits span -0.28 to +7.84 K/km between the 5th
and 95th percentiles. Group over the `Ti` axis at whatever resolution the decision needs:

```julia
using Statistics, Dates
Γ = derive_lapse_rate(acc, decoupling)
map(g -> median(filter(isfinite, g.lapse_rate)), groupby(Γ, Ti => month))  # seasonal cycle
map(g -> median(filter(isfinite, g.lapse_rate)), groupby(Γ, Ti => hour))   # diurnal cycle
```
"""
function derive_lapse_rate(acc, decoupling; min_cells::Int = _MIN_CELLS_DEFAULT)
    n = acc.n
    ti = Ti(acc.time)

    # `NaN` where no slope was fitted — nothing substituted, as with `k`.
    lapse_rate = fill(NaN, n)
    elevation_spread = fill(NaN, n)
    n_cells = zeros(Int, n)

    k = decoupling.decoupling_factor
    cα, cβ = decoupling.coef_alpha, decoupling.coef_beta

    S = acc.sums
    if S !== nothing
        @inbounds for t in 1:n
            n_cells[t] = S.s1[t]
            S.s1[t] >= min_cells || continue

            m = Float64(S.s1[t])

            # Elevation spread. Zero means every cell in the region sits at the same reanalysis
            # elevation, so no slope is identifiable however many cells there are.
            Szz = S.zz[t] - S.z[t] * S.z[t] / m
            Szz > 0 || continue
            elevation_spread[t] = sqrt(Szz / m)

            # Bring every cell to the same fully-on-glacier state before taking the slope. The
            # remaining correction for cell `i` is what is left of `k` after the reanalysis has
            # already applied its own `glm`-worth of it:
            #
            #     T_glacier_i = T_i + (k - 1)(1 - glm_i)·E_ambient(z_i)
            #
            # `E_ambient` is the ambient excess — the excess an ice-free cell at that elevation would
            # see — which is *not* the observed excess, since the observed one is already damped by
            # `glm_i`. It is exactly the ice-free end of the decoupling fit, `α + β·z`, so the
            # coefficients that fit came out of are reused here rather than re-derived.
            #
            # An unmeasurable timestep is corrected with `k = 1` here and only here: `c = 0` makes
            # every correction term vanish and the slope is taken on the raw temperatures, which is
            # the bit-exact no-op. Nothing is written back into the decoupling series.
            #
            # A `k` outside `(0, 1]` is treated the same way, and for the same reason: it is not a
            # decoupling factor, it is the ratio of two near-zero fitted excesses. Applying it scales
            # every cell by an arbitrary amount, and because the correction carries `(1 - glm_i)` and
            # `glm` rises with elevation, the damage is a smooth false gradient rather than scatter —
            # so it lands in the slope in full and no estimator can see it. Measured at Wrangell it
            # took a 7.0 K/km slope to 46.1 K/km. See the warning in the docstring.
            fit_ok = isfinite(k[t]) && 0.0 < k[t] <= 1.0 &&
                     isfinite(cα[t]) && isfinite(cβ[t])
            c = fit_ok ? k[t] - 1.0 : 0.0
            α, β = fit_ok ? (cα[t], cβ[t]) : (0.0, 0.0)

            # Σ(1-glm)E_ambient and Σz(1-glm)E_ambient, from the stored regressor sums: Σz(1-glm) is
            # `z - w` and Σz²(1-glm) is `zz - zw`, `w` being the `glm·z` interaction regressor.
            s_corr  = α * (m - S.g[t]) + β * (S.z[t] - S.w[t])
            sz_corr = α * (S.z[t] - S.w[t]) + β * (S.zz[t] - S.zw[t])

            sT  = S.T[t]  + c * s_corr
            szT = S.zT[t] + c * sz_corr

            # Simple linear regression of the on-glacier temperature on elevation.
            SzT = szT - S.z[t] * sT / m

            # K/m -> K/km, sign flipped so positive means cooling with height (the convention
            # `climate_adjust_for_elevation` expects). Raw — see the docstring on why.
            lapse_rate[t] = -1000.0 * (SzT / Szz)
        end
    end

    metadata = Dict{String,Any}("min_cells" => min_cells,
                                "n_fitted" => count(isfinite, lapse_rate),
                                "n_timesteps" => n)

    return DimStack((lapse_rate = DimArray(lapse_rate, ti),
                     n_cells = DimArray(n_cells, ti),
                     elevation_spread = DimArray(elevation_spread, ti)); metadata)
end

"""
    _solve_sym4(A, b) -> NTuple{4,Float64} or nothing

Solve the symmetric positive-semidefinite 4x4 system `A x = b` by Cholesky, where `A` is the upper
triangle in row-major order — `(a11, a12, a13, a14, a22, a23, a24, a33, a34, a44)`.

These are normal equations, so `A = XᵀX` is symmetric and PSD by construction and Cholesky is both
the cheapest factorization and its own rank test: a non-positive pivot means the regressors are
collinear, which returns `nothing` rather than a garbage fit. That happens for real regions — every
cell at the same elevation, every cell with the same glacier fraction, or a region small enough that
`glm·z` is indistinguishable from `z`.

The pivot test is scale-relative. Elevation is in metres, so `Σz²` runs to 1e8 per cell and an
absolute tolerance would reject well-conditioned systems.
"""
function _solve_sym4(A::NTuple{10,Float64}, b::NTuple{4,Float64})
    a11, a12, a13, a14, a22, a23, a24, a33, a34, a44 = A
    # Diagonal scale, for the relative pivot test below.
    scale = max(a11, a22, a33, a44)
    (isfinite(scale) && scale > 0) || return nothing
    tol = 1e-12 * scale

    # Cholesky: A = LLᵀ. Each pivot is tested before its square root, so a rank-deficient system
    # exits here rather than producing a NaN that propagates into the fit.
    a11 > tol || return nothing
    l11 = sqrt(a11)
    l21 = a12 / l11
    l31 = a13 / l11
    l41 = a14 / l11

    d22 = a22 - l21 * l21
    d22 > tol || return nothing
    l22 = sqrt(d22)
    l32 = (a23 - l31 * l21) / l22
    l42 = (a24 - l41 * l21) / l22

    d33 = a33 - l31 * l31 - l32 * l32
    d33 > tol || return nothing
    l33 = sqrt(d33)
    l43 = (a34 - l41 * l31 - l42 * l32) / l33

    d44 = a44 - l41 * l41 - l42 * l42 - l43 * l43
    d44 > tol || return nothing
    l44 = sqrt(d44)

    # Forward substitution L y = b, then back substitution Lᵀ x = y.
    y1 = b[1] / l11
    y2 = (b[2] - l21 * y1) / l22
    y3 = (b[3] - l31 * y1 - l32 * y2) / l33
    y4 = (b[4] - l41 * y1 - l42 * y2 - l43 * y3) / l44

    x4 = y4 / l44
    x3 = (y3 - l43 * x4) / l33
    x2 = (y2 - l32 * x3 - l42 * x4) / l22
    x1 = (y1 - l21 * x2 - l31 * x3 - l41 * x4) / l11

    all(isfinite, (x1, x2, x3, x4)) || return nothing
    return (x1, x2, x3, x4)
end

"""
    _r2(m, sy, syy, explained)

Coefficient of determination from accumulated sums: `1 - SS_res / SS_tot`, where
`SS_res = SS_tot - SS_reg` and `SS_reg = β̂ᵀXᵀy - (Σy)²/m` is the explained sum of squares.

`NaN` when the response has no variance across cells — for the decoupling fit that is every
timestep with all cells at or below melting, where `E` is identically zero and there is nothing for
a fit to explain.
"""
function _r2(m::Real, sy::Real, syy::Real, explained::Real)
    m > 0 || return NaN
    ss_tot = syy - sy * sy / m
    ss_tot > 0 || return NaN
    ss_reg = explained - sy * sy / m
    return clamp(ss_reg / ss_tot, 0.0, 1.0)
end

# A fitted series made safe to apply: gaps filled, out-of-domain values clamped, and both counted so
# the interval metadata can report how much of the series was measured rather than substituted.
#
# One helper for both series because the two differ only in whether the lower bound is closed. `k`'s
# domain `(0, 1]` is open at zero, so clamping to the bound itself would leave a value the validator
# still rejects — `nextfloat` puts it just inside. Keeping that subtlety in one place is the point:
# it has to be re-derived at every new call site otherwise, and the counts have to agree with the
# clamp about which values are out of domain.
#
# Returns a plain `Vector`, not whatever the caller passed: these arrive as `DimStack` layers, and
# `collect` on one preserves the `Ti` axis rather than dropping to a `Vector`.
function _make_applicable(raw, fill::Float64, (lo, hi)::Tuple{Float64,Float64};
                          open_lower::Bool, clamp_to_valid_domain::Bool)
    applied = Vector{Float64}(undef, length(raw))
    n_filled = 0
    @inbounds for (i, x) in enumerate(raw)
        if isfinite(x)
            applied[i] = x
        else
            applied[i] = fill
            n_filled += 1
        end
    end

    clamp_to_valid_domain || return (applied, n_filled, 0)

    below = open_lower ? <=(lo) : <(lo)
    lo_clamp = open_lower ? nextfloat(lo) : lo
    n_clamped = 0
    @inbounds for (i, x) in enumerate(applied)
        if below(x) || x > hi
            applied[i] = clamp(x, lo_clamp, hi)
            n_clamped += 1
        end
    end
    return (applied, n_filled, n_clamped)
end

# One cell's forcing at an interval center: lapse first, then decouple.
#
# That order is not incidental. `k` multiplies an *ambient* temperature at the glacier's elevation,
# so the elevation adjustment has to run first — `climate_adjust_for_glacier` says so explicitly,
# and `forcing_at_elevation` (the path the sweep itself uses) does the same. The two do not commute,
# since the decoupling damps only the excess above a fixed melting point and lapsing changes which
# timesteps are above it.
function _cell_forcing_at_interval(fd, delta_elevation, decoupling_factor, lapse_rate)
    lapsed = climate_adjust_for_elevation(fd, delta_elevation; lapse_rate)
    return _decouple_per_timestep(lapsed, decoupling_factor)
end

# Per-timestep decoupling. `climate_adjust_for_glacier` takes a scalar `k`, but the derived factor
# varies with the forcing's own time resolution, so apply the same increment directly:
#
#     T′  = T + (k - 1) * max(T - T_ref, 0)
#     LW′ = longwave re-derived at T′, vapor pressure unchanged
#
# written as an increment so `k = 1` is bit-exact — the property that makes a sub-freezing or
# unmeasurable timestep a true no-op. Skips the whole pass when every `k` is 1.
function _decouple_per_timestep(stack, decoupling_factor)
    all(==(1.0), decoupling_factor) && return stack

    T = stack[:temperature_air]
    e = stack[:vapor_pressure]
    LW = stack[:longwave_downward]

    T′ = similar(T)
    @inbounds for t in eachindex(T)
        k = decoupling_factor[t]
        T′[t] = T[t] + (k - 1) * max(T[t] - _DECOUPLING_REFERENCE_TEMPERATURE, 0.0)
    end
    LW′ = _adjust_longwave(LW, e, T, e, T′)

    meta = merge(copy(DimensionalData.metadata(stack)), Dict(
        "glacier_decoupling_factor_mean" => Statistics.mean(decoupling_factor),
        "temperature_air_mean" => Statistics.mean(T′)))
    # The same rebuild the upstream `climate_adjust_*` functions use, so this path cannot drift from
    # them on layer metadata or on the `rebuild` convention — and it throws on a layer name that is
    # not in the stack rather than silently adding one.
    return _rebuild_forcing(stack, meta;
                                                temperature_air = T′, longwave_downward = LW′)
end
