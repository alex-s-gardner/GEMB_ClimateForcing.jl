# Downscaling parameters: what a region's reanalysis forcing says about its own decoupling and
# lapse rate.
#
# Two fits, each a *relationship across* the grid cells of a region rather than a property of any
# one of them:
#
#   1. the local on-glacier temperature decoupling factor `k`, from how the warm excess varies with
#      the reanalysis glacier mask fraction at matched elevation ([`derive_decoupling_factor`](@ref));
#   2. the local ambient (ice-free) temperature lapse rate, from how temperature brought to the
#      ice-free state varies with elevation ([`derive_lapse_rate`](@ref)).
#
# `_cell_forcing_at_interval` applies both to one donor cell: ambient, lapse, then decouple.
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
# and the lapse rate out of temperatures brought to a common ice-free state. That is why they live
# here rather than with whatever assembles a region — the arithmetic is entirely about this
# package's own adjustment conventions.
#
# ## The `acc` argument
#
# Both fits take an accumulator of per-timestep cross-cell regression sums rather than a stack of
# forcing, so a region of hundreds of cells over decades of hourly steps never has to be held in
# memory at once. Its producer lives in the caller (GEMB_GlacierSims builds one by streaming a
# glacier elevation-class table); `_downscaling_sums` and `_accumulate_downscaling_cell!` below are
# the reference producer, and the contract is:
#
#   acc.n              — number of timesteps
#   acc.time           — the `Ti` axis values, `length(acc.time) == acc.n`
#   acc.mean_elevation — area-weighted mean of the contributing cells' reanalysis elevations (m),
#                        the reference `k` is reported at
#   acc.sums           — the sums below, or `nothing` if no cell had usable forcing
#
# `acc.sums` holds one length-`n` vector per field. With the design vector `x = [1, z, glm, glm*z]`
# per cell — `z` its reanalysis elevation (m), `glm` its glacier-mask fraction — `w = glm*z`, and
# the warm excess `E = max(T - T_ref, 0)`:
#
#   s1, z, zz, T, zT         — over every contributing cell: count (Int), Σz, Σz², ΣT, ΣzT
#   s1_warm                  — count of cells above melting (`E > 0`) at that timestep (Int)
#   z_warm, g_warm, w_warm   — Σz, Σglm, Σw over those cells
#   zz_warm, zg_warm, zw_warm, gg_warm, gw_warm, ww_warm
#                            — the upper triangle of XᵀX over those cells
#   E, zE, gE, wE            — Xᵀ against `E` (cells at or below melting add zero)
#   EE                       — ΣE², for the decoupling fit's `r2`
#
# The `_warm` sums restrict the decoupling fit to cells above melting, where the warm excess is the
# product the fit models; see `derive_decoupling_factor`.
#
# Duck-typed on purpose: a NamedTuple of `Vector{Float64}` (plus `Vector{Int}` for the counts) is
# all either fit reads, and a producer may carry further fields of its own.

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

# The closed range a donor cell's `k` must lie in before `_ambient_forcing` divides by it. The
# conversion multiplies the donor's observed warm excess, and with it the donor's residual from the
# fit, by `1 / (1 - glm·(1 - k)) <= 1 / k`, so the lower bound caps that amplification at 4. It sits
# well below any credible damping: the published on-glacier `k` of Shaw et al. (2025, the vendored
# table of 186,792 glaciers) never falls below 0.36, and a regression `k` under 0.25 at a donor is
# fit noise near the ambient-excess threshold rather than a measurement. A fitted `k` outside this
# range is excluded where the series is resolved, not clamped into it.
const _DONOR_DECOUPLING_FACTOR_LIMITS = (0.25, 1.0)

# The ten forcing layers `climate_forcing` returns. Named here rather than derived from a stack so
# a caller accumulating layer-wise can allocate before it has seen one. The last three are not
# area-averaged like the rest: an average of `solar_zenith_angle` weights `cos θ` by shortwave
# (see `shortwave_forcing_layers`).
const _FORCING_VARIABLES = (:temperature_air, :pressure_air, :vapor_pressure, :wind_speed,
                            :precipitation, :shortwave_downward, :longwave_downward,
                            :solar_zenith_angle, :shortwave_downward_diffuse,
                            :cloud_optical_thickness)

# Zeroed cross-cell sums for `n` timesteps, in the layout the `acc` contract above describes.
function _downscaling_sums(n::Integer)
    counts() = zeros(Int, n)
    v() = zeros(n)
    return (; s1 = counts(), z = v(), zz = v(), T = v(), zT = v(),
            s1_warm = counts(), z_warm = v(), g_warm = v(), w_warm = v(),
            zz_warm = v(), zg_warm = v(), zw_warm = v(), gg_warm = v(), gw_warm = v(),
            ww_warm = v(),
            E = v(), zE = v(), gE = v(), wE = v(), EE = v())
end

# Add one cell — temperature series `T` (K), reanalysis elevation `z` (m), glacier-mask fraction
# `glm` — into the sums `S`. Every vector of `S` must share `T`'s indices.
function _accumulate_downscaling_cell!(S, T, z::Real, glm::Real)
    z, glm = Float64(z), Float64(glm)
    gz = glm * z
    for t in eachindex(T, S.s1)
        Tt = T[t]
        S.s1[t] += 1
        S.z[t] += z;  S.zz[t] += z * z
        S.T[t] += Tt; S.zT[t] += z * Tt
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
    return S
end

"""
    derive_decoupling_factor(acc; min_cells = $(_MIN_CELLS_DEFAULT)) -> DimStack

Fit the on-glacier temperature decoupling factor `k` per timestep from cross-cell sums, as
accumulated over a region's grid cells.

`climate_adjust_for_glacier` damps the warm excess above melting, `T′ = T + (k - 1)·max(T - T_ref, 0)`,
so `E′ = k·E`. The estimator inverts that using the one place the reanalysis already varies its own
degree of decoupling: `glm`, the fraction of each cell its land-surface scheme treats as ice.

**Which end of `glm` is decoupled matters, and it is the high end.** A cell with `glm = 1` is one
ERA5-Land already runs as ice, so its near-surface temperature is *already* damped. A cell with
`glm = 0` is run as ice-free ground and its temperature is ambient. Mixing the two, a cell's
observed excess is

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

**Only cells above melting are fitted.** The model is a product of the damping factor and the
ambient excess, and the excess is truncated at zero: a cell at or below melting reports `E = 0`
whatever its ambient temperature, so it lies off the linear surface rather than on it. Where the
melting level falls inside the region, including those cells bends the fitted surface toward zero
at the cold end and biases `k` low (0.74 against a true 0.80 on a synthetic region whose top
cells are below melting). Under the model a cell is above melting exactly when its ambient
temperature is, since the damping factor is positive, so restricting the fit to cells with
`E > 0` loses nothing and the four-term surface is exact on what remains.

Returns a `DimStack` on the accumulated `Ti` axis:

| layer | meaning |
|:--|:--|
| `decoupling_factor` | `k` at the reference elevation, raw. `NaN` where the region supported no fit |
| `n_cells` | grid cells contributing to that timestep |
| `n_cells_above_melting` | of those, the cells above melting — the ones the fit uses |
| `r2` | coefficient of determination of the four-term fit, over the cells above melting |
| `ambient_excess` | the fitted ice-free excess `α + β·z̄` (K) — what the fit's precision hinges on |
| `coef_alpha`, `coef_beta`, `coef_gamma`, `coef_delta` | all four coefficients, so `k` is re-evaluable at any elevation — see [`decoupling_factor_at_elevation`](@ref) |

`metadata` carries `reference_elevation` (`z̄`), `min_cells`, `n_fitted`, and `n_timesteps`.

**`NaN` means "the forcing carried no information here", and nothing more.** A timestep reports `NaN`
when it has fewer than `max(min_cells, 4)` cells above melting (four parameters need four cells), no
spread in `glm` or elevation among them, collinear regressors, or an ambient excess below `$(_MIN_AMBIENT_EXCESS)` K. That last
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

`ambient_excess` is the honest measure of how well determined a given timestep's `k` is: the closer
it sits to `$(_MIN_AMBIENT_EXCESS)` K, the more of the ratio is noise. The coefficients are what
[`derive_lapse_rate`](@ref) uses to bring each cell to the ice-free state.

## Screening `k` downstream, and why it also screens the lapse rate

`k` is the ratio `[(α+γ) + (β+δ)·z̄] / (α + β·z̄)`, whose **denominator is exactly
`ambient_excess`**. `$(_MIN_AMBIENT_EXCESS)` K is only enough to keep the ratio finite, not enough to
make it meaningful, so an out-of-domain `k` is the expected outcome near the threshold rather than a
sign of a broken fit — `r2` is a median 0.97 for those timesteps, statistically indistinguishable from
the well-behaved ones. `r2` cannot be used to screen them; `ambient_excess` can.

Two years of ERA5-Land at Wrangell-St Elias (38 cells, 3427 timesteps with a measurable `k`), with
every cell in the fit, including those at or below melting, so the figures indicate scale rather
than the current fit's exact behavior:

| `_MIN_AMBIENT_EXCESS` | timesteps kept | with `k` outside `(0,1]` | max `k` |
|--:|--:|--:|--:|
| 0.5 (current) | 3427 (100%) | 788 (23.0%) | 7.35 |
| 1.0 | 3174 (93%) | 634 (20.0%) | 4.33 |
| 2.0 | 2770 (81%) | 456 (16.5%) | 2.63 |
| 3.0 | 2346 (68%) | 287 (12.2%) | 2.00 |

Raising the threshold trades coverage for tail control and never eliminates the tail, which is why it
stays at $(_MIN_AMBIENT_EXCESS) K here — the measurement guard — and the rest is a post-processing
decision. **Screen on `k ∈ (0, 1]` directly**, which is what
`climate_adjust_for_glacier` requires anyway. [`derive_lapse_rate`](@ref) applies the same screen:
at those timesteps it fits the observed temperatures without the ice-free correction.

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
gates whole regions, and the cost of raising it is measured mostly in regions dropped. The count above
melting, which is what the gate is applied to, also varies with time: a timestep whose melting level
leaves fewer than `min_cells` cells above it reports `NaN`. The subsampling and coverage figures below
were measured with every cell counted.

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
    n_cells_above_melting = zeros(Int, n)
    # The fitted ice-free excess at `z̄`, per timestep. `NaN` where nothing was fitted.
    ambient_excess = fill(NaN, n)
    # All four fitted coefficients per timestep, so `k` can be re-evaluated at an elevation other
    # than `z̄` — which is what every elevation interval needs, since `z̄` is a reanalysis surface
    # and the intervals are the glacier. See `decoupling_factor_at_elevation`.
    cα, cβ, cγ, cδ = fill(NaN, n), fill(NaN, n), fill(NaN, n), fill(NaN, n)

    S = acc.sums
    if S !== nothing
        for t in 1:n
            n_cells[t] = S.s1[t]
            m = S.s1_warm[t]
            n_cells_above_melting[t] = m
            # Only cells above melting are fitted, and four parameters need four of them before the
            # system is even determined.
            m >= max(min_cells, 4) || continue

            # XᵀX for x = [1, z, glm, glm*z] over the warm cells, upper triangle, and Xᵀy for y = E.
            A = (Float64(m), S.z_warm[t],  S.g_warm[t],  S.w_warm[t],
                             S.zz_warm[t], S.zg_warm[t], S.zw_warm[t],
                                           S.gg_warm[t], S.gw_warm[t],
                                                         S.ww_warm[t])
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
            r2[t] = _r2(m, S.E[t], S.EE[t],
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
                     n_cells_above_melting = DimArray(n_cells_above_melting, ti),
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
    for t in 1:n
        isfinite(α[t]) || continue
        a, b = α[t], β[t]
        # Evaluate at `z` if the ambient excess still clears the threshold there, and otherwise at
        # the edge of the half-line where it does: `E_ambient = a + b·z` is linear, so that edge is
        # `z* = (threshold - a) / b`, and a fitted timestep has `z̄` inside it. At the edge the
        # ambient excess is the threshold itself by construction, so it is used directly rather
        # than re-evaluated — rounding in `a + b·z*` could otherwise land a hair below the
        # threshold and drop a held timestep.
        if a + b * z >= _MIN_AMBIENT_EXCESS
            ambient = a + b * z
            z_eval = z
        else
            b == 0 && continue                     # constant and already below: nothing to hold onto
            ambient = _MIN_AMBIENT_EXCESS
            z_eval = (_MIN_AMBIENT_EXCESS - a) / b
            n_held += 1
        end
        k[t] = (ambient + (γ[t] + δ[t] * z_eval)) / ambient
    end

    return DimArray(k, ti;
                    metadata = Dict{String,Any}("elevation" => z, "in_range" => true,
                                                "n_held" => n_held,
                                                "n_fitted" => count(isfinite, k)))
end

"""
    derive_lapse_rate(acc, decoupling; min_cells = $(_MIN_CELLS_DEFAULT)) -> DimStack

Fit the **ambient** (ice-free) near-surface temperature lapse rate (K/km, positive for cooling with
height) per timestep from cross-cell sums, after bringing every cell to the ice-free state with the
decoupling fit.

The slope wanted is the ambient one because the forcing is applied ambient-first: a donor cell's
temperature is taken to the ice-free state, lapsed to the band elevation with this slope, and only
then decoupled with the band's own `k` ([`_cell_forcing_at_interval`](@ref)). Decoupling the lapsed
temperature once is then the whole of the on-glacier correction; a slope fitted to decoupled
temperatures would apply the damping to the vertical term twice.

Each cell's temperature is brought to the ice-free state with the four-term model
[`derive_decoupling_factor`](@ref) fits, `E_i = α + β·z_i + glm_i·(γ + δ·z_i)`, whose ice-free part
is `α + β·z_i`. Removing the glacier term leaves the cell's ambient excess plus its own residual from
the fit:

    T_ambient_i = T_i - glm_i·(γ + δ·z_i)

for a cell above melting, and `T_ambient_i = T_i` for a cell at or below it, where the reanalysis
applies no damping. An ice-free cell (`glm = 0`) is unchanged and a fully glaciated one gets the
whole correction. The correction is linear in the warm-cell sums already accumulated, so the fit
needs no second pass over the forcing. [`_ambient_forcing`](@ref) brings a donor cell to the same
ice-free state through the ratio form of the model, `E_ambient = E_obs / (1 - glm·(1 - k(z)))` with
`k(z) = 1 + (γ + δ·z)/(α + β·z)`. On forcing that follows the model the two are identical; they
differ only in how a cell's residual from the fit is carried (unscaled here, scaled by the ratio
there).

`decoupling` is the `DimStack` [`derive_decoupling_factor`](@ref) returns, read **raw**. A timestep
whose `k` is `NaN` *or* outside `(0, 1]` is fitted on the observed temperatures without correction.
For a `NaN` that is right in practice: `k` is unmeasurable where the region carries too little warm
excess for the damping to matter to the slope. An out-of-domain `k` comes from a ratio of two
near-zero fitted excesses, and its coefficients would add a correction proportional to `glm`. `glm`
rises with elevation in most regions, so that correction is a smooth, near-collinear false gradient
rather than scatter, and no slope estimator can reject it. Nothing substituted is written back into
either output. To know which timesteps were corrected, test `0 < decoupling_factor <= 1` on the
companion stack.

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

Two parameters rather than four makes this the less demanding of the two fits, so `min_cells` is set
by the decoupling fit's needs, not this one's. See [`derive_decoupling_factor`](@ref).

**Nothing is clamped.** A fitted slope of 200 K/km is real evidence that this region does not
constrain a lapse rate at that timestep, and clamping it into $(_LAPSE_RATE_LIMITS) K/km would hide
exactly that. The clamp lives in `GEMB_GlacierSims.elevation_interval_forcing`, the one place a
slope is applied and so the one place it has to be acceptable to `climate_adjust_for_elevation`.

Summarise with a **median**, not a mean, and expect real spread: a per-timestep slope over a handful
of cells is broad, and much of that breadth is physical — nocturnal valley inversions and steep
midday cooling are both in the record. Group over the `Ti` axis at whatever resolution the decision
needs:

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
    cγ, cδ = decoupling.coef_gamma, decoupling.coef_delta

    S = acc.sums
    if S !== nothing
        for t in 1:n
            n_cells[t] = S.s1[t]
            S.s1[t] >= min_cells || continue

            m = Float64(S.s1[t])

            # Elevation spread. Zero means every cell in the region sits at the same reanalysis
            # elevation, so no slope is identifiable however many cells there are.
            Szz = S.zz[t] - S.z[t] * S.z[t] / m
            Szz > 0 || continue
            elevation_spread[t] = sqrt(Szz / m)

            # Bring every cell to the ice-free state: `T_i - glm_i·(γ + δ·z_i)` for a cell above
            # melting, unchanged otherwise. Summed, that subtracts `γ·Σglm + δ·Σglm·z` from ΣT and
            # `γ·Σz·glm + δ·Σz·glm·z` from ΣzT, the sums running over the warm cells. Uncorrected
            # where `k` is unmeasured or out of domain — see the docstring.
            γ, δ = _corrects_to_ambient(k[t], cγ[t], cδ[t]) ? (cγ[t], cδ[t]) : (0.0, 0.0)
            sT  = S.T[t]  - (γ * S.g_warm[t]  + δ * S.w_warm[t])
            szT = S.zT[t] - (γ * S.zg_warm[t] + δ * S.zw_warm[t])

            # Simple linear regression of the ambient temperature on elevation.
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

# Whether a timestep's decoupling fit is used to bring the cells of the lapse fit to the ice-free
# state: `k` at the reference elevation measured and inside `(0, 1]`, and the glacier coefficients
# finite. Used by `derive_lapse_rate` and by GEMB_GlacierSims' lapse-rate uncertainty, so the slope
# and its standard error describe the same regression.
_corrects_to_ambient(k, γ, δ) =
    isfinite(k) && _DECOUPLING_FACTOR_LIMITS[1] < k <= _DECOUPLING_FACTOR_LIMITS[2] &&
    isfinite(γ) && isfinite(δ)

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

"""
    _cell_forcing_at_interval(fd, delta_elevation, decoupling_factor, lapse_rate,
                              glacier_fraction, donor_decoupling_factor) -> DimStack

One donor cell's forcing `fd`, delivered on glacier at an elevation `delta_elevation` (m) above the
cell's reanalysis surface:

1. **to ambient** — undo the damping the reanalysis applied to the cell's own warm excess
   ([`_ambient_forcing`](@ref)), with `glacier_fraction` the donor cell's `glm` and
   `donor_decoupling_factor` the per-timestep `k` at the **donor's** elevation;
2. **lapse** — [`climate_adjust_for_elevation`](@ref) with the ambient `lapse_rate` (K/km, scalar,
   monthly or per-timestep) that [`derive_lapse_rate`](@ref) fits;
3. **decouple** — the band is glacier, so the band's full per-timestep `decoupling_factor` at the
   **band's** elevation is applied to the lapsed ambient temperature
   ([`_decouple_per_timestep`](@ref)).

Step 1 and the two fits share one model, `E_obs = [1 - glm·(1 - k(z))]·E_ambient(z)`, so a region
whose forcing follows it gives the same band temperature `T_ref + k(z_band)·E_ambient(z_band)` from
every donor. The band's `k` is applied once, in step 3: step 1 removes the reanalysis's own damping
of the donor, and the slope of step 2 is the ambient one. The order is not incidental: the
decoupling damps only the excess above a fixed melting point, and lapsing changes which timesteps
are above it.

`decoupling_factor` must be finite and in `(0, 1]`, `donor_decoupling_factor` finite and in
`$(_DONOR_DECOUPLING_FACTOR_LIMITS)` (see [`_ambient_forcing`](@ref)), and `glacier_fraction` in
`[0, 1]`; anything else throws an `ArgumentError`. Both factors come from
[`decoupling_factor_at_elevation`](@ref), at the band's and the donor's elevation, with their gaps
resolved by the caller.
"""
function _cell_forcing_at_interval(fd, delta_elevation, decoupling_factor, lapse_rate,
                                   glacier_fraction, donor_decoupling_factor)
    _check_decoupling_series(decoupling_factor, fd, "decoupling_factor")
    ambient = _ambient_forcing(fd, glacier_fraction, donor_decoupling_factor)
    lapsed = climate_adjust_for_elevation(ambient, delta_elevation; lapse_rate)
    return _decouple_per_timestep(lapsed, decoupling_factor)
end

function _check_decoupling_series(k, fd, name)
    _check_series_length(k, fd, name)
    all(x -> isfinite(x) && 0 < x <= 1, k) || throw(ArgumentError(
        "$(name) must be finite and in (0, 1] at every timestep; fill and clamp a fitted series " *
        "before applying it"))
    return nothing
end

function _check_series_length(x, fd, name)
    length(x) == length(dims(fd, Ti)) || throw(DimensionMismatch(
        "$(name) has $(length(x)) values but the forcing has $(length(dims(fd, Ti))) timesteps"))
    return nothing
end

"""
    _ambient_forcing(fd, glacier_fraction, decoupling_factor) -> DimStack

A reanalysis cell's forcing brought to the ice-free (ambient) state by inverting the damping the
decoupling fit measures. [`derive_decoupling_factor`](@ref) models a cell with glacier-mask fraction
`glm` that is above melting as reporting `E_obs = [1 - glm·(1 - k)]·E_ambient`, with `k` the
decoupling factor at the cell's own elevation, so

    E_ambient = E_obs / D,    D = 1 - glm·(1 - k) = k + (1 - k)·(1 - glm),    E_obs = T - T_ref > 0

`decoupling_factor` is that `k` per timestep (from [`decoupling_factor_at_elevation`](@ref) at the
cell's elevation, gaps resolved). A step at or below melting (`E_obs <= 0`) is unchanged: the fit
uses only cells above melting, and under the model a cell is above melting exactly when its ambient
temperature is, since `D > 0`. Vapor pressure, longwave and every other layer are unchanged (warming
cannot bring air past saturation).

The division is bounded. For any `glm` in `[0, 1]`, `D` lies between `k` and 1, so with `k` in
`$(_DONOR_DECOUPLING_FACTOR_LIMITS)` the conversion multiplies the observed excess, and the cell's
residual from the fit with it, by at most `1/$(_DONOR_DECOUPLING_FACTOR_LIMITS[1])`. A `k` near zero
would make the inversion unbounded (a fully glaciated donor whose fit gives `k = 1e-7` and that
observes a 2 K residual would be warmed by 2e7 K); such a `k` is refused here, and excluded where a
fitted series is resolved. `D` is formed as `k + (1 - k)·(1 - glm)`, a sum of non-negative terms, so
`D >= k` also holds in floating point, and `D` at or above the lower limit is asserted at every
step.

Written as the increment `T + E_obs·glm·(1 - k)/D`, so `glm = 0` or `k = 1` returns the input bit
for bit; when that holds at every step the input stack itself is returned. `glacier_fraction`
outside `[0, 1]`, or a `decoupling_factor` that is non-finite, outside
`$(_DONOR_DECOUPLING_FACTOR_LIMITS)` or of the wrong length, throws.
"""
function _ambient_forcing(fd, glacier_fraction::Real, decoupling_factor)
    glm = Float64(glacier_fraction)
    0 <= glm <= 1 || throw(ArgumentError("glacier_fraction must be in [0, 1] (got $(glacier_fraction))"))
    _check_series_length(decoupling_factor, fd, "donor_decoupling_factor")
    k_min, k_max = _DONOR_DECOUPLING_FACTOR_LIMITS
    all(k -> isfinite(k) && k_min <= k <= k_max, decoupling_factor) || throw(ArgumentError(
        "donor_decoupling_factor must be finite and in [$(k_min), $(k_max)] at every timestep; " *
        "resolve a fitted series before applying it"))
    (glm == 0 || all(==(1.0), decoupling_factor)) && return fd

    T = fd[:temperature_air]
    T′ = similar(T)
    for t in eachindex(T, decoupling_factor)
        E = T[t] - _DECOUPLING_REFERENCE_TEMPERATURE
        if E > 0
            k = decoupling_factor[t]
            D = k + (1 - k) * (1 - glm)
            D >= k_min || error("ambient conversion denominator $(D) at step $(t) is below " *
                                "$(k_min) for k = $(k), glm = $(glm); this is a bug")
            T′[t] = T[t] + E * (glm * (1 - k)) / D
        else
            T′[t] = T[t]
        end
    end
    meta = merge(copy(DimensionalData.metadata(fd)),
                 Dict("temperature_air_mean" => Statistics.mean(T′)))
    return _rebuild_forcing(fd, meta; temperature_air = T′)
end

# Per-timestep decoupling. `climate_adjust_for_glacier` takes a scalar `k`, but the derived factor
# varies with the forcing's own time resolution, so apply the same correction directly:
#
#     T′ = T + (k - 1) * max(T - T_ref, 0)
#     e′ = min(e, eₛ(T′)) where T′ < T, e elsewhere
#
# with longwave and every other layer unchanged — the conventions of `climate_adjust_for_glacier`.
# Written as an increment so `k = 1` is bit-exact, the property that makes a sub-freezing or
# unmeasurable timestep a true no-op. Skips the whole pass when every `k` is 1.
function _decouple_per_timestep(stack, decoupling_factor)
    all(==(1.0), decoupling_factor) && return stack

    T = stack[:temperature_air]
    e = stack[:vapor_pressure]

    T′ = similar(T)
    for t in eachindex(T, decoupling_factor)
        T′[t] = T[t] + (decoupling_factor[t] - 1) * max(T[t] - _DECOUPLING_REFERENCE_TEMPERATURE, 0.0)
    end
    e′ = _cap_vapor_pressure.(e, T, T′)

    meta = merge(copy(DimensionalData.metadata(stack)), Dict(
        "glacier_decoupling_factor_mean" => Statistics.mean(decoupling_factor),
        "temperature_air_mean" => Statistics.mean(T′)))
    # The same rebuild the upstream `climate_adjust_*` functions use, so this path cannot drift from
    # them on layer metadata or on the `rebuild` convention — and it throws on a layer name that is
    # not in the stack rather than silently adding one.
    return _rebuild_forcing(stack, meta; temperature_air = T′, vapor_pressure = e′)
end
