# Correcting off-glacier air temperature to on-glacier air temperature

Research notes for a possible `climate_adjust_for_glacier` workflow in
GEMB_ClimateForcing.jl. Local copies of every paper cited here are in
`docs/references/`.

## The problem

A melting glacier surface is pinned near 0 °C. Warm air advected over it is cooled from
below, which (i) creates a stable near-surface layer that suppresses turbulent mixing and
(ii) drives a shallow density-driven katabatic ("glacier") wind. The result is a **glacier
boundary layer (GBL)** whose 2 m air temperature is *decoupled* from the ambient
free-atmosphere temperature at the same elevation.

Reanalysis products (ERA5-Land included) do not resolve this. Their grid cells are far
larger than a valley glacier and their land-surface schemes carry no glacier wind, so
their 2 m temperature is effectively an **ambient (off-glacier) temperature**, `T_amb`.
Feeding it straight into a surface energy balance model is a known bias source:

- Extrapolating an off-glacier record to on-glacier elevations with a constant lapse rate
  **overestimates melt**. Greuell & Böhm (1998) put the overestimate of the mass-balance
  change per +1 °C at **22 %**. Troxler et al. (2020) found up to **32 %** of the seasonal
  melt total at one McCall Glacier site.
- Kaushik et al. (2025, EGU abstract) report lapse-rate-based models overestimating point
  mass balance by **up to 92 %** on Chhota Shigri, cut to **8 %** using SM10 below.
- Getting temperature right matters most for the turbulent fluxes and for longwave
  parameterisations — exactly what GEMB computes.

Critically, the correction is **not a constant offset and not a modified lapse rate**. The
cooling switches on only above a threshold ambient temperature (when the katabatic layer
forms) and grows with distance travelled over ice. Below the threshold, on-glacier and
ambient temperatures track each other ~1:1.

## Four approaches, in increasing order of usefulness here

### 1. Linear lapse rate (the baseline — inadequate for warm conditions)

Well documented and consistently insufficient. From the distributed-logger studies:

| Condition | Observed on-glacier lapse rate | R² vs elevation |
|---|---|---|
| Cold ambient (P10) | −4.7 to −10.9 °C/km (near or steeper than ELR) | 0.83–0.97 |
| Warm ambient (P90) | −1.4 to −4.5 °C/km (shallow) | 0.05–0.28 (often not significant) |

(Shaw et al. 2017 Table 3; Troxler et al. 2020 Fig. 6; Shaw et al. 2021.)

So a lapse rate works acceptably under **cold** conditions and fails under **warm** ones —
precisely the hours that dominate melt. Note the elevation dependency doesn't just weaken,
it can invert: temperature sometimes *increases* toward the terminus.

### 2. Shea & Moore (2010) — "SM10", statistical, DEM-only

*The original paper is now held locally (`docs/references/Shea2010_JGR_distributed_Ta_vapor_pressure.pdf`),
so everything in this section is first-hand rather than reconstructed.*

The most practical per-glacier scheme. Piecewise-linear in ambient temperature, with a
katabatic-onset threshold `T*`:

```
T_gla = T1 + k_hi·(T_amb − T*)     for T_amb ≥ T*    (KBL developed, damped sensitivity)
T_gla = T1 − k_lo·(T*  − T_amb)    for T_amb <  T*    (no KBL, near 1:1)
T1    = k_lo · T*
```

> **Notation warning — the source paper is internally inconsistent.** Shea & Moore's Eq. 1
> as printed labels the *upper* branch slope `k1` and the *lower* branch `k2`, but their
> §4.1 text, Table 5 and Figure 7 all use the opposite convention: `k1` is the
> below-threshold slope (0.63–0.95, near 1:1) and `k2` the above-threshold damped slope
> (0.32–0.72). **The tables and figures are the usable convention**, and it is the one the
> whole downstream literature (Shaw 2017/2021, Carturan 2015) adopted. Above I use
> `k_lo`/`k_hi` to sidestep the collision; in the β table below, `k1` = `k_lo` and
> `k2` = `k_hi`, matching Table 5.

`k_lo`, `k_hi` are the temperature *sensitivities* (ratio of on-glacier to ambient
temperature change) below and above threshold. Both are functions of **flowpath length**
`FPL` (m) — the average flow distance to a point from an upslope summit/ridge. SM10 computed
it from LiDAR DEMs with the Quinn et al. (1991) multiple-flow-direction algorithm in SAGA
GIS (Shaw's later papers use MATLAB TopoToolbox instead):

```
k1 = β3 · exp(β4·FPL)            (SM10 Eq. 6)
k2 = β5 + β6 · exp(β7·FPL)       (SM10 Eq. 7)
```

Published coefficient sets (Shaw et al. 2021 Table 3; Shaw et al. 2017 Table 5):

| Set | β1 | β2 | β3 | β4 | β5 |
|---|---|---|---|---|---|
| SM10 original (Coast Mts, Canada) | 0.977 | −4.4×10⁻⁵ | 0.29 | 0.71 | 5.6×10⁻³ |
| Parlung (SE Tibet) | 0.894 | −2.972×10⁻⁵ | 0.349 | 0.624 | 4.4×10⁻³ |
| All glaciers, no terminus warming | 0.923 | −3.375×10⁻⁵ | 0.343 | 0.511 | 4.2×10⁻³ |
| SMopt (Tsanteleina, Italy) | 0.913 | −2.9×10⁻⁵ | 0.643 | 1.0 | −19.07×10⁻⁴ |

Note the sign convention on β5 differs between papers (`exp(−β5·DF)` in Shaw et al. 2021
Eq. 3 vs `exp(β5·DF)` in Shaw et al. 2017 Eq. 4) — the intent is always decay with
distance, so check the sign against a `k2→β3` limit at large `DF` before use. SM10's own
Figure 7 confirms the intended behaviour directly: `k1` falls roughly linearly from ~0.95 at
FPL ≈ 1000 m to ~0.63 at 10 000 m, and `k2` decays sharply from ~0.72 to an asymptote of
~0.3 reached by FPL ≈ 5000 m.

SM10 did not publish its β values in a table — they appear only as fitted curves in its
Figures 7 and 9. The per-site parameters *are* tabulated, and are the more useful anchor for
verifying an implementation (SM10 Table 5, "global" models; FPL from Table 2):

| Site | FPL (m) | T1 (°C) | T* (°C) | k1 | k2 | n | R² | RMSE (°C) |
|---|---|---|---|---|---|---|---|---|
| BM1 (Bridge) | 10075 | 4.08 | 6.08 | 0.63 | 0.32 | 1849 | 0.82 | 1.01 |
| PM1 | 2829 | 5.51 | 6.48 | 0.85 | 0.43 | 3429 | 0.90 | 1.04 |
| PM2 | 2168 | 4.76 | 5.60 | 0.89 | 0.48 | 7174 | 0.91 | 1.10 |
| PM3 | 1036 | 4.55 | 5.12 | 0.95 | 0.69 | 3854 | 0.94 | 0.97 |
| PM4 | 970 | 3.89 | 4.29 | 0.95 | 0.72 | 5075 | 0.93 | 1.25 |
| PM108 | 987 | 4.82 | 5.35 | 0.91 | 0.64 | 2355 | 0.96 | 0.86 |
| PM308 | 2076 | 4.82 | 5.12 | 0.93 | 0.62 | 2799 | 0.96 | 0.80 |
| WM1 (Weart) | 3522 | 3.70 | 4.90 | 0.83 | 0.35 | 3418 | 0.68 | 1.54 |
| WM2 | 1697 | 6.74 | 8.37 | 0.83 | 0.60 | 1410 | 0.87 | 1.17 |

Independent-site validation (Table 7) gives temperature RMSE 0.95–1.32 °C and MBE −0.53 to
+0.27 °C — i.e. predicting the coefficients from FPL costs essentially nothing versus fitting
them per site. Note WM2's `T*` = 8.37 °C is an outlier SM10 itself excluded from the `T*`
model, blaming undersampling at low temperatures.

For `T*`, SM10 originally used elevation (`T* = β1 + β2·Z`, its Eq. 5 — fitted on station
elevation, significant at the 5 % level), which does not transfer between climates. Carturan
et al. (2015) replaced it with a flowpath-length form that is now standard:

```
T* = C1·DF / (C2 + DF),   C1 = 6.61 °C,  C2 = 436.04 m
```

**Behaviour.** `k2` drops steeply over the first 2000–3000 m of flowpath, then flattens:
~0.7–0.8 for small glaciers (<1000 m flowpath), asymptoting to ~0.3 (as low as 0.2) at
>7000 m. So a small glacier is nearly coupled to ambient temperature; a long one is
strongly damped.

**Strengths.** Few parameters, needs only a DEM plus an off-glacier temperature series,
published coefficients validated on 11+ glaciers worldwide, explicitly represents the
threshold behaviour. Shaw et al. (2021) recommend it over a static lapse rate for
data-sparse glaciers, and Shaw et al. (2017) found melt-model RMSE improvements of
**28–36 %**.

**Weaknesses.** Purely statistical; no representation of warming on glacier tongues (see
below). Carturan et al. (2015) got RMSE = 0.73 °C, ME = −0.06 °C over all sites, but
performance degrades where valley or synoptic winds displace the katabatic flow.

### 3. Greuell & Böhm (1998) + Ayala et al. (2015) — "GB" / "ModGB", physical

Thermodynamic model of an air parcel descending the flowline, governed by adiabatic heating
and sensible heat exchange with the surface. Ayala et al. (2015) added the final term to
capture observed *warming* on lower tongues:

```
T_a(x) = (T0 − T_eq)·exp(−(x − x0)/L) + T_eq  +  K·(x − x0)/L
         └──────────── original GB ────────────┘  └── ModGB addition ──┘

L     = H·cos(α) / C_H          characteristic length scale (m)
T_eq  = b·L                     "equilibrium temperature"
b     = Γ_d · tan(α)            Γ_d = −0.0098 °C/m
```

with `x` = flowline distance, `T0` = temperature at the top of the flowline (`x0`, where
air enters the GBL), `α` = mean glacier slope, `C_H` ≈ 0.002 (bulk transfer coefficient),
`H` = boundary layer height (m), `K` (°C) = empirical warming factor.

`H` and `K` are tuning parameters fit to on-glacier observations. Published values:

| Glacier | H (m) | K (°C) |
|---|---|---|
| Haut Glacier d'Arolla (Ayala et al. 2015) | 4–8 | 6–8 (above the 8 °C threshold) |
| Tsanteleina (Shaw et al. 2017) | 2.8 | 6.6 |
| McCall (Troxler et al. 2020), median | 7.6 | 4.2, K/L = 1.2 °C/km |
| Pasterze (Greuell & Böhm 1998), original GB | x0 = 1440 m, L = 8340 m | — |

Applied only for warm conditions (`T0 ≥ 6 °C` in Shaw et al. 2017; `T_amb > 8 °C` in Ayala
et al. 2015); a linear lapse rate is used below that. Reported gains over a linear lapse
rate are real but modest: RMSE reduction of **0.2–0.5 °C**, better in 7 of 10 years at
McCall.

**Verdict for our purposes: not usable without on-glacier data.** `H` and `K` must be
calibrated against distributed on-glacier observations along the lower flowline, and
transferability is explicitly unresolved (Troxler et al. 2020; Shaw et al. 2021 note ModGB's
unknown parameters "can lead to high variability in `T_a` estimates on the lower glacier
ablation zone"). Shaw et al. 2021 §6.2 prefer SM10 for exactly this reason.

### 4. Shaw et al. (2025, Nature Climate Change) — global decoupling factor `k`

The newest and, for a global/regional forcing package, the most directly applicable. Built
from 169 glacier-summers (62 unique glaciers, 1994–2023), it defines a single **decoupling
factor**

```
k = (T_gla − β) / T_amb
```

(`β` = intercept of the `T_amb`–`T_gla` regression), then predicts `k` globally from a
multiple linear regression on five predictors:

```
k ≈ a0 + a1·Ta + a2·Q + a3·FF + a4·FPL + a5·ELE
```

| Parameter | Unit | Estimate [95 % CI] | rel. importance |
|---|---|---|---|
| a0 (intercept) | – | 0.9912 [0.91, 1.07] | — |
| a1 (off-glacier air temperature) | °C⁻¹ | −0.0072 [−0.010, −0.003] | 16.5 % |
| a2 (off-glacier specific humidity) | kg g⁻¹ | −37.474 [−49.8, −25.1] | 28.5 % |
| a3 (ERA5-Land synoptic wind speed) | m⁻¹ s | 0.021256 [0.005, 0.037] | 13.5 % |
| a4 (flowpath length) | m⁻¹ | 1.32×10⁻⁵ | 39.0 % |
| a5 (elevation) | m⁻¹ | 1.31×10⁻⁵ | 2.5 % |

R² = 0.60, RMSE = 0.105. `k` is clamped to [0.2, 1.0]; debris-covered bands are set to 1
with a linear reduction toward up-glacier clean ice. Mean cooling follows as
`MC = T_amb·k − T_amb`.

#### The published coefficients do not reproduce the published `k` — use the lookup table

I tested Table S2 against the authors' own `HISTORICAL_DECOUPLING_ESTIMATES_DATABASE.mat`,
which ships all five predictors (`TA`, `Q`, `FF`, `LEN`, `ELE`) alongside the resulting `K`
for 186 792 RGI v6 glaciers. Applying the table verbatim does **not** recover `K`:

| Variant of Table S2 | mean k | corr with published K | RMSE |
|---|---|---|---|
| Verbatim (`+a4·LEN +a5·ELE`) | 0.903 | 0.39 | 0.106 |
| **Signs on a4, a5 flipped** | 0.786 | **0.66** | **0.076** |
| Published K (truth) | 0.833 | — | — |

So **a4 and a5 are sign-typos in Table S2**: `k` *decreases* with glacier length and
elevation, it does not increase. A direct OLS refit on the authors' own predictors gives:

```
k = 0.9347 − 0.005960·TA − 20.033·Q + 0.021731·FF − 5.779e-6·LEN − 1.1375e-5·ELE
    (TA °C, Q kg/kg, FF m/s ×2.5 as shipped, LEN m, ELE m;  R² = 0.536, RMSE = 0.046)
```
Even this only reaches R² = 0.54 against their `K`, so the published fields were not generated
by a plain 5-term linear model on these columns alone — there is additional processing (the
debris adjustment, and per-elevation-band aggregation, are the likely culprits). The regression
form is therefore **not** reliably reproducible from the paper.

Two further discrepancies worth knowing:
- **`a2` is documented in "kg g⁻¹" but the shipped `Q` is in kg/kg** (mean 0.0047), and the
  −37.474 coefficient only makes physical sense against kg/kg. The README's claim that `Q` is
  in g kg⁻¹ contradicts the actual file contents. Do not rescale.
- **`FF` in the shipped file has already been multiplied by 2.5** ("to best match off-glacier
  wind speed observations close to glaciers", per the README). Any use of the regression with
  raw ERA5-Land wind must apply that factor first.
- **The paper's `LEN` is RGI glacier length, not the flowpath length `FPL`** used in the
  observations database and in SM10. These are different quantities; the observations file
  carries `AWS_FPL` separately.

**Conclusion: use the per-glacier lookup table, not the regression.** The `.mat` files are
HDF5 (MATLAB v7.3), so `HDF5.jl` reads them with no MATLAB dependency; `RGI_ID` is an array of
object references that dereference to strings like `"RGI60-01.00001"`. Regional mean `k`,
computed from the shipped table:

| RGI region | n | mean k | | RGI region | n | mean k |
|---|---|---|---|---|---|---|
| 01 Alaska | 27057 | 0.824 | | 11 Central Europe | 3925 | 0.829 |
| 02 W Canada/US | 18852 | 0.801 | | 12 Caucasus | 1549 | 0.833 |
| 03 Arctic Canada N | 3706 | 0.951 | | 13 Central Asia | 52735 | 0.819 |
| 04 Arctic Canada S | 6432 | 0.902 | | 14 South Asia W | 27693 | 0.835 |
| 06 Iceland | 389 | 0.843 | | 15 South Asia E | 12830 | 0.768 |
| 07 Svalbard | 942 | 0.891 | | 16 Low Latitudes | 2938 | 0.884 |
| 08 Scandinavia | 3412 | 0.785 | | 17 Southern Andes | 15248 | 0.925 |
| 09 Russian Arctic | 514 | 0.893 | | 18 New Zealand | 3537 | 0.841 |
| 10 North Asia | 5033 | 0.828 | | | | |

Region 05 (Greenland periphery) and 19 (Antarctic) are **absent** — the dataset explicitly
excludes ice sheets, ice caps and the large Antarctic/Greenland periphery glaciers. That is a
hard coverage limit for a GEMB forcing package, which is often pointed at exactly those.
Also note `COOL` in the file is *not* simply `TA·K − TA` (max abs difference 3.87 °C), despite
the README's definition — it is averaged over elevation bands before being written out.

The observations database (`GLACIER_DECOUPLING_OBSERVATIONS_DATABASE.mat`, 186 glacier-years,
415 AWS points with `AWS_K`, `AWS_R2`, `AWS_FPL`) is the real training set and is a genuinely
useful validation target. Fitting the 5-predictor model directly to it reproduces neither the
published coefficients nor R² = 0.60 (I get R² = 0.11–0.24 per AWS point, 0.14 per
glacier-year), so the paper's reported skill must come from a different aggregation than any
of the obvious ones. Notably, in the observations `k` correlates *negatively* with `AWS_FPL`
(r = −0.22) — consistent with SM10 and with the sign-flipped a4, and confirming Table S2's
positive a4 is wrong.

Headline numbers: global mean `k` ≈ **0.83** (0.73 ± 0.23 across the observations), i.e.
glacier boundary layers warm ~0.83 °C per 1 °C of ambient warming. Regionally averaged
cooling reaches **−0.40 °C** relative to elevation-adjusted ERA5-Land, and up to 6.5 °C
locally under extreme warmth. The authors explicitly recommend that modelling efforts adopt
"a dynamic parameterization for glacier cooling, such as presented here."

Note the sign of a4: `k` *increases* with flowpath length in this formulation, which runs
opposite to SM10's `k2` decay with `DF`. The two are not the same quantity — SM10's `k2` is
conditional on `T_amb ≥ T*` at a point, whereas Shaw et al.'s `k` is an all-hours regression
slope per station-season, and the multiple regression partials out temperature, humidity and
wind. Don't mix coefficients between the two.

**Data available.** Per-glacier `k` for all RGI v6 glaciers (present-day and to 2099 under
SSP2-4.5/SSP5-8.5) is published as MATLAB `.mat` files on Zenodo:
<https://doi.org/10.5281/zenodo.14044846>. Zenodo rate-limited our automated fetch; worth
downloading manually. This would let us look up `k` by RGI ID with no fitting at all.

## Vapour pressure and other variables

Temperature is not the only variable affected, and the literature is much thinner here.

### Vapour pressure — SM10's scheme, now confirmed

This was the main open question, and the paper answers it. **SM10 does *not* assume constant
relative humidity.** It switches on the *on-glacier* temperature sign, not on `T*`:

```
e_gla = j1·e_amb + j2      for T_gla ≥ 0 °C     (surface at 0 °C: moisture exchange active)
e_gla = j3·e_amb + j4      for T_gla <  0 °C     (near-identity: j3 ≈ 1, j4 ≈ 0)
```
(SM10 Eq. 4; `e` in hPa.) The physics: with a melting surface, the surface vapour pressure is
pinned at saturation over ice at 0 °C = **6.11 hPa**. If `e_amb > 6.11` the gradient drives
condensation *onto* the surface, drawing moisture out of the KBL; if `e_amb < 6.11` it drives
evaporation/sublimation, adding moisture. So the `T_gla ≥ 0` regression line has slope < 1 and
**crosses the 1:1 line near 6.11 hPa** — a pivot, not a scaling. Below freezing the surface
vapour pressure is no longer pinned, but both `e_amb` and `e_surf` are small, so exchange is
minimal and the relation collapses to identity.

`j1`, `j2` are functions of flowpath length (SM10 Eqs. 8–9); `j3`, `j4` showed **no**
significant relation to any topographic index and are taken as weighted means:

```
j1 = β8  · exp(β9·FPL)          decays ~0.73 (FPL 1 km) → ~0.43 (10 km)   [Fig. 9a]
j2 = β10 + β11 · ln(FPL)        rises  ~1.5 hPa (1 km) → ~2.7 hPa (10 km) [Fig. 9b]
```

Per-site values (SM10 Table 6) — note `j1` shrinks and `j2` grows with FPL, together holding
the 6.11 hPa pivot roughly fixed:

| Site | FPL (m) | j1 | j2 (hPa) | R² | RMSE | j3 | j4 (hPa) | R² |
|---|---|---|---|---|---|---|---|---|
| BM1 | 10075 | 0.43 | 2.92 | 0.65 | 0.45 | 1.21 | −1.69 | 0.73 |
| PM1 | 2829 | 0.63 | 2.05 | 0.74 | 0.52 | 1.04 | −0.42 | 0.62 |
| PM2 | 2168 | 0.67 | 1.89 | 0.75 | 0.57 | 0.83 | 0.68 | 0.66 |
| PM3 | 1036 | 0.74 | 1.67 | 0.71 | 0.68 | 1.03 | 0.01 | 0.79 |
| PM4 | 970 | 0.67 | 1.83 | 0.70 | 0.64 | 0.89 | 0.49 | 0.77 |
| PM108 | 987 | 0.72 | 1.68 | 0.79 | 0.50 | 0.91 | 0.45 | 0.80 |
| PM308 | 2076 | 0.74 | 1.47 | 0.82 | 0.50 | 0.85 | 0.60 | 0.80 |
| WM1 | 3522 | 0.61 | 2.22 | 0.74 | 0.51 | 0.63 | 1.41 | 0.39 |
| WM2 | 1697 | 0.79 | 1.34 | 0.68 | 0.64 | 1.15 | −0.69 | 0.40 |

Independent-site vapour-pressure RMSE is 0.54–0.71 hPa with MBE −0.38 to +0.02 hPa. SM10 also
specifies its saturation formula explicitly — Tetens/Bolton (1980) with a phase switch:

```
e_s = 6.108 · 10^(9.5·T/(T+265.5))    T >  0 °C   (over ice)
e_s = 6.108 · 10^(7.5·T/(T+237.3))    T ≤  0 °C
```
which is worth comparing against the Magnus/Buck form already in `src/utils.jl`.

**Implication for us: constant RH is the wrong interim assumption.** `temperature_adjust`'s
constant-RH propagation would scale `e` down with `T`, but SM10 shows the KBL *adds* moisture
whenever ambient air is drier than 6.11 hPa — the opposite sign. In a cold, dry setting
(much of Antarctica/Greenland, where the surface is also usually below freezing) SM10's
sub-freezing branch is near-identity anyway, so the cleanest defensible choice is to
**leave `vapor_pressure` untouched** by the glacier correction rather than propagate ΔT
through it, and document that a full SM10 `e` scheme needs an FPL raster.

- **Longwave.** Cooling the air lowers incoming longwave. `temperature_adjust` already
  handles this via the Konzelmann clear-sky emissivity with the cloud increment preserved.
- **Longwave.** Cooling the air lowers incoming longwave. `temperature_adjust` already
  handles this via the Konzelmann clear-sky emissivity with the cloud increment preserved.
- **Wind speed.** Genuinely problematic. A katabatic layer produces a low-level jet at
  1–5 m above the surface that reanalysis cannot represent, and Shaw et al. (2024) show
  glacier winds can *enhance* turbulent heat transfer even while the air is cooler —
  partially offsetting the melt reduction from cooling alone. None of the schemes above
  correct wind speed. Recommend leaving wind untouched and documenting the caveat.

## What is implemented, and what is not

**Implemented (Shaw et al. 2025 decoupling factor).** `climate_adjust_for_glacier` in
`src/glacier_adjustment.jl`, reading `k` from the published per-glacier lookup table via
`glacier_decoupling` in `src/glacier_decoupling.jl`.

The `.mat` is *not* read at runtime: `data/make_shaw2025_decoupling.jl` converts it once to
`data/shaw2025_glacier_decoupling.csv.gz` (2.9 MB, 186,792 glaciers), so HDF5.jl stays out of
`Project.toml` and the lookup needs no network. `K_CI` is carried through as `k_lower` — note
it is an absolute **lower bound**, not a half-width, and is unclamped so it can be negative.
The five-predictor regression is deliberately not re-implemented, for the reasons above
(Table S2 sign errors in `a4`/`a5`; a corrected refit still only reaches R² = 0.54 against the
authors' own `k`).

As applied:

```
T′ = T + (k − 1)·max(T − T_ref, 0)          T_ref = 273.15 K
```

written as an increment so `k = 1` is bit-exact, and gated at the melting point by default
(`apply_below_freezing=false`) — a bare `k` multiplier applied to sub-freezing air *warms* it,
which has no physical basis, so the gate reduces the scheme to a threshold form in the spirit
of SM10. Longwave is recomputed from the cooled temperature at unchanged `e` (Konzelmann
clear-sky emissivity, cloud increment preserved); `vapor_pressure` and `wind_speed` are left
alone, per the two caveats below. Metadata records a cumulative multiplicative
`glacier_decoupling_factor` plus the lookup provenance, and
`validate_climate_forcing_units` re-runs on the result.

**Not implemented:** the `:sm10` path, which is the only option where the `k` table has no
coverage (regions 05 and 19). Until then, those regions require passing `k` explicitly — a
regional mean, or a value from another scheme.

The `:sm10` scheme (Shea & Moore 2010 + Carturan `T*`) would add per-glacier flowline
behaviour and a genuine hourly threshold response. It needs a flowpath-length raster, which is
real work: SM10 used SAGA's Quinn et al. (1991) multiple-flow-direction algorithm, the Shaw
papers use MATLAB TopoToolbox, and Shaw et al. (2021) flag that flowline generation is *not*
standardised between studies and is DEM-quality dependent. That is the main implementation
cost. Verify any implementation against SM10's Table 5/Table 6 site values (both reproduced
above) before trusting it, and follow the same conventions as the implemented path:
scalar/lookup arguments, cumulative metadata so repeated calls compose, validation re-run on
the result, and an exact identity at `k = 1`.

Three caveats, all surfaced in the `climate_adjust_for_glacier` docstring:

- **Order of operations.** Elevation downscaling (`climate_adjust_for_elevation`) must come
  first — it moves the forcing to the target elevation to give `T_amb` at the glacier
  surface — and the glacier correction second. Applying them in the reverse order is wrong,
  since the GBL correction is defined relative to ambient temperature *at the on-glacier
  elevation*.
- **Seasonality.** All of these schemes are calibrated on ablation-season data. Outside the
  melt season, and whenever the surface is well below 0 °C, the physical basis disappears;
  the threshold form of SM10 handles this gracefully (`k1 ≈ 0.9–1.0`, near-identity),
  whereas a bare `k` multiplier does not. Hence the melting-point gate: the correction is
  applied only above `T_ref`, and `apply_below_freezing=true` opts into the ungated form
  (which warms sub-freezing air) for sensitivity testing only.
- **Do not reuse `temperature_adjust`'s vapour-pressure propagation here.** Constant RH is
  the wrong model for a glacier boundary layer: SM10 Eq. 4 shows `e_gla` pivots about
  6.11 hPa, so the KBL *adds* moisture when ambient air is drier than that — the opposite
  sign to constant-RH scaling. Leave `vapor_pressure` untouched unless the full SM10 `e`
  scheme (and hence an FPL raster) is implemented. Longwave should still be recomputed from
  the cooled temperature.

## Still missing

Both blocking items from the first research pass have since been supplied by the user and are
now held locally — Shea & Moore (2010) and the Shaw et al. (2025) Zenodo datasets. One
non-blocking gap remains.

1. **Ayala, A., Pellicciotti, F. & Shea, J. M. (2015)**, *Modeling 2 m air temperatures over
   mountain glaciers: Exploring the influence of katabatic cooling and external warming*,
   JGR Atmospheres 120, 3139–3157. doi:10.1002/2015JD023137
   - Closed at Wiley. An accepted manuscript is indexed at Northumbria Research Link but the
     host no longer resolves; the listed URL was
     <https://nrl.northumbria.ac.uk/id/eprint/22863/1/jgrd52075.pdf>
   - **Not needed.** ModGB's equations and parameter values are already in hand from Shaw et
     al. (2017) and Troxler et al. (2020), and we are not recommending ModGB.

2. **`FUTURE_DECOUPLING_ESTIMATES_DATABASE.mat`** — the third Zenodo file (SSP2-4.5 /
   SSP5-8.5 `k` and cooling to 2099, static and dynamic glacier geometry) was not among the
   supplied files. Only relevant if we ever want projected `k`; the historical file is
   sufficient for forcing reanalysis. Source: <https://doi.org/10.5281/zenodo.14044846>

3. **Greuell & Böhm (1998)**, J. Glaciol. 44, 9–20, doi:10.3189/S0022143000002306 — the
   original GB derivation. Paywalled, and not needed for the recommended path.

## References held locally

| File | Paper |
|---|---|
| `Shaw2025_NCC_glaciers_recouple.pdf` (+ `_supplement`) | Shaw et al. 2025, Nat. Clim. Chang. — global `k`, coefficient table in supplement Table S2 |
| `Shaw2021_TC_distributed_Ta_Tibet.pdf` | Shaw et al. 2021, The Cryosphere — SM10 equations + Table 3 coefficient sets, global `k1`/`k2` synthesis |
| `Carturan2015_TC_Ortles_Cevedale.pdf` | Carturan et al. 2015 — the `T*(DF)` reformulation; SM10 vs GB intercomparison |
| `Shaw2017_JoG_Tsanteleina.pdf` | Shaw et al. 2017 — SM10/ModGB equations, Table 5 parameters, melt-model impact |
| `Troxler2020_JoG_McCall.pdf` | Troxler et al. 2020 — ModGB equations 1–4, decadal transferability test |
| `Shaw2024_JGR_local_controls_cooling.pdf` | Shaw et al. 2024 — why size alone doesn't predict cooling; valley/synoptic wind intrusions |
| `Shaw2023_GRL_decaying_boundary_layer.pdf` | Shaw et al. 2023 — retreat weakens the GBL, raising sensitivity |
| `Sauter2026_RoG_glacier_atmosphere_review.pdf` | Sauter et al. 2026, Rev. Geophys. — current review; confirms reconstructing on-glacier climate without measurements remains an open problem |
| `Shea2010_JGR_distributed_Ta_vapor_pressure.pdf` | **Shea & Moore 2010** — the original SM10. Eq. 1 (T), Eq. 4 (vapour pressure), Eqs. 5–9 (FPL transfer functions), Tables 5–7 site parameters, Figs. 7/9 the β curves |

Datasets in `docs/references/shaw2025_decoupling/` (Shaw et al. 2025, Zenodo
[10.5281/zenodo.14044846](https://doi.org/10.5281/zenodo.14044846)):

| File | Contents |
|---|---|
| `HISTORICAL_DECOUPLING_ESTIMATES_DATABASE.mat` | 186 792 RGI v6 glaciers × (`RGI_ID`, `REGION`, `LAT`, `LON`, `ELE`, `TA`, `Q`, `FF`, `LEN`, `DEB`, `K`, `K_CI`, `COOL`, `COOL_CI`). **This is the table a `:decoupling` implementation should read.** MATLAB v7.3 = HDF5, so `HDF5.jl` suffices; `RGI_ID` is an object-reference array that dereferences to `"RGI60-01.00001"`-style strings |
| `GLACIER_DECOUPLING_OBSERVATIONS_DATABASE.mat` | 186 glacier-years / 415 AWS points — the observational training set, incl. `AWS_K`, `AWS_R2`, `AWS_K_UNC`, `AWS_FPL`. Good validation target |
| `DECOUPLING_COOLING_DATABASE_README.txt` | Authoritative variable definitions (and the source of the `FF`×2.5 and `Q` unit notes above) |
| `METADATA_TABLE_ON-GLACIER_TA_DATA.csv` + `METADATA_TABLE_REFERENCES.txt` | Supplementary Table S1 — per-glacier provenance and references for the on-glacier records |

Note `docs/references/.gitignore` excludes the PDFs and `.mat` files (~135 MB) from version
control; this file records what they are and where to get them.
