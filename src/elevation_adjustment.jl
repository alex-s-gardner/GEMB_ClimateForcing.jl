"""
Elevation adjustment of climate forcing.

Adjusts a climate-forcing `DimStack` for the difference in elevation between the
reanalysis grid cell and a desired target elevation. This is the "downscaling"
step required when a point of interest (e.g. an AWS, a glacier stake, or a
high-resolution DEM cell) sits at a different elevation than the coarse
reanalysis surface.

# Scientific basis

The scheme follows current best practice for physically-based downscaling of
*surface* reanalysis fields over glaciers and ice sheets. The per-variable
corrections closely follow Glover (1999, *J. Climate* 12, 551–563, Eqs. 15–20),
who established this exact set of elevation adjustments — temperature lapse,
hydrostatic pressure, constant-relative-humidity moisture, emissivity-based
longwave, and Clausius–Clapeyron precipitation — for downscaling GCM output to
the Greenland ice sheet. The individual adjustments are consistent with those
validated over the Greenland ice sheet by the RACMO2.3p2 studies of Noël et al.
(2018, *The Cryosphere* 12, 811–831; 2019, *Sci. Adv.* 5, eaaw0123). It is the
surface-field analogue of TopoSCALE (Fiddes & Gruber, 2014, *Geosci. Model Dev.*
7, 387–405); the full pressure-level ("3D") interpolation of TopoSCALE is not
possible from ERA5-Land, which exposes surface variables only.

Let `Δz = z_target − z_reanalysis` (metres; positive when the target is *above*
the reanalysis surface).

1. **Air temperature** — environmental lapse rate: `T′ = T − (Γ/1000)·Δz`, with
   `Γ` the lapse rate in K/km (positive = cooling with height). Over melting
   glaciers/ice sheets observed near-surface lapse rates are markedly shallower
   than the free-air 6.5 K/km (a surface pinned near 0 °C plus persistent
   inversions/katabatic flow decouple near-surface air from the free atmosphere;
   Fausto et al. 2009; Gardner et al. 2009), except the Antarctic interior which
   is steeper. Region-specific monthly tables are provided as named constants
   (see below) to pass via `lapse_rate`.

2. **Air pressure** — hydrostatic (barometric) formula with the mean layer
   temperature: `P′ = P·exp(−g·Δz / (R_d·T̄))`, `T̄ = ½(T + T′)`. Noël et al.
   (2018) show uncorrected surface-pressure bias scales directly with the
   grid/target elevation difference, motivating this correction. (Glover 1999
   Eq. 16 uses the equivalent first-order linear form.)

3. **Vapor pressure** — *constant relative humidity*: RH is diagnosed at the
   reanalysis level, preserved with elevation, and vapor pressure is recomputed
   from the saturation curve at the adjusted temperature. This matches Glover
   (1999, Eq. 20: `q′ = q·q_sat(T′)/q_sat(T)`) and the way RACMO2 diagnoses 2 m
   specific humidity from T, P and RH (Curry & Webster 1999, per Noël et al. 2018).

4. **Longwave downward** — recomputed from the adjusted temperature and vapor
   pressure using the Konzelmann et al. (1994) clear-sky emissivity (developed on
   the Greenland ice sheet), while *preserving the reanalysis cloud emissivity
   increment* Δε. Noël et al. (2019) show polar melt is governed by cloud-driven
   longwave, so the cloud signal must be retained rather than reduced to clear
   sky (the TopoSCALE longwave scheme). Glover (1999, Eqs. 17–18) uses the
   simpler total-emissivity rescaling `LW′ = [LW/(σT⁴)]·σ(T′)⁴`.

5. **Shortwave downward** and **wind speed** — left unchanged. Solar
   elevation-attenuation is ≲ a few %/km and requires solar-geometry/terrain
   information not available here; RACMO-style SMB downscaling likewise carries
   shortwave through unmodified.

6. **Precipitation** — left unchanged by default (`precip_scaling_method =
   nothing`). This follows the RACMO2.3 statistical-downscaling studies of Noël
   et al. (2016, 2020, 2022, 2025), which interpolate total precipitation onto
   the fine grid *without* an elevation correction, relying instead on the
   orographic response resolved by the host model. With
   `precip_scaling_method = :clausius_clapeyron`, precipitation is scaled by the
   saturation ratio `eₛ(T′)/eₛ(T)` (Glover 1999, Eq. 19), so it decreases as air
   cools with elevation — a physically-based representation of the ice-sheet
   "elevation-desert" effect.

   Note: although *total* precipitation is preserved by default, its rain/snow
   *phase* is affected downstream. Because temperature is elevation-cooled, a
   temperature-based partitioning (as GEMB applies) will assign a larger snow
   fraction at higher, colder target elevations.

See [`climate_adjust_for_elevation`](@ref).
"""

using DimensionalData
using Dates
using Statistics

# ----------------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------------
const _G0 = 9.80665          # standard gravity, m/s²
const _R_DRY = 287.05        # specific gas constant of dry air, J/(kg·K)
const _SIGMA_SB = 5.67e-8    # Stefan–Boltzmann constant, W/m²/K⁴

# ----------------------------------------------------------------------------
# Region-specific monthly near-surface temperature lapse rates (K/km, POSITIVE =
# cooling with height), index 1..12 = Jan..Dec. Pass one of these to the
# `lapse_rate` keyword of `climate_adjust_for_elevation`. All are for snow/ice
# (glacier / ice-sheet) surfaces.
# ----------------------------------------------------------------------------

"""
Greenland ice-sheet mean monthly near-surface lapse rates (K/km, positive =
cooling with height), from Fausto et al. (2009), *J. Glaciol.* 55(189), 95–105,
Table 4 (annual mean 6.8 K/km). Shallow in summer (melting surface pinned near
0 °C), steep in winter. Fausto's table is signed positive = warming with height;
these are negated to the package convention (positive = cooling). The
parameterization does not represent the coastal inversion layer. Pass to
`climate_adjust_for_elevation`'s `lapse_rate`.
"""
const GREENLAND_LAPSE_RATE = [7.9, 8.9, 7.9, 7.3, 5.9, 4.7, 4.6, 5.7, 6.9, 7.3, 6.5, 7.6]

"""
Arctic-glacier monthly near-surface lapse rates (K/km), Gardner et al. (2009),
*J. Climate* 22, 4281–4298. Ablation-season mean ≈ 4.9 K/km, winter ≈ 3.2 K/km.
"""
const ARCTIC_LAPSE_RATE = [3.2, 3.2, 3.5, 4.0, 4.5, 4.9, 4.9, 4.9, 4.5, 4.0, 3.5, 3.2]

"""
Antarctic **interior/plateau (>1500 m a.s.l.)** near-surface lapse rate (K/km),
Fortuin & Oerlemans (1990), *Ann. Glaciol.* 14, 78–84, Table II (A = −14.285 ±
0.645 K km⁻¹). Intense radiative cooling makes this super-adiabatic and roughly
year-round constant.

!!! warning
    This value applies to the high interior only. Fortuin & Oerlemans give very
    different regimes elsewhere: escarpment (200–1500 m) ≈ 5.1 K/km, and ice
    shelves (<200 m) ≈ 0 (temperature is latitude-controlled, not
    elevation-controlled). Do **not** apply 14.3 K/km to coastal/escarpment
    Antarctica — it will badly over-extrapolate. Prefer [`empirical_lapse_rate`](@ref)
    there.
"""
const ANTARCTICA_LAPSE_RATE = fill(14.3, 12)

"""
    saturation_vapor_pressure(T_K::Real; over_ice::Bool) -> Float64

Saturation vapor pressure (Pa) at temperature `T_K` (Kelvin), using the Buck
(1981) coefficients over water or over ice.

Over water: `eₛ = 611.21·exp(17.502·Tc/(240.97+Tc))`
Over ice:   `eₛ = 611.15·exp(22.452·Tc/(272.55+Tc))`

with `Tc` in °C. If `over_ice` is not supplied it defaults to `T_K < 273.15`.
These are the exact Buck (1981) coefficients as tabulated by Liston & Elder
(2006, MicroMet, Eq. 4).

# References
- Buck, A. L. (1981). New equations for computing vapor pressure and enhancement
  factor. *J. Appl. Meteorol.* 20, 1527–1532.
- Liston, G. E. & Elder, K. (2006). MicroMet. *J. Hydrometeorol.* 7, 217–234 (Eq. 4).
"""
function saturation_vapor_pressure(T_K::Real; over_ice::Bool=T_K < 273.15)
    Tc = Float64(T_K) - 273.15
    if over_ice
        return 611.15 * exp(22.452 * Tc / (272.55 + Tc))
    else
        return 611.21 * exp(17.502 * Tc / (240.97 + Tc))
    end
end

"""
    konzelmann_clear_sky_emissivity(e_Pa::Real, T_K::Real) -> Float64

Clear-sky atmospheric emissivity from the Konzelmann et al. (1994)
parameterization, developed on the Greenland ice sheet:

`ε_cs = 0.23 + 0.484·(e/T)^(1/8)`

**Vapor pressure `e_Pa` must be in Pascals** (the most common implementation
error is passing hPa/kPa; a physical check at T≈263 K, e≈400 Pa gives ε_cs≈0.74,
whereas hPa would give an unphysical ≈0.52). `T_K` in Kelvin. The coefficient
0.484 is the original Konzelmann et al. (1994) Greenland-fitted value; some later
implementations use 0.443.

# References
- Konzelmann, T. et al. (1994). Parameterization of global and longwave incoming
  radiation for the Greenland ice sheet. *Glob. Planet. Change* 9, 143–164.
"""
function konzelmann_clear_sky_emissivity(e_Pa::Real, T_K::Real)
    return 0.23 + 0.484 * (Float64(e_Pa) / Float64(T_K))^(1 / 8)
end

"""
    empirical_lapse_rate(values::AbstractVector, elevations::AbstractVector) -> Float64

Estimate a *local* near-surface lapse rate (K/km, positive = cooling with
height) as the negative ordinary-least-squares slope of `values` (e.g. air
temperatures) regressed on `elevations` (metres): fitting
`value = a·elevation + b`, this returns `−1000·a`.

This is the gold-standard downscaling approach used by RACMO2.3's statistical
downscaling (Noël et al. 2016; Noël et al. 2025, *Nat. Commun.* Eq. 1), which
derives elevation gradients per grid cell from the cell and its ≥6 neighbouring
cells rather than assuming a fixed climatological rate. Where such neighbouring
reanalysis samples are available, pass the returned value to
[`climate_adjust_for_elevation`](@ref) via `lapse_rate` for a locally-fitted
gradient.

At least two samples with non-zero elevation spread are required; otherwise an
`ArgumentError` is thrown.

# Example
```julia
# Temperatures (K) and elevations (m) at a grid cell and its 8 neighbours.
Γ = empirical_lapse_rate(neighbour_T2m, neighbour_elevations)   # K/km, positive
adjusted = climate_adjust_for_elevation(stack, Δz; lapse_rate=Γ)
```

# References
- Noël, B. et al. (2016). Downscaled Greenland SMB (1958–2015). *The Cryosphere* 10, 2361–2377.
- Noël, B. et al. (2025). Poleward shift of subtropical highs drives Patagonian
  glacier mass loss. *Nat. Commun.* 16, 3795.
"""
function empirical_lapse_rate(values::AbstractVector, elevations::AbstractVector)
    length(values) == length(elevations) ||
        throw(ArgumentError("values and elevations must have equal length"))
    n = length(values)
    n >= 2 || throw(ArgumentError("need at least 2 samples to fit a lapse rate"))
    z̄ = Statistics.mean(elevations)
    v̄ = Statistics.mean(values)
    Szz = zero(Float64)
    Szv = zero(Float64)
    @inbounds for i in eachindex(values, elevations)
        dz = elevations[i] - z̄
        Szz += dz * dz
        Szv += dz * (values[i] - v̄)
    end
    Szz > 0 || throw(ArgumentError("elevations have zero spread; cannot fit a lapse rate"))
    # Regression slope is K/m (negative when cooling with height); return the
    # conventional lapse rate in K/km, positive for cooling with height.
    return -1000.0 * (Szv / Szz)
end

# Physically plausible bounds for a near-surface lapse rate (K/km, positive =
# cooling with height). The upper bound comfortably admits the Antarctic-interior
# super-adiabatic value (≈14.3 K/km; Fortuin & Oerlemans 1990) with headroom; the
# negative bound admits strong near-surface temperature inversions (temperature
# increasing with height), which are common over polar/glacier surfaces
# (Gardner et al. 2009 report instantaneous values down to ≈ −12 K/km). Values
# outside this range almost always indicate a unit error (e.g. K/m instead of
# K/km) or a sign mistake.
const _LAPSE_RATE_MIN = -30.0   # K/km
const _LAPSE_RATE_MAX =  25.0   # K/km

"""
    _resolve_lapse_rate(lapse_rate, months, n) -> Vector{Float64}

Expand the user `lapse_rate` argument into an `n`-element per-timestep vector
(K/km). A scalar broadcasts to all steps; a length-12 vector is treated as
monthly, assumed ordered January (index 1) to December (index 12) and indexed by
`months`; a length-`n` vector is used per-timestep.

If the record itself has length 12 the length-12 case is interpreted as monthly.

All resolved values are checked to lie within physically plausible bounds
(`$(_LAPSE_RATE_MIN)` to `$(_LAPSE_RATE_MAX)` K/km); an out-of-range value raises
an `ArgumentError` (most often caused by passing K/m instead of K/km, or a wrong
sign).
"""
function _resolve_lapse_rate(lapse_rate, months::AbstractVector{<:Integer}, n::Integer)
    Γ = if lapse_rate isa Real
        fill(Float64(lapse_rate), n)
    else
        L = length(lapse_rate)
        if L == 12
            monthly = Float64.(collect(lapse_rate))
            monthly[months]
        elseif L == n
            Float64.(collect(lapse_rate))
        else
            throw(ArgumentError(
                "lapse_rate must be a scalar, a 12-element monthly vector, or a " *
                "vector matching the $(n)-step climate record (got length $(L))"))
        end
    end

    # Guard against unit (K/m vs K/km) and sign errors.
    lo, hi = extrema(Γ)
    if lo < _LAPSE_RATE_MIN || hi > _LAPSE_RATE_MAX
        throw(ArgumentError(
            "lapse_rate out of physical range: got [$(lo), $(hi)] K/km, expected " *
            "within [$(_LAPSE_RATE_MIN), $(_LAPSE_RATE_MAX)] K/km. Note lapse_rate " *
            "is in K/km (positive = cooling with height); a value near zero often " *
            "means K/m was passed by mistake."))
    end
    return Γ
end

"""
    climate_adjust_for_elevation(climate_forcing_original::DimStack,
                                 delta_elevation::Real;
                                 lapse_rate=6.5,
                                 precip_scaling_method=nothing) -> DimStack

Adjust a climate-forcing `DimStack` for a difference in elevation between the
reanalysis grid cell and the desired target elevation (downscaling).

# Arguments
- `climate_forcing_original::DimStack`: forcing as returned by
  [`climate_forcing`](@ref), with a `Ti` (time) dimension and variables
  `temperature_air`, `pressure_air`, `vapor_pressure`, `wind_speed`,
  `precipitation`, `shortwave_downward`, `longwave_downward`.
- `delta_elevation::Real`: `z_target − z_reanalysis` in metres. Positive when the
  target elevation is *above* the reanalysis surface (variables are cooled,
  pressure reduced). Negative brings the forcing *down* to a lower target.

# Keyword arguments
- `lapse_rate=6.5`: near-surface temperature lapse rate in **K/km**, positive for
  cooling with height. Accepts:
  - a **scalar** (applied to every time step; default `6.5`, the free-air value);
  - a **length-12 vector** — treated as monthly values, assumed ordered from
    January (index 1) to December (index 12), indexed by the month of each time step;
  - a **vector matching the length of the climate record** — used as a
    per-time-step lapse rate.

  (If the record itself has exactly 12 steps, a length-12 vector is interpreted
  as monthly.) Region-specific monthly tables are provided as named constants:
  [`GREENLAND_LAPSE_RATE`](@ref), [`ARCTIC_LAPSE_RATE`](@ref),
  [`ANTARCTICA_LAPSE_RATE`](@ref). Use [`empirical_lapse_rate`](@ref) to fit the
  rate locally from neighbouring grid cells (the RACMO/gold-standard approach).
  All resolved values must lie within physically plausible bounds (−30 to 25
  K/km); an out-of-range value raises an `ArgumentError` (usually a K/m-vs-K/km
  unit slip or a sign error).
- `precip_scaling_method=nothing`: precipitation elevation treatment.
  - `nothing` (default) — precipitation unchanged, matching the RACMO2.3
    downscaling studies (Noël et al. 2016–2025).
  - `:clausius_clapeyron` — scale precipitation by `eₛ(T′)/eₛ(T)` (Glover 1999,
    Eq. 19): precipitation decreases as air cools with elevation (the ice-sheet
    "elevation-desert" effect).

# Returns
- A new `DimStack` on the same time dimension with adjusted `temperature_air`,
  `pressure_air`, `vapor_pressure`, `longwave_downward`, and (if a precipitation
  scaling is requested) `precipitation`. `wind_speed` and `shortwave_downward`
  are carried through unchanged. Metadata gains `delta_elevation`, updated
  elevation, and recomputed means; physical-range validation is re-run.

# Method
For `Δz = delta_elevation` and lapse rate `Γ` (K/km):
- `T′ = T − (Γ/1000)·Δz`
- `P′ = P·exp(−g·Δz/(R_d·T̄))`, `T̄ = ½(T+T′)`, `g=9.80665`, `R_d=287.05`
- `RH = e/eₛ(T)`; `e′ = clamp(RH, 0, 1)·eₛ(T′)` (over-ice curve below 0 °C)
- `Δε = LW/(σ·T⁴) − ε_cs(e,T)`; `LW′ = (ε_cs(e′,T′) + Δε)·σ·T′⁴`, `σ=5.67e-8`,
  `ε_cs` from Konzelmann et al. (1994)
- `precip′ = precip·eₛ(T′)/eₛ(T)` if `precip_scaling_method=:clausius_clapeyron`,
  else unchanged

`Δz = 0` reproduces the input exactly.

# Example
```julia
stack = climate_forcing(:era5land, 72.58, -38.46;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

# Target sits 250 m above the grid cell; use Greenland monthly lapse rates.
adjusted = climate_adjust_for_elevation(stack, 250.0; lapse_rate=GREENLAND_LAPSE_RATE)

# Or a single locally-observed rate, plus elevation-desert precipitation scaling.
adjusted = climate_adjust_for_elevation(stack, 250.0;
                                        lapse_rate=5.5,
                                        precip_scaling_method=:clausius_clapeyron)
```

# References
- Glover, R. W. (1999). Influence of spatial resolution and treatment of
  orography on GCM estimates of the surface mass balance of the Greenland ice
  sheet. *J. Climate* 12, 551–563. (Eqs. 15–20: the adjustment scheme.)
- Fiddes, J. & Gruber, S. (2014). TopoSCALE v.1.0. *Geosci. Model Dev.* 7, 387–405.
- Fausto, R. S. et al. (2009). *J. Glaciol.* 55(189), 95–105.
- Gardner, A. S. et al. (2009). *J. Climate* 22, 4281–4298.
- Fortuin, J. P. F. & Oerlemans, J. (1990). *Ann. Glaciol.* 14, 78–84.
- Konzelmann, T. et al. (1994). *Glob. Planet. Change* 9, 143–164.
- Noël, B. et al. (2016). *The Cryosphere* 10, 2361–2377.
- Noël, B. et al. (2018). *The Cryosphere* 12, 811–831.
- Noël, B. et al. (2019). *Sci. Adv.* 5, eaaw0123.
- Noël, B. et al. (2020). *Nat. Commun.* 11, 4597.
- Noël, B. et al. (2022). *Geophys. Res. Lett.* 49, e2021GL095697.
- Noël, B. et al. (2025). *Nat. Commun.* 16, 3795.
"""
function climate_adjust_for_elevation(
    climate_forcing_original::DimStack,
    delta_elevation::Real;
    lapse_rate=6.5,
    precip_scaling_method::Union{Symbol,Nothing}=nothing,
)
    Δz = Float64(delta_elevation)

    # Validate precipitation scaling method (only snow/ice-relevant option).
    if !(precip_scaling_method === nothing || precip_scaling_method === :clausius_clapeyron)
        throw(ArgumentError(
            "precip_scaling_method must be nothing or :clausius_clapeyron (got $(precip_scaling_method))"))
    end

    # Per-timestep month index from the time dimension.
    times = collect(dims(climate_forcing_original, Ti))
    months = month.(times)
    n = length(times)

    # Resolve lapse rate to a per-timestep vector (K/km, positive = cooling up).
    Γ = _resolve_lapse_rate(lapse_rate, months, n)

    # Source variables as DimArrays (already Float64 from the loader). Broadcasting
    # over them preserves the time dimension, so there is no need to strip with
    # parent() and rebuild the dims by hand.
    T  = climate_forcing_original[:temperature_air]
    P  = climate_forcing_original[:pressure_air]
    e  = climate_forcing_original[:vapor_pressure]
    LW = climate_forcing_original[:longwave_downward]
    precip = climate_forcing_original[:precipitation]

    # 1. Temperature: environmental lapse (Γ in K/km, Δz in m ⇒ /1000).
    T′ = @. T - (Γ / 1000.0) * Δz

    # 2. Pressure: hydrostatic with mean layer temperature.
    P′ = @. P * exp(-_G0 * Δz / (_R_DRY * 0.5 * (T + T′)))

    # 3. Vapor pressure: constant relative humidity (over-ice curve below 0 °C).
    #    Using the same saturation formula to diagnose and reconstruct RH makes
    #    Δz = 0 an exact identity. Clamp RH to (0, 1] to avoid supersaturation.
    es_T  = saturation_vapor_pressure.(T)
    es_T′ = saturation_vapor_pressure.(T′)
    RH = @. clamp(e / es_T, 0.0, 1.0)
    e′ = @. RH * es_T′

    # 4. Longwave down: recompute from adjusted T, e; preserve cloud increment Δε.
    LW′ = _adjust_longwave(LW, e, T, e′, T′)

    # 5. Precipitation: unchanged by default; Clausius–Clapeyron scaling if asked
    #    (Glover 1999, Eq. 19). eₛ ratio proxies q_sat(T′)/q_sat(T); exact at Δz=0.
    precip′ = precip_scaling_method === :clausius_clapeyron ? (@. precip * (es_T′ / es_T)) : precip

    src_meta = metadata(climate_forcing_original)
    z_reanalysis = get(src_meta, "elevation", nothing)
    new_meta = merge(copy(src_meta), Dict(
        "delta_elevation" => Δz,
        # Cumulative elevation offset relative to the source reanalysis, so
        # repeated adjustments compose. Surfaced downstream (via `initialize_forcing`)
        # in GEMB output metadata and diagnostic plots.
        "elevation_offset" => get(src_meta, "elevation_offset", 0.0) + Δz,
        "precip_scaling_method" => string(precip_scaling_method),
        "temperature_air_mean" => Statistics.mean(T′),
        "precipitation_mean" => Statistics.mean(precip′) * 8760.0,
    ))
    if !isnothing(z_reanalysis)
        new_meta["elevation_reanalysis"] = z_reanalysis
        new_meta["elevation"] = z_reanalysis + Δz
    end

    # Wind speed and shortwave are carried through untouched.
    adjusted = _rebuild_forcing(climate_forcing_original, new_meta;
                                temperature_air=T′, pressure_air=P′, vapor_pressure=e′,
                                precipitation=precip′, longwave_downward=LW′)

    # Re-validate physical ranges after adjustment.
    validate_climate_forcing_units(adjusted)

    return adjusted
end
