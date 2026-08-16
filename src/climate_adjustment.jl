"""
Direct climate perturbations of a forcing `DimStack`.

Where [`climate_adjust_for_elevation`](@ref) perturbs a forcing for a change in
*elevation* (and derives the temperature change from a lapse rate), the two
functions here apply a perturbation that is prescribed directly:

1. [`temperature_adjust`](@ref) — a uniform temperature offset `ΔT` (K), as used
   for warming/cooling sensitivity experiments (`+2 K`) or for bias-correcting a
   reanalysis against an in-situ record. The offset is propagated into the
   variables that physically depend on air temperature: vapor pressure (holding
   relative humidity fixed) and downwelling longwave irradiance (through the
   Konzelmann et al. (1994) clear-sky emissivity and `σT⁴`).

2. [`precipitation_adjust`](@ref) — a fractional rescaling of precipitation
   (`1.15` = 15 % wetter). Nothing else in the stack depends on the precipitation
   rate, so this is a single-variable change.

The two are deliberately independent: temperature perturbations do *not* apply a
Clausius–Clapeyron precipitation response, so a combined scenario is expressed by
composing them, and the size of the precipitation response stays an explicit
choice rather than an implicit consequence of `ΔT`.

Both share the conventions of the elevation path: the identity perturbation
(`ΔT = 0`, `scaling = 1`) reproduces the input exactly, unperturbed variables are
carried through untouched, physical-range validation is re-run on the result, and
the perturbation is recorded cumulatively in the stack metadata so repeated calls
compose.
"""

using DimensionalData
using Statistics

# ----------------------------------------------------------------------------
# Shared helpers, all reused rather than redefined: `_SIGMA_SB`,
# `saturation_vapor_pressure` and `konzelmann_clear_sky_emissivity` from
# elevation_adjustment.jl; `_rebuild_forcing` and `_adjust_longwave` from
# utils.jl; `validate_climate_forcing_units` from datasets/era5_land.jl.
# ----------------------------------------------------------------------------

"""
    temperature_adjust(climate_forcing_original::DimStack,
                       delta_temperature::Real) -> DimStack

Apply a uniform air-temperature offset to a climate-forcing `DimStack`, and
propagate it through the variables that depend on temperature.

Use this for temperature sensitivity experiments (e.g. a `+2 K` warming scenario)
or to bias-correct a reanalysis against an observed temperature record. For a
temperature change that arises from a *difference in elevation*, use
[`climate_adjust_for_elevation`](@ref) instead — it derives `ΔT` from a lapse rate
and additionally corrects surface pressure.

# Arguments
- `climate_forcing_original::DimStack`: forcing as returned by
  [`climate_forcing`](@ref), with a `Ti` (time) dimension and variables
  `temperature_air`, `pressure_air`, `vapor_pressure`, `wind_speed`,
  `precipitation`, `shortwave_downward`, `longwave_downward`.
- `delta_temperature::Real`: temperature offset in **K**, added to every time
  step. Positive = warming, negative = cooling.

# Returns
- A new `DimStack` on the same time dimension with adjusted `temperature_air`,
  `vapor_pressure` and `longwave_downward`. `pressure_air`, `wind_speed`,
  `precipitation` and `shortwave_downward` are carried through unchanged.
  Metadata gains `delta_temperature` and a cumulative `temperature_offset`, and
  `temperature_air_mean` is recomputed; physical-range validation is re-run.

# Method
For `ΔT = delta_temperature`:
- `T′ = T + ΔT`
- `RH = e/eₛ(T)`; `e′ = clamp(RH, 0, 1)·eₛ(T′)` — **relative humidity is held
  constant** (over-ice saturation curve below 0 °C), so specific humidity rises
  with warming
- `Δε = LW/(σ·T⁴) − ε_cs(e,T)`; `LW′ = (ε_cs(e′,T′) + Δε)·σ·T′⁴`, `σ=5.67e-8`,
  `ε_cs` from Konzelmann et al. (1994) — the cloud/aerosol emissivity increment
  `Δε` diagnosed from the input is preserved, so only the clear-sky response and
  the `T⁴` scaling change

Surface pressure is *not* adjusted: it is set by elevation and synoptic mass
distribution, not by a uniform warming. Wind speed and shortwave irradiance are
likewise left untouched.

`ΔT = 0` reproduces the input exactly.

!!! note
    A large negative `ΔT` can push `temperature_air` below the 180 K validation
    floor, or `longwave_downward` below 50 W/m². That raises an `ArgumentError`
    from `validate_climate_forcing_units` by design — it means the perturbation
    has left the physically plausible range of the parameterizations.

# Example
```julia
stack = climate_forcing(:era5land, 72.58, -38.46;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

# +2 K warming scenario (vapor pressure and longwave respond).
warmer = temperature_adjust(stack, 2.0)

# Combine with a wetter precipitation scenario.
warmer_wetter = precipitation_adjust(warmer, 1.15)
```

# References
- Buck, A. L. (1981). New equations for computing vapor pressure and enhancement
  factor. *J. Appl. Meteorol.* 20, 1527–1532. (Saturation vapor pressure.)
- Konzelmann, T. et al. (1994). Parameterization of global and longwave incoming
  radiation for the Greenland ice sheet. *Glob. Planet. Change* 9, 143–164.
- Glover, R. W. (1999). *J. Climate* 12, 551–563. (Same constant-RH and
  emissivity treatment, in an elevation-downscaling context.)

See also [`precipitation_adjust`](@ref), [`climate_adjust_for_elevation`](@ref).
"""
function temperature_adjust(climate_forcing_original::DimStack, delta_temperature::Real)
    ΔT = Float64(delta_temperature)

    # Source variables as DimArrays. Broadcasting preserves the time dimension.
    T  = climate_forcing_original[:temperature_air]
    e  = climate_forcing_original[:vapor_pressure]
    LW = climate_forcing_original[:longwave_downward]

    # 1. Temperature: uniform offset.
    T′ = @. T + ΔT

    # 2. Vapor pressure: constant relative humidity (over-ice curve below 0 °C).
    #    Diagnosing and reconstructing RH with the same saturation formula makes
    #    ΔT = 0 an exact identity. Clamp RH to [0, 1] to avoid supersaturation.
    es_T  = saturation_vapor_pressure.(T)
    es_T′ = saturation_vapor_pressure.(T′)
    RH = @. clamp(e / es_T, 0.0, 1.0)
    e′ = @. RH * es_T′

    # 3. Longwave down: recompute from adjusted T, e; preserve cloud increment Δε.
    LW′ = _adjust_longwave(LW, e, T, e′, T′)

    src_meta = metadata(climate_forcing_original)
    new_meta = merge(copy(src_meta), Dict(
        "delta_temperature" => ΔT,
        # Cumulative offset relative to the source reanalysis, so repeated
        # adjustments compose (mirrors "elevation_offset").
        "temperature_offset" => get(src_meta, "temperature_offset", 0.0) + ΔT,
        "temperature_air_mean" => Statistics.mean(T′),
    ))

    # Pressure, wind, precipitation and shortwave are carried through untouched.
    adjusted = _rebuild_forcing(climate_forcing_original, new_meta;
                                temperature_air=T′, vapor_pressure=e′,
                                longwave_downward=LW′)

    # Re-validate physical ranges after adjustment.
    validate_climate_forcing_units(adjusted)

    return adjusted
end

"""
    precipitation_adjust(climate_forcing_original::DimStack,
                         scaling::Real) -> DimStack

Rescale precipitation in a climate-forcing `DimStack` by a fractional factor,
leaving every other variable unchanged.

Use this for accumulation sensitivity experiments (`0.8` = 20 % drier, `1.15` =
15 % wetter) or to bias-correct a reanalysis precipitation rate against an
observed accumulation record. No other forcing variable depends on the
precipitation rate, so this is a single-variable change; the phase (snow vs rain)
is determined downstream from air temperature, so a scaling applies to whichever
phase the temperature implies.

# Arguments
- `climate_forcing_original::DimStack`: forcing as returned by
  [`climate_forcing`](@ref), with a `Ti` (time) dimension and the seven standard
  forcing variables.
- `scaling::Real`: dimensionless multiplicative factor applied to every time
  step. Must be non-negative; `1.0` is the identity, `0.0` removes all
  precipitation.

# Returns
- A new `DimStack` on the same time dimension with `precipitation` scaled and all
  other variables carried through unchanged. Metadata gains a cumulative
  `precipitation_scaling` and a recomputed `precipitation_mean` (annual,
  kg/m²/yr); physical-range validation is re-run.

# Method
`precip′ = scaling · precip`. The scaling is uniform in time — it changes the
total and the amplitude but not the timing or intermittency of events. Metadata
`precipitation_scaling` accumulates *multiplicatively*, so two successive calls
with `1.1` record `1.21`.

`scaling = 1.0` reproduces the input exactly.

!!! note
    A large `scaling` can push the hourly rate past the 100 kg/m²/hr validation
    ceiling, raising an `ArgumentError` from `validate_climate_forcing_units`.

# Example
```julia
stack = climate_forcing(:era5land, 72.58, -38.46;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

# 15 % more precipitation.
wetter = precipitation_adjust(stack, 1.15)

# Dry, warm scenario.
scenario = precipitation_adjust(temperature_adjust(stack, 3.0), 0.85)
```

See also [`temperature_adjust`](@ref), [`climate_adjust_for_elevation`](@ref).
"""
function precipitation_adjust(climate_forcing_original::DimStack, scaling::Real)
    s = Float64(scaling)
    if !(s >= 0)
        throw(ArgumentError(
            "scaling must be a non-negative fractional factor (got $(scaling)); " *
            "1.0 is the identity, 1.15 means 15% wetter"))
    end

    precip = climate_forcing_original[:precipitation]
    precip′ = @. precip * s

    src_meta = metadata(climate_forcing_original)
    new_meta = merge(copy(src_meta), Dict(
        # Cumulative and multiplicative, so repeated adjustments compose.
        "precipitation_scaling" => get(src_meta, "precipitation_scaling", 1.0) * s,
        # Hourly mean → annual total, matching the loader's convention.
        "precipitation_mean" => Statistics.mean(precip′) * 8760.0,
    ))

    # Precipitation is the only variable this touches.
    adjusted = _rebuild_forcing(climate_forcing_original, new_meta; precipitation=precip′)

    # Re-validate physical ranges after adjustment.
    validate_climate_forcing_units(adjusted)

    return adjusted
end
