# Off-glacier (ambient) → on-glacier air temperature correction.
#
# A melting glacier surface pinned at 0 °C cools the air above it, and the resulting stable
# layer suppresses turbulent mixing with the free atmosphere and drives a katabatic glacier
# wind. Reanalysis 2 m temperature carries none of this: ERA5-Land grid cells dwarf a valley
# glacier and its land-surface scheme has no glacier boundary layer. Feeding ambient
# temperature straight to a surface energy balance model therefore overestimates melt — 22 %
# of the mass-balance change per +1 °C (Greuell & Böhm 1998), up to 32 % of seasonal melt
# (Troxler et al. 2020).
#
# This applies the Shaw et al. (2025) decoupling factor `k` from the published per-glacier
# lookup table (see glacier_decoupling.jl). The paper's five-predictor regression is
# deliberately *not* re-implemented: its Table S2 has sign errors in a4/a5, and even a
# corrected refit only reaches R² = 0.54 against the authors' own published `k`.
#
# Full derivation, coefficient tables and the alternatives considered (SM10, ModGB):
# docs/on_glacier_temperature_correction.md

using DimensionalData
using Statistics

# `_rebuild_forcing` / `_cap_vapor_pressure` / `_wrap_longitude` (utils.jl) and
# `validate_climate_forcing_units` (datasets/era5_land.jl) are shared with the sibling
# adjustment functions and reused here.

"""
    climate_adjust_for_glacier(stack::DimStack, decoupling_factor::Real; kwargs...) -> DimStack
    climate_adjust_for_glacier(stack::DimStack; rgi_id=nothing, kwargs...) -> DimStack

Correct an *ambient* (off-glacier) climate forcing to *on-glacier* conditions, accounting for
the glacier boundary layer that cools the air over a melting surface.

The first form takes the decoupling factor `k` directly. The second looks it up in the
Shaw et al. (2025) per-glacier table — by RGI v6 identifier if `rgi_id` is given, otherwise
by nearest glacier centroid to the `"latitude"`/`"longitude"` in the stack metadata.

!!! warning "Apply the elevation correction first"
    `k` multiplies an ambient temperature *at the glacier's elevation*. Run
    [`climate_adjust_for_elevation`](@ref) first to move the reanalysis forcing to the target
    elevation, then this second. The reverse order is wrong.

# Arguments
- `stack::DimStack`: forcing as returned by [`climate_forcing`](@ref), with a `Ti` (time)
  dimension and the seven standard forcing variables.
- `decoupling_factor::Real`: unitless `k` in `(0, 1]`. `1.0` is the identity (fully coupled,
  no glacier boundary layer); Shaw et al. clamp their estimates to `[0.2, 1.0]` with a global
  mean of ≈0.83.

# Keywords
- `rgi_id::Union{Nothing,AbstractString}=nothing`: RGI v6 identifier (e.g. `"RGI60-11.02810"`)
  for the table lookup. Ignored when `k` is passed positionally.
- `max_distance::Real=10.0`: km, bound on the nearest-centroid lookup when `rgi_id` is not
  given. See [`glacier_decoupling`](@ref).
- `reference_temperature::Real=273.15`: K, the temperature about which `k` pivots. The melting
  point is the physically meaningful choice and matches the °C basis of the published
  regression; exposed mainly for sensitivity testing.
- `apply_below_freezing::Bool=false`: by default the correction is applied only where ambient
  temperature exceeds `reference_temperature`, leaving colder steps untouched. See the note
  below — do not set this to `true` casually.

# Returns
- A new `DimStack` on the same time dimension with `temperature_air` cooled and
  `vapor_pressure` capped at saturation of the cooled air. `longwave_downward`,
  `pressure_air`, `wind_speed`, `precipitation`, `shortwave_downward` and any further layer
  are carried through **unchanged**. Metadata gains a
  cumulative `glacier_decoupling_factor` (multiplicative, so repeated calls compose), the
  provenance of the lookup when one was used, and a recomputed `temperature_air_mean`;
  physical-range validation is re-run.

# Method
In °C the Shaw et al. relation is `T_gla = k·T_amb`, i.e. cooling proportional to how far
ambient temperature sits above melting. Written to be an exact identity at `k = 1`:

```
T′ = T + (k − 1)·(T − T_ref)                      apply_below_freezing = true
T′ = T + (k − 1)·max(T − T_ref, 0)                apply_below_freezing = false (default)
```

**Downwelling longwave is left unchanged.** The decoupling is a shallow surface layer, a few
meters to tens of meters thick, while downwelling longwave is emitted by the whole lower
atmosphere and cloud base above it; cooling that layer changes little of the emitting column.
What `k` measures is the 2 m temperature response to the reanalysis glacier mask, and ERA5-Land's
downwelling longwave does not depend on that mask, so the data carry no signal for adjusting it.

**Vapor pressure is left untouched except for one physical cap**, and this is *not* the
constant-relative-humidity propagation used by [`temperature_adjust`](@ref):

```
e′ = min(e, eₛ(T′))   where T′ < T;   e′ = e elsewhere
```

with `eₛ` the saturation vapor pressure over ice below 273.16 K and over water above
([`saturation_vapor_pressure_reference`](@ref)). Cooled air
cannot hold more vapor than saturation at its new temperature; the excess is taken to condense.
The cap is the only clamp, and it applies only where the correction cooled the air, so `k = 1`
and every sub-freezing step under the default gating return the input bit for bit.

Beyond the cap, `e` is not propagated with temperature: Shea & Moore (2010, Eq. 4) show `e_gla` is a linear function of `e_amb` that pivots about 6.11 hPa — saturation over ice at
0 °C — so the boundary layer *adds* moisture whenever ambient air is drier than that, the
opposite sign to constant-RH scaling. Their sub-freezing branch is near-identity. Applying the
right correction needs the flowpath-length-dependent `j1`/`j2` coefficients, hence a flowpath
raster; without one, leaving `e` alone is the defensible choice. Wind speed is likewise
untouched: no published scheme corrects it, and Shaw et al. (2024) show the katabatic jet can
*enhance* turbulent heat exchange even while the air is cooler, partly offsetting the cooling.

!!! note "Melt-season gating"
    The underlying regression is calibrated on ablation-season data. A bare `k` multiplier
    applied to a sub-freezing ambient temperature *warms* it (`k·(−10) > −10`), which has no
    physical basis — the boundary layer that produces the decoupling requires a melting
    surface. Hence `apply_below_freezing=false` by default, which reduces the scheme to a
    threshold form (identity below melting) in the spirit of Shea & Moore (2010).

!!! note "Coverage"
    RGI regions 05 (Greenland periphery) and 19 (Antarctic and Subantarctic) are absent from
    the Shaw et al. table, so the lookup forms will fail there. Pass `k` explicitly — a
    regional mean, or a value from another scheme — for those regions.

# Example
```julia
stack = climate_forcing(:era5land, 45.97, 7.53;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

# Elevation correction first, to the glacier surface, then the glacier correction.
z_glacier = 2900.0
at_glacier = climate_adjust_for_elevation(stack, z_glacier - metadata(stack)["elevation"])

# Look up k for Haut Glacier d'Arolla by RGI id...
on_glacier = climate_adjust_for_glacier(at_glacier; rgi_id="RGI60-11.02810")

# ...or by nearest centroid to the forcing location, implicit from the metadata.
on_glacier = climate_adjust_for_glacier(at_glacier)

# ...or with k supplied directly (e.g. the lower CI bound, for a sensitivity test).
row = glacier_decoupling("RGI60-11.02810")
low = climate_adjust_for_glacier(at_glacier, row.k_lower)
```

# References
- Shaw, T. E. et al. (2025). Decoupling of air temperature over mountain glaciers.
  *Nature Climate Change*. Dataset: https://doi.org/10.5281/zenodo.14044846
- Shea, J. M. & Moore, R. D. (2010). Prediction of glacier boundary layer temperature and
  vapor pressure. *J. Geophys. Res.* 115, D23107.
- Greuell, W. & Böhm, R. (1998). 2 m temperatures along melting mid-latitude glaciers.
  *J. Glaciol.* 44, 9–20.
- Konzelmann, T. et al. (1994). *Glob. Planet. Change* 9, 143–164. (Longwave emissivity.)

See also [`glacier_decoupling`](@ref), [`climate_adjust_for_elevation`](@ref),
[`temperature_adjust`](@ref).
"""
function climate_adjust_for_glacier(stack::DimStack, decoupling_factor::Real;
                                    reference_temperature::Real=273.15,
                                    apply_below_freezing::Bool=false)
    k = Float64(decoupling_factor)
    if !(0 < k <= 1)
        throw(ArgumentError(
            "decoupling_factor must be in (0, 1] (got $(decoupling_factor)); 1.0 is the " *
            "identity and Shaw et al. (2025) clamp their estimates to [0.2, 1.0]. A value " *
            "above 1 would warm the glacier surface layer, which no published scheme supports."))
    end
    T_ref = Float64(reference_temperature)

    T = stack[:temperature_air]
    e = stack[:vapor_pressure]

    # 1. Temperature. Written as an increment so k = 1 is bit-exact.
    T′ = apply_below_freezing ? (@. T + (k - 1) * (T - T_ref)) :
                                (@. T + (k - 1) * max(T - T_ref, 0.0))

    # 2. Vapor pressure: capped at saturation of the cooled air, and nowhere else.
    e′ = _cap_vapor_pressure.(e, T, T′)

    src_meta = metadata(stack)
    new_meta = merge(copy(src_meta), Dict(
        "decoupling_factor" => k,
        # Cumulative and multiplicative: above melting two successive factors compose as
        # their product, mirroring "precipitation_scaling".
        "glacier_decoupling_factor" => get(src_meta, "glacier_decoupling_factor", 1.0) * k,
        "glacier_decoupling_reference_temperature" => T_ref,
        "glacier_decoupling_below_freezing" => apply_below_freezing,
        "temperature_air_mean" => Statistics.mean(T′),
    ))
    # Longwave, pressure, wind, precipitation and shortwave are carried through untouched — see
    # the docstring for why longwave and vapor pressure are not re-derived.
    adjusted = _rebuild_forcing(stack, new_meta; temperature_air=T′, vapor_pressure=e′)

    validate_climate_forcing_units(adjusted)

    return adjusted
end

function climate_adjust_for_glacier(stack::DimStack;
                                    rgi_id::Union{Nothing,AbstractString}=nothing,
                                    max_distance::Real=10.0, kwargs...)
    row = if rgi_id !== nothing
        glacier_decoupling(rgi_id)
    else
        src_meta = metadata(stack)
        lat = get(src_meta, "latitude", nothing)
        lon = get(src_meta, "longitude", nothing)
        if lat === nothing || lon === nothing
            throw(ArgumentError(
                "cannot look up a decoupling factor: no rgi_id given and the stack " *
                "metadata has no \"latitude\"/\"longitude\". Pass rgi_id=..., or supply " *
                "the factor directly as climate_adjust_for_glacier(stack, k)."))
        end
        glacier_decoupling(lat, lon; max_distance=max_distance)
    end

    @info "Applying Shaw et al. (2025) glacier decoupling" rgi_id=row.rgi_id k=row.k match_distance_km=round(row.distance, digits=2)

    adjusted = climate_adjust_for_glacier(stack, row.k; kwargs...)

    # Record where k came from. Pure metadata, so it is attached afterwards rather than
    # threaded through the core method as a private keyword.
    return rebuild(adjusted; metadata=merge(metadata(adjusted), Dict(
        "glacier_decoupling_source" => "Shaw et al. (2025) lookup table",
        "glacier_decoupling_rgi_id" => row.rgi_id,
        "glacier_decoupling_k_lower" => row.k_lower,
        "glacier_decoupling_match_distance" => row.distance,
    )))
end
