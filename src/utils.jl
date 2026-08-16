"""
    get_cds_api_key()

Get the CDS (Copernicus Climate Data Store) API key from environment variable or .cdsapirc file.

# Returns
- String containing the CDS API key

# Description
Searches for the API key in the following order:
1. Environment variable `CDS_API_KEY`
2. `~/.cdsapirc` file (standard location for CDS API credentials)

# Throws
- `ErrorException` if no API key is found

# Example
```julia
token = get_cds_api_key()
forcing_data = climate_forcing(:era5land, lat, lon; time_range=..., token=token)
```

# See Also
- Get your CDS API key from: https://cds.climate.copernicus.eu/api-how-to
"""
function get_cds_api_key()
    # First check environment variable
    if haskey(ENV, "CDS_API_KEY")
        return ENV["CDS_API_KEY"]
    end

    # Try reading from .cdsapirc file
    cdsapirc = joinpath(homedir(), ".cdsapirc")
    if isfile(cdsapirc)
        for line in readlines(cdsapirc)
            if startswith(line, "key:")
                return String(strip(split(line, ":", limit=2)[2]))
            end
        end
    end

    error("CDS API key not found. Set ENV[\"CDS_API_KEY\"] or create ~/.cdsapirc file")
end

"""
    dewpoint_to_vapor_pressure(T_dewpoint::AbstractVector)

Convert dewpoint temperature (K) to vapor pressure (Pa) using the Magnus formula.

# Arguments
- `T_dewpoint`: Dewpoint temperature in Kelvin

# Returns
- Vapor pressure in Pascals (Pa)

# Formula
Uses the Magnus formula:
```
e = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
```

where T is in Kelvin and e is in Pascals.

# References
- Magnus, G. (1844). "Versuche über die Spannkräfte des Wasserdampfs"
- Commonly used approximation for meteorological applications
"""
function dewpoint_to_vapor_pressure(T_dewpoint)
    # Magnus formula: convert K to °C, compute vapor pressure
    # e = 611.2 * exp(17.67 * T_celsius / (T_celsius + 243.5))
    # Equivalent form: e = 611.2 * exp(17.67 * (T_K - 273.15) / (T_K - 29.65))
    #
    # Promote to Float64 inside the fused kernel so a native Float32 input can be passed
    # directly (no throwaway Float64 copy) while keeping full double precision in exp().
    return @. 611.2 * exp(17.67 * (Float64(T_dewpoint) - 273.15) / (Float64(T_dewpoint) - 29.65))
end

"""
    vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)

Calculate relative humidity [%] from vapor pressure [Pa] and air temperature [K].

Uses Tetens' formula for saturation vapor pressure and clamps the result to [0, 100].

Matches MATLAB's `vapor_pressure_to_relative_humidity.m`.
"""
function vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)
    Tc = temperature_air .- 273.15
    A = 610.78
    B = 17.27
    C = 237.3
    es = A .* exp.((B .* Tc) ./ (Tc .+ C))
    relative_humidity = (vapor_pressure ./ es) .* 100.0
    return clamp.(relative_humidity, 0.0, 100.0)
end

"""
    relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)

Estimate actual vapor pressure [Pa] from air temperature [K] and relative humidity [%].

Uses Tetens' formula for saturation vapor pressure.

Matches MATLAB's `relative_humidity_to_vapor_pressure.m`.
"""
function relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)
    Tc = temperature_air .- 273.15
    A = 610.78   # Pa
    B = 17.27    # dimensionless
    C = 237.3    # degrees Celsius
    es = A .* exp.((B .* Tc) ./ (Tc .+ C))
    return es .* (relative_humidity ./ 100.0)
end

# ----------------------------------------------------------------------------
# Shared helpers for the forcing-adjustment functions
# (elevation_adjustment.jl, climate_adjustment.jl, glacier_adjustment.jl).
# ----------------------------------------------------------------------------

"""
    _rebuild_forcing(stack, new_metadata; changed_layers...) -> DimStack

Rebuild a forcing `DimStack` with new stack metadata, replacing only the named layers.

Every adjustment function changes a different subset of the seven forcing variables; this
carries the rest through untouched rather than making each call site re-enumerate them.
`rebuild` preserves each layer's own metadata, so broadcast results (which drop it) come
back correctly annotated without threading it by hand.

Values may be `DimArray`s or plain arrays — only the underlying data is used.
"""
function _rebuild_forcing(stack::DimStack, new_metadata; changed_layers...)
    changed = map(parent, NamedTuple(changed_layers))
    unknown = setdiff(keys(changed), keys(stack))
    isempty(unknown) ||
        throw(ArgumentError("not a layer of this forcing stack: $(join(unknown, ", "))"))
    return rebuild(stack; data=merge(NamedTuple(stack), changed), metadata=new_metadata)
end

"""
    _adjust_longwave(LW, e, T, e′, T′)

Recompute downwelling longwave irradiance [W/m²] for adjusted temperature `T′` and vapor
pressure `e′`, preserving the cloud/aerosol emissivity increment diagnosed from the input.

`Δε = LW/(σT⁴) − ε_cs(e, T)` is the departure of the observed irradiance from the
Konzelmann et al. (1994) clear-sky value; holding it fixed means only the clear-sky
response and the `σT⁴` scaling change. Exact identity when `(e′, T′) == (e, T)`.
"""
function _adjust_longwave(LW, e, T, e′, T′)
    ε_cs  = konzelmann_clear_sky_emissivity.(e, T)
    ε_cs′ = konzelmann_clear_sky_emissivity.(e′, T′)
    Δε = @. LW / (_SIGMA_SB * T^4) - ε_cs
    return @. (ε_cs′ + Δε) * _SIGMA_SB * T′^4
end

# Normalize longitude to the −180…180°E convention (accepts 0–360°E, or any multiple
# wrap: 540 → 180, −190 → 170). Used by the geoid grid and the glacier lookup table,
# both of which are −180…180 while the ERA5-Land grids are 0–360.
_wrap_longitude(lon::Real) = mod(lon + 180.0, 360.0) - 180.0
