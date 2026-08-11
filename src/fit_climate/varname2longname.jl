"""
    varname2longname(varname::String)

Return a descriptive long name for a given variable name.

# Arguments
- `varname::String`: Short variable name (e.g., "temperature_air")

# Returns
- `String`: Long descriptive name with units (e.g., "screen level air temperature [K]")

# Details
Maps standard GEMB variable names to human-readable descriptions with units.
Used for plotting labels and documentation.

Supported variables:
- `"longwave_downward"`: downward longwave radiation [W m⁻²]
- `"shortwave_downward"`: downward shortwave radiation [W m⁻²]
- `"relative_humidity"`: screen level relative humidity [%]
- `"temperature_air"`: screen level air temperature [K]
- `"precipitation"`: precipitation [kg m⁻²]
- `"wind_speed"`: screen level wind speed [m s⁻¹]
- `"vapor_pressure"`: vapor pressure [Pa]
- `"pressure_air"`: screen level air pressure [Pa]

Matches MATLAB's `varname2longname.m`.

# References
Original MATLAB implementation in GEMB/src/fit_simulated_climate_to_data/

# Example
```julia
julia> varname2longname("temperature_air")
"screen level air temperature [K]"
```
"""
function varname2longname(varname::String)
    # Define unicode superscripts for units
    sup_minus = "⁻"
    sup_1 = "¹"
    sup_2 = "²"

    # Construct unit strings
    unit_wm2 = " [W m" * sup_minus * sup_2 * "]"  # [W m⁻²]
    unit_kgm2 = " [kg m" * sup_minus * sup_2 * "]" # [kg m⁻²]
    unit_ms1 = " [m s" * sup_minus * sup_1 * "]"  # [m s⁻¹]

    # Variable mapping
    var_defs = Dict(
        "longwave_downward" => "downward longwave radiation" * unit_wm2,
        "shortwave_downward" => "downward shortwave radiation" * unit_wm2,
        "relative_humidity" => "screen level relative humidity [%]",
        "temperature_air" => "screen level air temperature [K]",
        "precipitation" => "precipitation" * unit_kgm2,
        "wind_speed" => "screen level wind speed" * unit_ms1,
        "vapor_pressure" => "vapor pressure [Pa]",
        "pressure_air" => "screen level air pressure [Pa]"
    )

    if haskey(var_defs, varname)
        return var_defs[varname]
    else
        error("Unknown variable name: $varname")
    end
end
