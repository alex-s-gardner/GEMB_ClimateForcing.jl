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
function dewpoint_to_vapor_pressure(T_dewpoint::AbstractVector{<:Real})
    # Magnus formula: convert K to °C, compute vapor pressure
    # e = 611.2 * exp(17.67 * T_celsius / (T_celsius + 243.5))
    # Equivalent form: e = 611.2 * exp(17.67 * (T_K - 273.15) / (T_K - 29.65))
    #
    # Promote to Float64 inside the fused kernel so a native Float32 input can be passed
    # directly (no throwaway Float64 copy) while keeping full double precision in exp().
    return @. 611.2 * exp(17.67 * (Float64(T_dewpoint) - 273.15) / (Float64(T_dewpoint) - 29.65))
end
