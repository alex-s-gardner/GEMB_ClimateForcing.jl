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
    return 611.2 .* exp.(17.67 .* (T_dewpoint .- 273.15) ./ (T_dewpoint .- 29.65))
end
