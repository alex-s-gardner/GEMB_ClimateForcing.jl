"""
    simulate_coeffs_disp(coeffs::NamedTuple, struct_name::String="coeffs")

Display coefficient structure fields as executable Julia code.

# Arguments
- `coeffs::NamedTuple`: Structure containing coefficient values (scalars or vectors)
- `struct_name::String="coeffs"`: Name to use in the output code

# Details
Iterates through the fields of a coefficient NamedTuple and prints them to stdout.
The output is formatted as valid Julia assignment statements, making it easy to
copy, paste, and reuse parameter sets for debugging or configuration files.

Matches MATLAB's `simulate_coeffs_disp.m`.

# References
Original MATLAB implementation in GEMB/src/fit_simulated_climate_to_data/

# Example
```julia
# Fit temperature parameters
coeffs_temp = fit_air_temperature(dec_year, temp, lat, elev)

# Display as copy-pasteable code
simulate_coeffs_disp(coeffs_temp, "coeffs.temperature_air")

# Output (example):
# coeffs.temperature_air.mean_offset = 2.3456
# coeffs.temperature_air.lat_scale = 1.0234
# ...
```
"""
function simulate_coeffs_disp(coeffs::NamedTuple, struct_name::String="coeffs")
    keys_list = keys(coeffs)

    for key in keys_list
        value = getfield(coeffs, key)

        if value isa AbstractVector
            # Handle vectors
            if length(value) == 1
                # Single element vector
                println("$struct_name.$key = [$(round(value[1], digits=4))]")
            else
                # Multi-element vector
                vals_str = join([string(round(v, digits=4)) for v in value], " ")

                # Print as Julia vector
                println("$struct_name.$key = [$vals_str]")
            end
        else
            # Handle scalars
            println("$struct_name.$key = $(round(value, digits=4))")
        end
    end

    println("")  # Blank line for readability
end
