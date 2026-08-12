"""
GEMB_ClimateForcing

Load climate forcing data from various reanalysis datasets and return DimStack
with climate variables. Supports ERA5-Land, ERA5, MERRA-2, and other datasets.

# Example
```julia
using GEMB_ClimateForcing
using GEMB  # Extension provides DimStack → ClimateForcing conversion

# Load ERA5-Land data for Summit, Greenland
forcing_data = climate_forcing(
    :era5land, 72.58, -38.46;
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"]
)

# Convert to GEMB.ClimateForcing (requires GEMB.jl)
cf = GEMB.ClimateForcing(forcing_data)

# Use with GEMB
output = gemb(profile, cf, mp)
```
"""
module GEMB_ClimateForcing

using Dates
using DimensionalData
using Random
using Statistics
using Zarr
using HTTP
using OpenSSL

# Export main interface
export climate_forcing, climate_chunk_map, get_cds_api_key, climate_adjust_for_elevation,
    empirical_lapse_rate, climate_model_invariant, ERA5_LAND_INVARIANT_PARAMETERS,
    GREENLAND_LAPSE_RATE, ARCTIC_LAPSE_RATE, ANTARCTICA_LAPSE_RATE

# Export synthetic forcing and climate fitting
export simulate_climate_forcing, simulation_parameter_sets,
    datetime2decyear, decyear2datetime,
    dewpoint_to_vapor_pressure, vapor_pressure_to_relative_humidity,
    relative_humidity_to_vapor_pressure,
    fit_air_temperature, fit_precipitation, fit_longwave_irradiance_delta,
    fit_seasonal_daily_noise, varname2longname, simulate_coeffs_disp,
    geopotential2height, geoid_undulation, GEOID_MODELS

# Include submodules
include("utils.jl")
include("authenticated_http_store.jl")
include("interface.jl")
include("datasets/era5_land.jl")
include("datasets/copernicus_dem.jl")
include("geoid.jl")
include("invariant.jl")
include("elevation_adjustment.jl")

# Synthetic forcing generation
include("simulate/simulate_climate_forcing.jl")

# Climate fitting functions
include("fit_climate/fit_air_temperature.jl")
include("fit_climate/fit_precipitation.jl")
include("fit_climate/fit_longwave_irradiance_delta.jl")
include("fit_climate/fit_seasonal_daily_noise.jl")
include("fit_climate/varname2longname.jl")
include("fit_climate/simulate_coeffs_disp.jl")

end # module
