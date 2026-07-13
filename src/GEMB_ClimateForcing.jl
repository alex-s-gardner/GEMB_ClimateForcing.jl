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
using Statistics
using Zarr
using HTTP
using OpenSSL

# Export main interface
export climate_forcing

# Include submodules
include("utils.jl")
include("authenticated_http_store.jl")
include("interface.jl")
include("datasets/era5_land.jl")

end # module
