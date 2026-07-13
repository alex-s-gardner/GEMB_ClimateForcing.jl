# Quick Start Guide

## Installation

```bash
cd ~/Documents/GitHub/GEMB_ClimateForcing.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Get a CDS API Key

1. Register at https://cds.climate.copernicus.eu/
2. Go to your profile page and copy your API key
3. Set environment variable:
   ```bash
   export CDS_API_KEY="your-api-key-here"
   ```

## Basic Usage

### Load ERA5-Land Forcing

```julia
using GEMB_ClimateForcing
using Dates

# Load 1 month of data for Summit, Greenland
cf = climate_forcing(
    :era5land, 72.58, -38.46;  # lat, lon
    time_range=(DateTime(2020,1,1), DateTime(2020,2,1)),
    token=ENV["CDS_API_KEY"]
)
```

### Run GEMB with ERA5-Land

```julia
using GEMB

# Initialize model
mp = ModelParameters(output_frequency=:daily)
profile = initialize_profile(mp, cf)

# Run simulation
output = gemb(profile, cf, mp)

# Analyze results
using Statistics
temp_surface = surface_timeseries(parent(output[:temperature]))
println("Mean surface temp: ", round(mean(temp_surface) - 273.15, digits=1), "°C")
```

## Run the Complete Example

```bash
cd ~/Documents/GitHub/GEMB_ClimateForcing.jl
export CDS_API_KEY="your-api-key-here"
julia --project=. examples/era5_land_example.jl
```

This will:
1. Load 1 year of ERA5-Land data for Summit, Greenland
2. Create climatological forcing for spinup
3. Spin up GEMB for 50 years
4. Run 1-year simulation with ERA5 forcing
5. Display results summary

## Run Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests include:
- Input validation
- Error handling
- Integration tests (if CDS_API_KEY is set)

## Performance Notes

- First load may take 30-60 seconds (depends on network)
- Geo-chunked strategy (default) is fastest for time-series
- Data is lazy-loaded - only downloads requested time/location
- Subsequent requests benefit from HTTP caching

## Troubleshooting

### "token is required for ERA5-Land access"
Set `CDS_API_KEY` environment variable before running.

### HTTP 401/403 errors
Check that your CDS API key is valid and not expired.

### Slow initial load
This is normal - xarray needs to load metadata from all 4 Zarr stores.

### Missing Python packages
PythonCall.jl automatically installs Python and xarray via CondaPkg.
First run may take longer while setting up the Python environment.
