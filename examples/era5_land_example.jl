# Example: Load ERA5-Land forcing data and run GEMB simulation
#
# This example demonstrates pure Julia ERA5-Land data access:
# 1. Authenticated HTTPS access to ERA5-Land ARCO Zarr stores (no Python!)
# 2. Loading climate forcing using Zarr.jl + AuthenticatedHTTPStore
# 3. Running GEMB with ERA5 forcing
# 4. Extracting and visualizing results
#
# Requirements:
# - CDS API key from https://cds.climate.copernicus.eu/
# - Set environment variable: export CDS_API_KEY="your-token-here"
# - Pure Julia - no PythonCall or xarray needed!

using GEMB
using GEMB_ClimateForcing
using Dates
using Statistics

## Configuration
# Summit Station, Greenland (ice sheet summit)
lat = 72.58
lon = -38.46

# Time range: 1 year of data
time_range = (DateTime(2020, 1, 1), DateTime(2020, 12, 31))

# Get CDS API key from environment
token = get(ENV, "CDS_API_KEY", nothing)
if isnothing(token)
    error("CDS_API_KEY environment variable not set. Get a key from https://cds.climate.copernicus.eu/")
end

## Load ERA5-Land forcing data
println("="^60)
println("Loading ERA5-Land Climate Forcing Data")
println("="^60)

forcing_data = climate_forcing(
    :era5land, lat, lon;
    time_range=time_range,
    token=token,
    chunk_strategy=:geo  # Optimized for time-series extraction
)

# Convert DimStack to GEMB.ClimateForcing (via extension)
cf = GEMB.ClimateForcing(forcing_data)

println("\nForcing data summary:")
println("  Time steps: ", length(cf.temperature_air))
println("  Mean temperature: ", round(cf.temperature_air_mean - 273.15, digits=1), "°C")
println("  Mean wind speed: ", round(cf.wind_speed_mean, digits=1), " m/s")
println("  Mean precipitation: ", round(cf.precipitation_mean, digits=0), " kg/m²/yr")

## Initialize GEMB

println("\n", "="^60)
println("Initializing GEMB")
println("="^60)

# Model parameters (daily output)
mp = GEMB.ModelParameters(output_frequency=:daily)

# Initialize vertical profile
profile = GEMB.initialize_profile(mp, cf)

## Create climatological forcing for spinup

println("\nCreating climatological forcing for spinup...")
cf_spinup = GEMB.forcing_climatology(cf)

# Spin up for 50 years to reach equilibrium
println("Running 50-year spinup...")
mp_spinup = GEMB.ModelParameters(output_frequency=:last)
profile_spunup = GEMB.gemb_spinup(profile, cf_spinup, mp_spinup, 50)

println("  Spinup complete!")
println("  Final surface density: ", round(profile_spunup.density[1], digits=1), " kg/m³")
println("  Final surface temperature: ", round(profile_spunup.temperature[1] - 273.15, digits=1), "°C")

## Run GEMB with ERA5 forcing

println("\n", "="^60)
println("Running GEMB Simulation")
println("="^60)

output = GEMB.gemb(profile_spunup, cf, mp)

println("\nSimulation complete!")
println("  Output time steps: ", size(output[:temperature], 2))
println("  Vertical layers: ", size(output[:temperature], 1))

## Analyze results

println("\n", "="^60)
println("Results Summary")
println("="^60)

# Surface time series
temp_surface = GEMB.surface_timeseries(parent(output[:temperature]))
albedo_surface = parent(output[:albedo_surface])
density_surface = GEMB.surface_timeseries(parent(output[:density]))

println("\nSurface conditions:")
println("  Mean temperature: ", round(mean(temp_surface) - 273.15, digits=1), "°C")
println("  Temperature range: ", round(minimum(temp_surface) - 273.15, digits=1), " to ",
        round(maximum(temp_surface) - 273.15, digits=1), "°C")
println("  Mean albedo: ", round(mean(albedo_surface), digits=3))
println("  Mean surface density: ", round(mean(density_surface), digits=1), " kg/m³")
println("  Firn air content: ", round(mean(parent(output[:firn_air_content])), digits=2), " m")

# Annual totals
total_precip = sum(cf.precipitation)
total_melt = sum(parent(output[:melt]))
total_runoff = sum(parent(output[:runoff]))
total_refreezing = sum(parent(output[:refreezing]))

println("\nAnnual mass balance:")
println("  Total precipitation: ", round(total_precip, digits=0), " kg/m²")
println("  Total melt: ", round(total_melt, digits=0), " kg/m²")
println("  Total runoff: ", round(total_runoff, digits=0), " kg/m²")
println("  Total refreezing: ", round(total_refreezing, digits=0), " kg/m²")
println("  Net accumulation: ", round(total_precip - total_runoff, digits=0), " kg/m²")

println("\n", "="^60)
println("Example Complete!")
println("="^60)
