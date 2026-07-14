"""
ERA5-Land ARCO (Analysis-Ready, Cloud-Optimized) Zarr data loader.

Pure Julia implementation using Zarr.jl with authenticated HTTP access.
Loads data from ECMWF's cloud-optimized Zarr stores via HTTPS with Bearer token.
"""

using Zarr
using DimensionalData
using Dates
using Statistics

# Import task-based parallelism primitives
import Base.Threads: @spawn

# ERA5-Land ARCO Zarr store base URL
const ERA5_LAND_BASE = "https://arco.datastores.ecmwf.int"

# Mapping of variable groups to ARCO store IDs
const ERA5_LAND_STORES = Dict(
    "sfc-2m-temperature" => "007",
    "sfc-pressure-precipitation" => "009",
    "sfc-wind" => "008",
    "sfc-radiation-heat" => "010"
)

"""
    era5_land_url(variable_group::String, chunk_strategy::Symbol) -> String

Construct ERA5-Land ARCO Zarr store URL.

# Arguments
- `variable_group`: Variable group name (e.g., "sfc-2m-temperature")
- `chunk_strategy`: :geo or :time

# Returns
- URL string for the Zarr store
"""
function era5_land_url(variable_group::String, chunk_strategy::Symbol)
    store_id = ERA5_LAND_STORES[variable_group]
    chunk_str = chunk_strategy == :geo ? "geo" : "time"
    return "$(ERA5_LAND_BASE)/cadl-arco-$(chunk_str)-$(store_id)/arco/reanalysis_era5_land/$(variable_group)/$(chunk_str)Chunked.zarr"
end

"""
    find_nearest_index(values::AbstractVector, target) -> Int

Find index of nearest value in sorted or unsorted array.
"""
function find_nearest_index(values::AbstractVector, target)
    _, idx = findmin(abs.(values .- target))
    return idx
end

"""
    validate_climate_forcing_units(stack::DimStack)

Validate that climate forcing variables in a DimStack have physically reasonable values.

# Expected Units and Ranges
- `temperature_air`: Kelvin (K), range [180, 330]
- `pressure_air`: Pascal (Pa), range [30000, 110000]
- `precipitation`: kg/m², range [0, 100] per hour
- `wind_speed`: m/s, range [0, 100]
- `shortwave_downward`: W/m², range [0, 1500]
- `longwave_downward`: W/m², range [50, 500]
- `vapor_pressure`: Pascal (Pa), range [0, 10000]

# Arguments
- `stack::DimStack`: DimStack with climate forcing variables

# Throws
- `ArgumentError` if any variable has values outside expected physical ranges
"""
function validate_climate_forcing_units(stack::DimStack)
    # Extract data vectors from DimStack
    temperature_air = parent(stack[:temperature_air])
    pressure_air = parent(stack[:pressure_air])
    precipitation = parent(stack[:precipitation])
    wind_speed = parent(stack[:wind_speed])
    shortwave_downward = parent(stack[:shortwave_downward])
    longwave_downward = parent(stack[:longwave_downward])
    vapor_pressure = parent(stack[:vapor_pressure])
    errors = String[]

    # Temperature (K): should be in reasonable range for Earth's surface
    t_min, t_max = extrema(temperature_air)
    if t_min < 180.0 || t_max > 330.0
        push!(errors, "temperature_air: expected [180, 330] K, got [$(t_min), $(t_max)] K")
    end

    # Pressure (Pa): surface pressure range
    p_min, p_max = extrema(pressure_air)
    if p_min < 30000.0 || p_max > 110000.0
        push!(errors, "pressure_air: expected [30000, 110000] Pa, got [$(p_min), $(p_max)] Pa")
    end

    # Precipitation (kg/m² per hour): should be non-negative and not extreme
    pr_min, pr_max = extrema(precipitation)
    if pr_min < -1e-6  # Allow small numerical errors
        push!(errors, "precipitation: expected ≥ 0 kg/m², got minimum $(pr_min) kg/m²")
    end
    if pr_max > 100.0
        push!(errors, "precipitation: expected ≤ 100 kg/m²/hr, got maximum $(pr_max) kg/m²/hr (possible unit error: should be kg/m², not m)")
    end

    # Wind speed (m/s): should be non-negative
    ws_min, ws_max = extrema(wind_speed)
    if ws_min < -1e-6
        push!(errors, "wind_speed: expected ≥ 0 m/s, got minimum $(ws_min) m/s")
    end
    if ws_max > 100.0
        push!(errors, "wind_speed: expected ≤ 100 m/s, got maximum $(ws_max) m/s (hurricane-force winds)")
    end

    # Shortwave radiation (W/m²): should be non-negative and below solar constant
    sw_min, sw_max = extrema(shortwave_downward)
    if sw_min < -1e-6
        push!(errors, "shortwave_downward: expected ≥ 0 W/m², got minimum $(sw_min) W/m²")
    end
    if sw_max > 1500.0
        push!(errors, "shortwave_downward: expected ≤ 1500 W/m², got $(sw_max) W/m² (possible unit error: should be W/m², not J/m²)")
    end

    # Longwave radiation (W/m²): should be in thermal radiation range
    lw_min, lw_max = extrema(longwave_downward)
    if lw_min < 50.0
        push!(errors, "longwave_downward: expected ≥ 50 W/m², got $(lw_min) W/m² (possible unit error: should be W/m², not J/m²)")
    end
    if lw_max > 500.0
        push!(errors, "longwave_downward: expected ≤ 500 W/m², got $(lw_max) W/m²")
    end

    # Vapor pressure (Pa): should be non-negative and below saturation
    vp_min, vp_max = extrema(vapor_pressure)
    if vp_min < -1e-6
        push!(errors, "vapor_pressure: expected ≥ 0 Pa, got minimum $(vp_min) Pa")
    end
    if vp_max > 10000.0
        push!(errors, "vapor_pressure: expected ≤ 10000 Pa, got $(vp_max) Pa (exceeds saturation at typical temperatures)")
    end

    # Throw error if any validation failed
    if !isempty(errors)
        error_msg = "Climate forcing unit validation failed:\n  " * join(errors, "\n  ")
        throw(ArgumentError(error_msg))
    end

    return true
end

"""
    load_era5_land(lat::Real, lon::Real; kwargs...) -> DimStack

Load ERA5-Land reanalysis data using pure Julia (Zarr.jl) and return a DimStack with climate variables.

# Keyword Arguments
- `time_range::Tuple{DateTime,DateTime}`: Time range to extract (required)
- `token::Union{String,Nothing}`: CDS API key (required)
- `chunk_strategy::Symbol=:geo`: :geo (time-series) or :time (spatial)
- `cache_path::Union{String,Nothing}=nothing`: Path for persistent disk cache (uses Zarr.CachingStore)

# Returns
- `DimStack`: Stack with climate forcing variables as DimArrays:
  - `temperature_air`: Air temperature (K)
  - `pressure_air`: Surface pressure (Pa)
  - `vapor_pressure`: Vapor pressure (Pa)
  - `wind_speed`: Wind speed magnitude (m/s)
  - `precipitation`: Precipitation rate (kg/m²/hr)
  - `shortwave_downward`: Downward shortwave radiation (W/m²)
  - `longwave_downward`: Downward longwave radiation (W/m²)

  Metadata includes `latitude`, `longitude`, `dataset`, `chunk_strategy`

# Notes
- ERA5-Land is on 0.1° grid (approximately 9 km)
- Uses nearest neighbor selection for lat/lon
- Hourly temporal resolution
- Requires free CDS API key from https://cds.climate.copernicus.eu/
- Pure Julia implementation using Zarr.jl + DimensionalData.jl

# Performance
- Uses task-based parallelism for concurrent variable group loading
- First load: ~10-25 seconds (1 year, network dependent)
- Subsequent loads: faster due to HTTP caching
- Parallel loading provides ~1.5-2x speedup over sequential loading
"""
function load_era5_land(
    lat::Real,
    lon::Real;
    time_range::Tuple{DateTime,DateTime},
    token::Union{String,Nothing}=nothing,
    chunk_strategy::Symbol=:geo,
    cache_path::Union{String,Nothing}=nothing,
    kwargs...
)
    # Validate token
    if isnothing(token)
        throw(ArgumentError("token is required for ERA5-Land access. Get one from https://cds.climate.copernicus.eu/"))
    end

    start_time, end_time = time_range

    println("Loading ERA5-Land data from cloud-optimized Zarr stores (pure Julia)...")
    println("  Location: $(lat)°N, $(lon)°E")
    println("  Time range: $(start_time) to $(end_time)")
    println("  Chunk strategy: $(chunk_strategy)")

    try
        # Create authenticated stores for each variable group
        println("  Creating authenticated HTTP stores...")

        url_temp = era5_land_url("sfc-2m-temperature", chunk_strategy)
        url_precip = era5_land_url("sfc-pressure-precipitation", chunk_strategy)
        url_wind = era5_land_url("sfc-wind", chunk_strategy)
        url_rad = era5_land_url("sfc-radiation-heat", chunk_strategy)

        # Create authenticated HTTP stores
        store_temp = AuthenticatedHTTPStore(url_temp; token=token)
        store_precip = AuthenticatedHTTPStore(url_precip; token=token)
        store_wind = AuthenticatedHTTPStore(url_wind; token=token)
        store_rad = AuthenticatedHTTPStore(url_rad; token=token)

        # Wrap with Zarr.CachingStore for persistent disk caching if cache_path provided
        if !isnothing(cache_path)
            mkpath(cache_path)
            store_temp = Zarr.CachingStore(store_temp, Zarr.DirectoryStore(joinpath(cache_path, "sfc-2m-temperature")))
            store_precip = Zarr.CachingStore(store_precip, Zarr.DirectoryStore(joinpath(cache_path, "sfc-pressure-precipitation")))
            store_wind = Zarr.CachingStore(store_wind, Zarr.DirectoryStore(joinpath(cache_path, "sfc-wind")))
            store_rad = Zarr.CachingStore(store_rad, Zarr.DirectoryStore(joinpath(cache_path, "sfc-radiation-heat")))
            println("  Using disk cache: $(cache_path)")
        end

        # Open Zarr groups in parallel (highest impact for performance)
        println("  Opening Zarr groups in parallel...")
        zarr_groups = @sync begin
            t1 = @spawn Zarr.zopen(store_temp, consolidated=true, fill_as_missing=false)
            t2 = @spawn Zarr.zopen(store_precip, consolidated=true, fill_as_missing=false)
            t3 = @spawn Zarr.zopen(store_wind, consolidated=true, fill_as_missing=false)
            t4 = @spawn Zarr.zopen(store_rad, consolidated=true, fill_as_missing=false)
            (fetch(t1), fetch(t2), fetch(t3), fetch(t4))
        end
        zg_temp, zg_precip, zg_wind, zg_rad = zarr_groups

        # Access individual arrays
        println("  Accessing variables...")
        t2m_zarr = zg_temp["t2m"]
        d2m_zarr = zg_temp["d2m"]
        sp_zarr = zg_precip["sp"]
        tp_zarr = zg_precip["tp"]
        u10_zarr = zg_wind["u10"]
        v10_zarr = zg_wind["v10"]
        ssrd_zarr = zg_rad["ssrd"]
        strd_zarr = zg_rad["strd"]

        # Get coordinate arrays
        println("  Reading coordinates...")
        # ERA5-Land coordinates: time, latitude, longitude
        # Dimensions are (time, latitude, longitude)

        lat_values = zg_temp["latitude"][:]
        lon_values = zg_temp["longitude"][:]
        time_values_raw = zg_temp["time"][:]

        # Convert time from hours since epoch to DateTime
        # ERA5-Land uses "hours since 1900-01-01" typically
        time_units = zg_temp["time"].attrs["units"]
        if occursin("hours since", time_units)
            epoch_match = match(r"hours since (\d{4}-\d{2}-\d{2})", time_units)
            if !isnothing(epoch_match)
                epoch = DateTime(epoch_match.captures[1])
                time_values = [epoch + Hour(Int(h)) for h in time_values_raw]
            else
                error("Could not parse time units: $(time_units)")
            end
        else
            error("Unexpected time units: $(time_units)")
        end

        # Find nearest lat/lon indices
        println("  Finding nearest grid point...")
        lat_idx = find_nearest_index(lat_values, lat)
        lon_idx = find_nearest_index(lon_values, lon)

        selected_lat = lat_values[lat_idx]
        selected_lon = lon_values[lon_idx]

        println("    Requested: $(lat)°N, $(lon)°E")
        println("    Nearest: $(selected_lat)°N, $(selected_lon)°E")

        # Find time range indices
        time_start_idx = findfirst(t -> t >= start_time, time_values)
        time_end_idx = findlast(t -> t <= end_time, time_values)

        if isnothing(time_start_idx) || isnothing(time_end_idx)
            error("Requested time range not available in dataset")
        end

        time_slice = time_start_idx:time_end_idx
        selected_times = time_values[time_slice]

        println("  Extracting $(length(selected_times)) time steps in parallel...")

        # Extract data for selected point and time range in parallel
        # ERA5-Land ARCO Zarr actual storage order: (longitude, latitude, time)
        # Even though _ARRAY_DIMENSIONS says ["time", "latitude", "longitude"]
        data = @sync begin
            t1 = @spawn Float64.(t2m_zarr[lon_idx, lat_idx, time_slice])
            t2 = @spawn Float64.(d2m_zarr[lon_idx, lat_idx, time_slice])
            t3 = @spawn Float64.(sp_zarr[lon_idx, lat_idx, time_slice])
            t4 = @spawn Float64.(tp_zarr[lon_idx, lat_idx, time_slice])
            t5 = @spawn Float64.(u10_zarr[lon_idx, lat_idx, time_slice])
            t6 = @spawn Float64.(v10_zarr[lon_idx, lat_idx, time_slice])
            t7 = @spawn Float64.(ssrd_zarr[lon_idx, lat_idx, time_slice])
            t8 = @spawn Float64.(strd_zarr[lon_idx, lat_idx, time_slice])
            (fetch(t1), fetch(t2), fetch(t3), fetch(t4),
             fetch(t5), fetch(t6), fetch(t7), fetch(t8))
        end
        temperature_air, temperature_dewpoint, pressure_air, precipitation_raw,
            u10, v10, shortwave_raw, longwave_raw = data

        println("  Converting units...")

        # Convert dewpoint to vapor pressure
        vapor_pressure = dewpoint_to_vapor_pressure(temperature_dewpoint)

        # Convert precipitation: m -> kg/m²
        precipitation = precipitation_raw .* 1000.0
        precipitation[precipitation .< 0] .= 0.0  # Remove numerical noise

        # Wind speed from components
        wind_speed = hypot.(u10, v10)

        # Convert radiation: J/m² (hourly accumulation) -> W/m²
        shortwave_downward = shortwave_raw ./ 3600.0
        shortwave_downward[shortwave_downward .< 0] .= 0.0  # Remove numerical noise
        longwave_downward = longwave_raw ./ 3600.0

        println("  Data loaded successfully!")
        println("  Number of time steps: $(length(selected_times))")

        # Create DimArrays with time dimension
        time_dim = Ti(selected_times)

        # Build DimStack with all climate variables
        stack = DimStack((
            temperature_air = DimArray(temperature_air, (time_dim,);
                                      metadata=Dict("units" => "K", "long_name" => "2m air temperature")),
            pressure_air = DimArray(pressure_air, (time_dim,);
                                   metadata=Dict("units" => "Pa", "long_name" => "surface pressure")),
            vapor_pressure = DimArray(vapor_pressure, (time_dim,);
                                     metadata=Dict("units" => "Pa", "long_name" => "vapor pressure")),
            wind_speed = DimArray(wind_speed, (time_dim,);
                                 metadata=Dict("units" => "m/s", "long_name" => "10m wind speed")),
            precipitation = DimArray(precipitation, (time_dim,);
                                    metadata=Dict("units" => "kg/m²/hr", "long_name" => "total precipitation")),
            shortwave_downward = DimArray(shortwave_downward, (time_dim,);
                                         metadata=Dict("units" => "W/m²", "long_name" => "surface solar radiation downward")),
            longwave_downward = DimArray(longwave_downward, (time_dim,);
                                        metadata=Dict("units" => "W/m²", "long_name" => "surface thermal radiation downward"))
        ); metadata = Dict(
            "latitude" => selected_lat,
            "longitude" => selected_lon,
            "dataset" => "ERA5-Land",
            "chunk_strategy" => string(chunk_strategy),
            "temperature_air_mean" => Statistics.mean(temperature_air),
            "wind_speed_mean" => Statistics.mean(wind_speed),
            "precipitation_mean" => Statistics.mean(precipitation) * 8760.0,  # hourly to annual
            "temperature_observation_height" => 2.0,  # ERA5-Land 2m temperature
            "wind_observation_height" => 10.0  # ERA5-Land 10m wind
        ))

        # Validate units before returning
        println("  Validating units...")
        validate_climate_forcing_units(stack)
        println("  ✓ All units validated")

        println("  DimStack created successfully!")
        return stack

    catch e
        if isa(e, HTTP.Exceptions.StatusError) || (isa(e, ErrorException) && occursin("HTTP", e.msg))
            error("HTTP error accessing ERA5-Land. Check your token and network connection.\nError: $e")
        else
            rethrow()
        end
    end
end
