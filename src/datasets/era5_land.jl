"""
ERA5-Land ARCO (Analysis-Ready, Cloud-Optimized) Zarr data loader.

Pure Julia implementation using Zarr.jl with authenticated HTTP access.
Loads data from ECMWF's cloud-optimized Zarr stores via HTTPS with Bearer token.
"""

using Zarr
using DimensionalData
using Dates
using Statistics

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
    validate_climate_forcing_units(;
        temperature_air, pressure_air, precipitation, wind_speed,
        shortwave_downward, longwave_downward, vapor_pressure
    )

Validate that climate forcing variables have physically reasonable values.

# Expected Units
- `temperature_air`: Kelvin (K), range [180, 330]
- `pressure_air`: Pascal (Pa), range [30000, 110000]
- `precipitation`: kg/m², range [0, 100] per hour
- `wind_speed`: m/s, range [0, 100]
- `shortwave_downward`: W/m², range [0, 1500]
- `longwave_downward`: W/m², range [50, 500]
- `vapor_pressure`: Pascal (Pa), range [0, 10000]

# Throws
- `ArgumentError` if any variable has values outside expected physical ranges
"""
function validate_climate_forcing_units(;
    temperature_air::AbstractVector,
    pressure_air::AbstractVector,
    precipitation::AbstractVector,
    wind_speed::AbstractVector,
    shortwave_downward::AbstractVector,
    longwave_downward::AbstractVector,
    vapor_pressure::AbstractVector
)
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
    load_era5_land(lat::Real, lon::Real; kwargs...) -> ClimateForcing

Load ERA5-Land reanalysis data using pure Julia (Zarr.jl) and return GEMB-compatible ClimateForcing struct.

# Keyword Arguments
- `time_range::Tuple{DateTime,DateTime}`: Time range to extract (required)
- `token::Union{String,Nothing}`: CDS API key (required)
- `chunk_strategy::Symbol=:geo`: :geo (time-series) or :time (spatial)

# Returns
- `ClimateForcing`: Struct with all required forcing variables

# Notes
- ERA5-Land is on 0.1° grid (approximately 9 km)
- Uses nearest neighbor selection for lat/lon
- Hourly temporal resolution
- Requires free CDS API key from https://cds.climate.copernicus.eu/
- Pure Julia implementation using Zarr.jl + DimensionalData.jl
"""
function load_era5_land(
    lat::Real,
    lon::Real;
    time_range::Tuple{DateTime,DateTime},
    token::Union{String,Nothing}=nothing,
    chunk_strategy::Symbol=:geo,
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

        store_temp = AuthenticatedHTTPStore(url_temp; token=token)
        store_precip = AuthenticatedHTTPStore(url_precip; token=token)
        store_wind = AuthenticatedHTTPStore(url_wind; token=token)
        store_rad = AuthenticatedHTTPStore(url_rad; token=token)

        # Open Zarr groups
        println("  Opening Zarr groups...")
        zg_temp = Zarr.zopen(store_temp, consolidated=true, fill_as_missing=false)
        zg_precip = Zarr.zopen(store_precip, consolidated=true, fill_as_missing=false)
        zg_wind = Zarr.zopen(store_wind, consolidated=true, fill_as_missing=false)
        zg_rad = Zarr.zopen(store_rad, consolidated=true, fill_as_missing=false)

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

        println("    Requested: $(lat)°N, $(lon)°E")
        println("    Nearest: $(lat_values[lat_idx])°N, $(lon_values[lon_idx])°E")

        # Find time range indices
        time_start_idx = findfirst(t -> t >= start_time, time_values)
        time_end_idx = findlast(t -> t <= end_time, time_values)

        if isnothing(time_start_idx) || isnothing(time_end_idx)
            error("Requested time range not available in dataset")
        end

        time_slice = time_start_idx:time_end_idx
        selected_times = time_values[time_slice]

        println("  Extracting $(length(selected_times)) time steps...")

        # Extract data for selected point and time range
        # ERA5-Land ARCO Zarr actual storage order: (longitude, latitude, time)
        # Even though _ARRAY_DIMENSIONS says ["time", "latitude", "longitude"]
        temperature_air = Float64.(t2m_zarr[lon_idx, lat_idx, time_slice])
        temperature_dewpoint = Float64.(d2m_zarr[lon_idx, lat_idx, time_slice])
        pressure_air = Float64.(sp_zarr[lon_idx, lat_idx, time_slice])
        precipitation_raw = Float64.(tp_zarr[lon_idx, lat_idx, time_slice])
        u10 = Float64.(u10_zarr[lon_idx, lat_idx, time_slice])
        v10 = Float64.(v10_zarr[lon_idx, lat_idx, time_slice])
        shortwave_raw = Float64.(ssrd_zarr[lon_idx, lat_idx, time_slice])
        longwave_raw = Float64.(strd_zarr[lon_idx, lat_idx, time_slice])

        println("  Converting units...")

        # Convert dewpoint to vapor pressure
        vapor_pressure = GEMB.dewpoint_to_vapor_pressure(temperature_dewpoint)

        # Convert precipitation: m -> kg/m²
        precipitation = precipitation_raw .* 1000.0
        precipitation[precipitation .< 0] .= 0.0  # Remove numerical noise

        # Wind speed from components
        wind_speed = hypot.(u10, v10)

        # Convert radiation: J/m² (hourly accumulation) -> W/m²
        shortwave_downward = shortwave_raw ./ 3600.0
        longwave_downward = longwave_raw ./ 3600.0

        println("  Data loaded successfully!")
        println("  Number of time steps: $(length(selected_times))")

        # Validate units before creating ClimateForcing
        println("  Validating units...")
        validate_climate_forcing_units(
            temperature_air=temperature_air,
            pressure_air=pressure_air,
            precipitation=precipitation,
            wind_speed=wind_speed,
            shortwave_downward=shortwave_downward,
            longwave_downward=longwave_downward,
            vapor_pressure=vapor_pressure
        )
        println("  ✓ All units validated")

        # Create ClimateForcing struct using GEMB's initialize_forcing
        cf = GEMB.initialize_forcing(
            selected_times,
            temperature_air,
            pressure_air,
            precipitation,
            wind_speed,
            shortwave_downward,
            longwave_downward,
            vapor_pressure;
            temperature_air_mean=Statistics.mean(temperature_air),
            wind_speed_mean=Statistics.mean(wind_speed),
            precipitation_mean=Statistics.mean(precipitation) * 8760.0,  # hourly to annual
            temperature_observation_height=2.0,  # ERA5-Land 2m temperature
            wind_observation_height=10.0  # ERA5-Land 10m wind
        )

        println("  ClimateForcing struct created successfully!")
        return cf

    catch e
        if isa(e, HTTP.Exceptions.StatusError) || (isa(e, ErrorException) && occursin("HTTP", e.msg))
            error("HTTP error accessing ERA5-Land. Check your token and network connection.\nError: $e")
        else
            rethrow()
        end
    end
end
