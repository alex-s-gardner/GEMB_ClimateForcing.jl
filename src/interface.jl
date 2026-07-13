"""
    climate_forcing(dataset::Symbol, lat::Real, lon::Real; kwargs...) -> ClimateForcing

Load climate forcing data from specified dataset and return GEMB-compatible ClimateForcing struct.

# Arguments
- `dataset::Symbol`: Dataset identifier
  - `:era5land` - ERA5-Land reanalysis (0.1° resolution, 1950-present)
  - `:era5` - ERA5 reanalysis (0.25° resolution) [future]
  - `:merra2` - MERRA-2 reanalysis [future]
- `lat::Real`: Latitude in degrees [-90, 90]
- `lon::Real`: Longitude in degrees [-180, 180] or [0, 360]

# Keyword Arguments
- `time_range::Tuple{DateTime,DateTime}`: Time range to extract (required)
- `token::Union{String,Nothing}=nothing`: API token/key for authentication (required for ERA5-Land)
- `chunk_strategy::Symbol=:geo`: Chunking strategy
  - `:geo` - Geo-chunked (optimized for time-series at a point)
  - `:time` - Time-chunked (optimized for spatial maps)
- `kwargs...`: Dataset-specific keyword arguments

# Returns
- `ClimateForcing`: Struct ready for `gemb(profile, cf, mp)`

# Examples
```julia
# Load ERA5-Land for Summit, Greenland
cf = climate_forcing(
    :era5land, 72.58, -38.46;
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"]
)

# Use with GEMB
using GEMB
mp = ModelParameters(output_frequency=:daily)
profile = initialize_profile(mp, cf)
output = gemb(profile, cf, mp)
```

# Notes
- For ERA5-Land, obtain a free API key from https://cds.climate.copernicus.eu/
- Store token in environment variable: `export CDS_API_KEY="your-token-here"`
- Geo-chunked strategy is recommended for single-point time-series extraction
"""
function climate_forcing(
    dataset::Symbol,
    lat::Real,
    lon::Real;
    time_range::Union{Tuple{DateTime,DateTime},Nothing}=nothing,
    token::Union{String,Nothing}=nothing,
    chunk_strategy::Symbol=:time,
    kwargs...
)
    # Validate required arguments
    if isnothing(time_range)
        throw(ArgumentError("time_range keyword argument is required"))
    end

    # Validate time_range
    start_time, end_time = time_range
    if start_time >= end_time
        throw(ArgumentError("time_range start must be before end"))
    end

    # Validate lat/lon ranges
    if lat < -90 || lat > 90
        throw(ArgumentError("lat must be in range [-90, 90]"))
    end
    if lon < -180 || lon > 360
        throw(ArgumentError("lon must be in range [-180, 180] or [0, 360]"))
    end

    # Validate chunk_strategy
    if !(chunk_strategy in (:geo, :time))
        throw(ArgumentError("chunk_strategy must be :geo or :time"))
    end

    # Dispatch to dataset-specific loader
    if dataset == :era5land
        return load_era5_land(lat, lon; time_range=time_range, token=token, chunk_strategy=chunk_strategy, kwargs...)
    elseif dataset == :era5
        throw(ArgumentError("ERA5 dataset not yet implemented"))
    elseif dataset == :merra2
        throw(ArgumentError("MERRA-2 dataset not yet implemented"))
    else
        throw(ArgumentError("Unsupported dataset: $dataset. Supported: :era5land"))
    end
end
