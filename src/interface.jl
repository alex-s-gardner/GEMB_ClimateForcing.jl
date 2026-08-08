"""
    climate_forcing(dataset::Symbol, lat::Real, lon::Real; kwargs...) -> DimStack

Load climate forcing data from specified dataset and return a DimStack with climate variables.

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
- `cache_path::Union{String,Nothing}=nothing`: Path for persistent disk cache (Zarr.CachingStore)
- `kwargs...`: Dataset-specific keyword arguments

# Returns
- `DimStack`: Stack with climate forcing variables as DimArrays
  - `temperature_air`, `pressure_air`, `vapor_pressure`, `wind_speed`,
    `precipitation`, `shortwave_downward`, `longwave_downward`
  - Metadata includes location, dataset info, and observation heights

# Examples
```julia
using GEMB_ClimateForcing
using GEMB  # Automatically provides DimStack → ClimateForcing conversion

# Load ERA5-Land for Summit, Greenland
forcing_data = climate_forcing(
    :era5land, 72.58, -38.46;
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"]
)

# Convert to GEMB.ClimateForcing (requires GEMB.jl)
cf = GEMB.ClimateForcing(forcing_data)

# Use with GEMB
mp = GEMB.ModelParameters(output_frequency=:daily)
profile = GEMB.initialize_profile(mp, cf)
output = GEMB.gemb(profile, cf, mp)
```

# Notes
- For ERA5-Land, obtain a free API key from https://cds.climate.copernicus.eu/
- Store token in environment variable: `export CDS_API_KEY="your-token-here"`
- Geo-chunked strategy is recommended for single-point time-series extraction
- DimStack can be used directly or converted to GEMB.ClimateForcing via extension
"""
function climate_forcing(
    dataset::Symbol,
    lat::Real,
    lon::Real;
    time_range::Union{Tuple{DateTime,DateTime},Nothing}=nothing,
    token::Union{String,Nothing}=nothing,
    chunk_strategy::Symbol=:geo,  # :geo is optimized for point time-series (counter-intuitively!)
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

"""
    climate_chunk_map(dataset::Symbol; chunk_strategy, token, variable_group, cache_path) -> Raster{Int64}

Return a global `Raster{Int64}` where the integer value at each grid cell is the unique
spatial chunk ID for that cell in the given Zarr store layout.

Cells sharing the same chunk ID occupy the same Zarr chunk and are read together in a
single network request when that location is queried. Use this to visualize download
locality before submitting batch queries.

No climate data values are downloaded — only coordinate vectors and chunk metadata are
fetched.

# Arguments
- `dataset::Symbol`: Dataset identifier. Currently only `:era5land`.

# Keyword Arguments
- `chunk_strategy::Symbol=:geo`: `:geo` (geoChunked.zarr) or `:time` (timeChunked.zarr)
- `token::Union{String,Nothing}=nothing`: Bearer token (required for ERA5-Land)
- `variable_group::String="sfc-2m-temperature"`: ERA5-Land variable group to read chunk
  metadata from. All groups share the same spatial grid, so the default works for all cases.
- `cache_path::Union{String,Nothing}=nothing`: Zarr.CachingStore cache directory

# Returns
- `Raster{Int64}` with `X` (longitude, 0–359.9°E) and `Y` (latitude, 90→-90°N) dims.
  Chunk IDs are zero-based integers in `[0, n_lon_chunks*n_lat_chunks - 1]`. The
  `metadata` dict contains `"chunk_strategy"`, `"lon_chunk_size"`, `"lat_chunk_size"`,
  `"n_lon_chunks"`, `"n_lat_chunks"`, and `"total_spatial_chunks"`.

# Examples
```julia
using GEMB_ClimateForcing
token = ENV["CDS_API_KEY"]

geo_map  = climate_chunk_map(:era5land; chunk_strategy=:geo,  token=token)
time_map = climate_chunk_map(:era5land; chunk_strategy=:time, token=token)

# Crop to Greenland (0–360 lon convention: -38°E → 322°E)
greenland = geo_map[X=310..360, Y=59..84]
```
"""
function climate_chunk_map(
    dataset::Symbol;
    chunk_strategy::Symbol = :geo,
    token::Union{String,Nothing} = nothing,
    variable_group::String = "sfc-2m-temperature",
    cache_path::Union{String,Nothing} = nothing,
)
    if !(chunk_strategy in (:geo, :time))
        throw(ArgumentError("chunk_strategy must be :geo or :time"))
    end

    if dataset == :era5land
        return era5_land_chunk_map(;
            chunk_strategy = chunk_strategy,
            token          = token,
            variable_group = variable_group,
            cache_path     = cache_path,
        )
    elseif dataset == :era5
        throw(ArgumentError("ERA5 dataset not yet implemented"))
    elseif dataset == :merra2
        throw(ArgumentError("MERRA-2 dataset not yet implemented"))
    else
        throw(ArgumentError("Unsupported dataset: $dataset. Supported: :era5land"))
    end
end
