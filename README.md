# GEMB_ClimateForcing.jl

Load climate forcing data from various reanalysis datasets and return a `DimStack` with climate variables. Seamlessly converts to [GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl) `ClimateForcing` via package extension.

## Features

- **Pure Julia**: No Python dependencies - native Zarr.jl with authenticated HTTP access
- **Unified Interface**: Single `climate_forcing()` function for all datasets
- **Cloud-Optimized**: Direct access to ARCO (Analysis-Ready, Cloud-Optimized) Zarr stores via HTTPS
- **Authenticated Access**: Bearer token authentication following Zarr.jl's GCStore pattern
- **No Downloads**: Lazy loading of only requested time ranges and locations
- **DimensionalData Integration**: Zarr arrays automatically work with DimStack indexing
- ⚡ **Parallel Loading**: Concurrent access to multiple variable groups for optimal performance
- 🚀 **Chunk Caching**: Optional persistent disk cache via Zarr.CachingStore
- **Extensible**: Easy to add new datasets (ERA5, MERRA-2, JRA-55, etc.)

## Supported Datasets

| Dataset | Symbol | Resolution | Coverage | Status |
|---------|--------|------------|----------|--------|
| ERA5-Land | `:era5land` | 0.1° (~9 km) | 1950-present | ✅ Implemented |
| ERA5 | `:era5` | 0.25° (~25 km) | 1940-present | 🔜 Planned |
| MERRA-2 | `:merra2` | 0.5° × 0.625° | 1980-present | 🔜 Planned |
| JRA-55 | `:jra55` | 1.25° | 1958-present | 🔜 Planned |

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/GEMB_ClimateForcing.jl")
```

## Quick Start

### 1. Get a CDS API Key

ERA5-Land requires a free API key from the Copernicus Climate Data Store:

1. Register at [https://cds.climate.copernicus.eu/](https://cds.climate.copernicus.eu/)
2. Go to your profile page
3. Copy your API key
4. Set environment variable:
   ```bash
   export CDS_API_KEY="your-token-here"
   ```

### 2. Load Climate Forcing

```julia
using GEMB_ClimateForcing
using GEMB  # Extension automatically provides DimStack → ClimateForcing conversion
using Dates

# Load ERA5-Land for Summit, Greenland (returns DimStack)
forcing_data = climate_forcing(
    :era5land,           # Dataset
    72.58, -38.46,       # Latitude, Longitude
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"]
)

# Convert to GEMB.ClimateForcing (extension method)
cf = GEMB.ClimateForcing(forcing_data)
```

### 3. Run GEMB

```julia
# Initialize model
mp = GEMB.ModelParameters(output_frequency=:daily)
profile = GEMB.initialize_profile(mp, cf)

# Run simulation
output = GEMB.gemb(profile, cf, mp)
```

## API Reference

### `climate_forcing(dataset, lat, lon; kwargs...)`

Load climate forcing data from specified dataset.

**Arguments:**
- `dataset::Symbol` - Dataset identifier (`:era5land`, `:era5`, etc.)
- `lat::Real` - Latitude in degrees [-90, 90]
- `lon::Real` - Longitude in degrees [-180, 180] or [0, 360]

**Keyword Arguments:**
- `time_range::Tuple{DateTime,DateTime}` - Time range to extract (required)
- `token::Union{String,Nothing}` - API token for authentication (required for ERA5-Land)
- `chunk_strategy::Symbol=:geo` - Chunking strategy (`:geo` or `:time`)
  - `:geo` - Optimized for time-series extraction at a point (recommended)
  - `:time` - Optimized for spatial maps
- `cache_path::Union{String,Nothing}=nothing` - Path for persistent disk cache (Zarr.CachingStore)

**Returns:**
- `DimStack` - Stack with climate forcing variables as DimArrays:
  - `temperature_air`, `pressure_air`, `vapor_pressure`, `wind_speed`,
    `precipitation`, `shortwave_downward`, `longwave_downward`
  - Metadata includes location info and observation heights

**Example:**
```julia
using GEMB_ClimateForcing
using GEMB

# Load data (returns DimStack)
forcing_data = climate_forcing(
    :era5land, 72.58, -38.46;
    time_range=(DateTime(2020,1,1), DateTime(2021,1,1)),
    token=ENV["CDS_API_KEY"],
    chunk_strategy=:geo
)

# Convert to GEMB.ClimateForcing (requires GEMB.jl loaded)
cf = GEMB.ClimateForcing(forcing_data)
```

## ERA5-Land Details

### Variables

GEMB_ClimateForcing automatically loads and converts the following ERA5-Land variables:

| ERA5-Land Variable | GEMB Variable | Units | Conversion |
|-------------------|---------------|-------|------------|
| `t2m` | `temperature_air` | K | Direct |
| `d2m` | `vapor_pressure` | Pa | Via dewpoint formula |
| `sp` | `pressure_air` | Pa | Direct |
| `tp` | `precipitation` | kg/m² | Multiply by 1000 |
| `u10`, `v10` | `wind_speed` | m/s | Magnitude |
| `ssrd` | `shortwave_downward` | W/m² | Divide by 3600 |
| `strd` | `longwave_downward` | W/m² | Divide by 3600 |

### Data Source

ERA5-Land data is accessed from ECMWF's cloud-optimized ARCO Zarr stores:
```
https://arco.datastores.ecmwf.int/cadl-arco-{geo|time}-{store-id}/arco/reanalysis_era5_land/...
```

### Pure Julia Implementation

GEMB_ClimateForcing uses pure Julia components:
- **Zarr.jl** (v0.10+) for reading cloud-optimized Zarr arrays
- **AuthenticatedHTTPStore** - Custom store following GCStore pattern for Bearer token auth
- **HTTP.jl** + **OpenSSL.jl** for HTTPS requests with custom headers
- **DimensionalData.jl** for dimension-aware array indexing (already used by GEMB.jl)

No PythonCall or xarray required!

### Performance

**Characteristics:**
- First data load: ~10-25 seconds (1 year of hourly data)
- Parallel loading of 4 variable groups provides 1.5-2x speedup
- **Disk caching**: With `cache_path`, subsequent loads skip network entirely
- Memory: only requested time/location downloaded (lazy loading)

**Tips:**
- **Use geo-chunked strategy** (default) for extracting time-series at single points
- **Use time-chunked strategy** when extracting spatial maps
- **Enable disk caching** with `cache_path` for repeated access to the same data across sessions
- Cache persists on disk - subsequent Julia sessions benefit without re-downloading

## Examples

See `examples/era5_land_example.jl` for a complete working example including:
- Loading ERA5-Land forcing
- Running GEMB with spinup
- Analyzing results

Run with:
```bash
export CDS_API_KEY="your-token-here"
julia --project=. examples/era5_land_example.jl
```

## Adding New Datasets

To add support for a new dataset:

1. Create `src/datasets/your_dataset.jl`
2. Implement `load_your_dataset(lat, lon; kwargs...)` function
3. Add to dispatch in `src/interface.jl`
4. Update README with dataset details

See `src/datasets/era5_land.jl` for a template implementation.

## Citation

If you use ERA5-Land data, please cite:

> Muñoz Sabater, J., (2019): ERA5-Land hourly data from 1950 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: [10.24381/cds.e2161bac](https://doi.org/10.24381/cds.e2161bac)

## License

MIT License - see LICENSE file for details

## Related Projects

- [GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl) - Glacier Energy and Mass Balance model
- [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land) - Dataset information
