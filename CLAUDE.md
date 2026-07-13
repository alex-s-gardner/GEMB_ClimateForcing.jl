# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB_ClimateForcing.jl loads climate forcing data from reanalysis datasets and returns [GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl)-compatible `ClimateForcing` structs. The package uses **pure Julia** (no Python dependencies) to access cloud-optimized ARCO Zarr stores via authenticated HTTPS.

## Development Commands

### Setup
```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Set up CDS API key for ERA5-Land access
export CDS_API_KEY="your-token-here"
```

### Testing
```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run tests with CDS API key (includes integration tests)
export CDS_API_KEY="your-token-here"
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Running Examples
```bash
# Complete ERA5-Land + GEMB simulation example
export CDS_API_KEY="your-token-here"
julia --project=. examples/era5_land_example.jl

# Test authentication only
julia --project=. examples/test_authentication.jl
```

## Architecture

### Core Components

1. **`src/interface.jl`** - Unified entry point
   - `climate_forcing(dataset::Symbol, lat, lon; kwargs...)` - Main API
   - Input validation (lat/lon ranges, time_range, chunk_strategy)
   - Dispatches to dataset-specific loaders

2. **`src/authenticated_http_store.jl`** - Custom Zarr store
   - `AuthenticatedHTTPStore <: Zarr.AbstractStore`
   - Follows Zarr.jl's `GCStore` pattern for Bearer token authentication
   - Passes `Authorization: Bearer <token>` header in HTTP requests
   - Read-only store for HTTPS access to cloud-optimized Zarr

3. **`src/datasets/era5_land.jl`** - ERA5-Land implementation
   - `load_era5_land()` - Loads from ECMWF ARCO Zarr stores
   - Accesses 4 variable groups (temperature, pressure/precip, wind, radiation)
   - Nearest-neighbor lat/lon selection
   - Unit conversions (J/m² → W/m², m → kg/m², dewpoint → vapor pressure)
   - Returns `GEMB.ClimateForcing` struct with all required variables

### Data Flow

```
climate_forcing(:era5land, lat, lon; time_range, token)
    ↓
Create AuthenticatedHTTPStore for each ERA5-Land variable group
    ↓
Zarr.zopen(store, consolidated=true) - load Zarr metadata
    ↓
Find nearest lat/lon indices in coordinate arrays
    ↓
Slice time range and extract data at point
    ↓
Convert units and compute derived variables (e.g., wind_speed from u10/v10)
    ↓
GEMB.initialize_forcing() - create ClimateForcing struct
    ↓
Return ClimateForcing with DimArray variables
```

### Pure Julia Design

This package deliberately avoids PythonCall/xarray in favor of pure Julia:
- **Zarr.jl** (v0.10+) for reading cloud Zarr stores
- **AuthenticatedHTTPStore** implements Bearer token auth (pattern from GCStore)
- **HTTP.jl + OpenSSL.jl** for HTTPS requests with custom headers
- **DimensionalData.jl** for dimension-aware indexing (already used by GEMB.jl)

### ERA5-Land Specifics

**Variable Groups and ARCO Store IDs:**
- `sfc-2m-temperature` (store 007): `t2m`, `d2m`
- `sfc-pressure-precipitation` (store 009): `sp`, `tp`
- `sfc-wind` (store 008): `u10`, `v10`
- `sfc-radiation-heat` (store 010): `ssrd`, `strd`

**Chunk Strategies:**
- `:geo` (default) - Geo-chunked stores, optimized for time-series at a point
- `:time` - Time-chunked stores, optimized for spatial maps

**Authentication:**
Requires free CDS API key from https://cds.climate.copernicus.eu/

## Adding New Datasets

To add support for a new reanalysis dataset (e.g., ERA5, MERRA-2):

1. Create `src/datasets/your_dataset.jl`
2. Implement `load_your_dataset(lat, lon; time_range, token, kwargs...)`
3. Return a `GEMB.ClimateForcing` struct with required variables:
   - `temperature_air`, `pressure_air`, `vapor_pressure`
   - `wind_speed`, `precipitation`
   - `shortwave_downward`, `longwave_downward`
4. Add dispatch case in `src/interface.jl`:
   ```julia
   elseif dataset == :yourdataset
       return load_your_dataset(lat, lon; time_range=time_range, token=token, kwargs...)
   ```
5. Update README and tests

Use `src/datasets/era5_land.jl` as a template.

## Important Implementation Notes

### Coordinate System Handling
- ERA5-Land uses lat/lon coordinates (not projected)
- Latitude: -90 to 90 (south to north)
- Longitude: -180 to 180 or 0 to 360 (both accepted, normalized internally)
- Nearest-neighbor selection used for point extraction

### Unit Conversions
ERA5-Land variables require conversion for GEMB compatibility:
- `tp` (m) → `precipitation` (kg/m²): multiply by 1000
- `ssrd`, `strd` (J/m²) → `shortwave/longwave_downward` (W/m²): divide by 3600
- `d2m` (K) → `vapor_pressure` (Pa): via dewpoint formula
- `u10`, `v10` (m/s) → `wind_speed` (m/s): magnitude √(u² + v²)

### Performance Characteristics
- First load: ~30-60 seconds (network dependent, loads metadata from 4 stores)
- Subsequent loads: faster due to HTTP caching
- Memory: only requested time/location downloaded (lazy loading)
- Geo-chunked is 2-5x faster than time-chunked for point time-series

### Testing Strategy
Tests use conditional integration testing:
- Basic tests (input validation, error handling) run without credentials
- Integration tests (actual data loading) only run if `CDS_API_KEY` is set
- Use 1-day time ranges for fast integration tests

## JuliaGeo Integration

While this package doesn't currently use JuliaGeo directly, it follows JuliaGeo patterns:
- Uses DimensionalData.jl for dimension-aware arrays
- Ready for GeoInterface.jl integration if spatial domains are added
- Coordinate reference system: WGS84 (EPSG:4326) for ERA5-Land
