# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB_ClimateForcing.jl loads climate forcing data from reanalysis datasets and returns a `DimStack` with climate variables. The package uses **pure Julia** (no Python dependencies) to access cloud-optimized ARCO Zarr stores via authenticated HTTPS.

**Architecture**: GEMB_ClimateForcing has **no dependency on GEMB**. Instead, GEMB.jl provides a package extension (`ext/GEMBClimateForcing.jl`) that converts the DimStack to `GEMB.ClimateForcing`. This allows the forcing data to be used by other models or tools without requiring GEMB.

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
   - Returns `DimStack` with all climate variables as DimArrays
   - Physical validation via `validate_climate_forcing_units()`

4. **`src/utils.jl`** - Utility functions
   - `dewpoint_to_vapor_pressure()` - Magnus formula for vapor pressure calculation

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
Create DimStack with all variables + metadata
    ↓
Validate physical ranges (validate_climate_forcing_units)
    ↓
Return DimStack
    
[Optional: if GEMB.jl is loaded]
    ↓
GEMB.ClimateForcing(dimstack) - extension method
    ↓
GEMB-specific validation and conversion
    ↓
Return GEMB.ClimateForcing struct
```

### Pure Julia Design

This package deliberately avoids PythonCall/xarray in favor of pure Julia:
- **Zarr.jl** (v0.10+) for reading cloud Zarr stores
- **AuthenticatedHTTPStore** implements Bearer token auth (pattern from GCStore)
- **HTTP.jl + OpenSSL.jl** for HTTPS requests with custom headers
- **DimensionalData.jl** for dimension-aware indexing and DimStack output
- **No GEMB dependency** - GEMB.jl provides conversion via package extension

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
3. Return a `DimStack` with required variables:
   - `temperature_air`, `pressure_air`, `vapor_pressure`
   - `wind_speed`, `precipitation`
   - `shortwave_downward`, `longwave_downward`
   - Metadata with: `latitude`, `longitude`, `temperature_air_mean`, 
     `wind_speed_mean`, `precipitation_mean`, 
     `temperature_observation_height`, `wind_observation_height`
4. Call `validate_climate_forcing_units(stack)` before returning
5. Add dispatch case in `src/interface.jl`:
   ```julia
   elseif dataset == :yourdataset
       return load_your_dataset(lat, lon; time_range=time_range, token=token, kwargs...)
   ```
6. Update README and tests

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
- First load: ~10-25 seconds (network dependent, uses parallel loading)
- Subsequent loads: faster due to HTTP caching
- Memory: only requested time/location downloaded (lazy loading)
- Geo-chunked is 2-5x faster than time-chunked for point time-series
- **Parallel loading:** 4 variable groups loaded concurrently for 1.5-2x speedup

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
