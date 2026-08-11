# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB_ClimateForcing.jl provides climate forcing for the GEMB surface mass/energy balance model as a `DimStack` of climate variables. It covers four largely independent workflows, all pure Julia (no Python):

1. **Reanalysis loading** — point extraction from cloud-optimized ARCO Zarr stores (ERA5-Land) via authenticated HTTPS (`climate_forcing`).
2. **Elevation downscaling** — physically-based per-variable correction of a forcing stack to a target elevation (`climate_adjust_for_elevation`).
3. **Static/invariant fields** — lazy `Raster`s for land-sea mask, geopotential, vegetation, and a 30 m global DEM (`climate_model_invariant`).
4. **Synthetic forcing + fitting** — generate stochastic forcing from parameter sets, and fit those parameters from observed data (`simulate_climate_forcing`, `fit_*`). These are Julia translations of GEMB's MATLAB `simulate_*` / `fit_*` functions and validate against them.

**No GEMB dependency**: GEMB.jl provides a package extension that converts the `DimStack` to `GEMB.ClimateForcing`. The conversion code lives in GEMB.jl, *not* here — there is no `ext/` directory in this repo. This lets the forcing be used by other models without requiring GEMB.

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

# Run invariant-parameter download tests (no CDS token needed, fetches ~50 MB/file)
GEMB_TEST_INVARIANT=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

**Test gating.** `test/runtests.jl` always runs offline tests (input validation, unit
conversions, elevation physics, synthetic forcing, chunk-map formulas, DEM tile
geometry). Network-dependent tests are opt-in via environment variables so CI stays
offline by default:
- `CDS_API_KEY` — enables ERA5-Land Zarr integration tests
- `GEMB_TEST_INVARIANT=1` — enables ERA5-Land invariant NetCDF downloads
- `GEMB_TEST_CHUNK_MAP=1` — enables `climate_chunk_map` integration
- `GEMB_TEST_COPERNICUS_DEM=1` — enables Copernicus DEM `/vsicurl/` reads

There is no per-file test runner; run a single suite by editing `include`s in
`test/runtests.jl` or invoking the file directly with the package loaded.

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
   - `climate_chunk_map(dataset; chunk_strategy, ...)` - returns a global `Raster{Int64}` of
     Zarr chunk IDs (downloads only coordinate/chunk metadata) for visualizing download
     locality before batch queries

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

4. **`src/utils.jl`** - Vapor-pressure conversions shared across all workflows
   - `dewpoint_to_vapor_pressure()`, `vapor_pressure_to_relative_humidity()`,
     `relative_humidity_to_vapor_pressure()` — Magnus/Buck formulas

5. **`src/elevation_adjustment.jl`** - Elevation downscaling (workflow 2)
   - `climate_adjust_for_elevation(stack, Δz; lapse_rate, precip_scaling_method)` — per-variable
     correction of a forcing `DimStack` to a target elevation (Glover 1999 scheme). `Δz = 0`
     reproduces the input; physical-range validation re-runs on the result.
   - `empirical_lapse_rate()` fits a local rate from neighbouring grid cells (RACMO-style).
   - Exported climatological monthly tables: `GREENLAND_LAPSE_RATE`, `ARCTIC_LAPSE_RATE`,
     `ANTARCTICA_LAPSE_RATE` (K/km, Jan→Dec). See README for the full per-variable table and references.

6. **`src/invariant.jl`** - Time-invariant fields (workflow 3)
   - `climate_model_invariant(; model, parameter, ...)` — returns lazy `Raster`/`RasterStack`.
   - Two dispatch paths: **file-registry models** (`:era5_land`) download single-field NetCDFs
     keyed by GRIB shortName in `ERA5_LAND_INVARIANT_PARAMETERS`; **extent-based models**
     (`:copernicus_dem_30m`, in `_INVARIANT_EXTENT_MODELS`) delegate to a tiled loader.
   - Invariant grid is 0–359.9°E, latitude **descending** (90→−90°N). Geopotential `z` is
     m²/s²; orography = `z / 9.80665`.

7. **`src/datasets/copernicus_dem.jl`** - Copernicus GLO-30 DEM loader
   - Reads 1°×1° Cloud-Optimized GeoTIFF tiles from AWS Open Data via GDAL `/vsicurl/` byte-range
     reads (tiles never downloaded in full); mosaics covering tiles into an on-the-fly VRT.
   - Tile geometry is analytical (derived from tile id + published `tileList.txt`), so building
     the global mosaic opens zero tiles.
   - `_configure_gdal_http()` points GDAL's curl at Julia's CA bundle (`NetworkOptions.ca_roots_path()`)
     — required or `/vsicurl/` TLS handshakes fail on macOS.

8. **`src/simulate/simulate_climate_forcing.jl`** - Synthetic forcing (workflow 4)
   - `simulate_climate_forcing(set_id, time_step_hours=0)` — generates a full stochastic forcing
     `DimStack` from a named parameter set (`simulation_parameter_sets`, e.g. `"test_1"`), seeded
     RNG (`MersenneTwister`) for MATLAB-matching reproducibility.
   - Time↔decimal-year helpers `datetime2decyear` / `decyear2datetime` live here.
   - `src/forcing_climatology.jl` — `forcing_climatology(ds[, range])` averages complete years
     (drops leap day 366 and partial years) into a one-year climatological forcing.

9. **`src/fit_climate/`** - Parameter fitting (workflow 4, inverse of simulate)
   - `fit_air_temperature`, `fit_precipitation`, `fit_longwave_irradiance_delta`,
     `fit_seasonal_daily_noise` — estimate `simulate_*` coefficients from observed series.
   - `simulate_coeffs_disp` prints a fitted coefficient NamedTuple as copy-pasteable Julia.
   - `varname2longname` maps variable symbols to human-readable names.

### MATLAB parity

Workflow-4 functions are direct ports of GEMB's MATLAB `simulate_*` / `fit_*` code and are
expected to reproduce its numerics. When editing them, preserve the RNG seeding order and the
exact coefficient math — divergence from MATLAB is a regression, not a cleanup opportunity.

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

- **DimensionalData.jl** — all forcing output is a `DimStack`; the `Ti` (time) dimension is
  required by `climate_adjust_for_elevation` and `forcing_climatology`.
- **Rasters.jl** — invariant fields and the DEM are returned as lazy `Raster`/`RasterStack`
  (backends: `RastersNCDatasetsExt` via NCDatasets, `RastersArchGDALExt` via ArchGDAL).
- **GDAL** (via ArchGDAL) — `/vsicurl/` remote COG reads for the Copernicus DEM.
- Coordinate reference system: WGS84 (EPSG:4326) throughout.

### Grid-orientation gotcha
Both the ERA5-Land ARCO grid and the invariant NetCDFs use **0–359.9°E longitude** and
**descending latitude (90→−90°N)**. When cropping a `Raster`, `X = a .. b` / `Y = a .. b`
takes `min .. max` regardless of storage order — but verify axis order after any reprojection
and never assume ascending latitude.
