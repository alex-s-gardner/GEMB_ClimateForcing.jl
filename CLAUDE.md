# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB_ClimateForcing.jl provides climate forcing for the GEMB surface mass/energy balance model as a `DimStack` of climate variables. It covers five largely independent workflows, all pure Julia (no Python):

1. **Reanalysis loading** — point extraction from cloud-optimized ARCO Zarr stores (ERA5-Land) via authenticated HTTPS (`climate_forcing`).
2. **Elevation downscaling + climate perturbation** — physically-based per-variable correction of a forcing stack to a target elevation (`climate_adjust_for_elevation`), plus direct temperature/precipitation perturbations for scenario experiments (`temperature_adjust`, `precipitation_adjust`).
3. **Static/invariant fields** — lazy `Raster`s for land-sea mask, geopotential, vegetation, and a 30 m global DEM (`climate_model_invariant`).
4. **Synthetic forcing + fitting** — generate stochastic forcing from parameter sets, and fit those parameters from observed data (`simulate_climate_forcing`, `fit_*`). These are Julia translations of GEMB's MATLAB `simulate_*` / `fit_*` functions and validate against them.
5. **Satellite observations** — 10-daily C3S surface albedo, *ordered* from the CDS Retrieve API rather than read from a cloud store (`satellite_albedo`), plus its glacier bare-ice reduction (`compute_glacier_ice_albedo`).

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
geometry, CDS client error/retry classification and job-concurrency scheduling, the
glacier-decoupling lookup + adjustment against the vendored table, and the bare-ice albedo
reduction kernel / QC mask / QFLAG legend parsing).
Network-dependent tests are opt-in via environment variables so CI stays
offline by default:
- `CDS_API_KEY` — enables ERA5-Land Zarr integration tests
- `GEMB_TEST_INVARIANT=1` — enables ERA5-Land invariant NetCDF downloads
- `GEMB_TEST_CHUNK_MAP=1` — enables `climate_chunk_map` integration
- `GEMB_TEST_COPERNICUS_DEM=1` — enables Copernicus DEM `/vsicurl/` reads
- `GEMB_TEST_GEOID=1` — enables streaming the EGM96 geoid grid (`/vsicurl/`) for
  `geopotential2height`
- `GEMB_TEST_SATELLITE_ALBEDO=1` (plus `CDS_API_KEY`) — enables `satellite_albedo` CDS
  job orders. Slow (jobs queue server-side for minutes) and requires the dataset licence
  to have been accepted in the CDS web UI.

There is no per-file test runner; run a single suite by editing `include`s in
`test/runtests.jl` or invoking the file directly with the package loaded.

**Standalone scripts** (not included by `runtests.jl`, each needs a CDS token):
```bash
# Benchmark the load hot path. --check gates numerics against a golden snapshot
# (test/.golden_forcing.jls, gitignored) — run this before/after any optimization.
julia --project=. -t auto test/benchmark_load.jl            # run + record baseline
julia --project=. -t auto test/benchmark_load.jl --check    # run + assert vs golden
julia --project=. -t auto test/benchmark_load.jl --cpu-only # no network

julia --project=. test/verify_chunk_alignment.jl  # confirm chunk-ID formula vs live store
julia --project=. test/test_caching.jl            # Zarr.CachingStore cold/warm timing
```

### Running Examples
```bash
# Complete ERA5-Land + GEMB simulation example
export CDS_API_KEY="your-token-here"
julia --project=. examples/era5_land_example.jl

# Test authentication only
julia --project=. examples/test_authentication.jl

# Glacier bare-ice albedo from C3S satellite albedo. SLOW on a cold cache (CDS job
# queue latency, 30–90 min/year) and needs the satellite-albedo licence accepted.
julia --project=. examples/glacier_ice_albedo_example.jl
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

4. **`src/utils.jl`** - Credentials + vapor-pressure conversions shared across all workflows
   - `get_cds_api_key()` — resolves the CDS token from `ENV["CDS_API_KEY"]`, else the
     `key:` line of `~/.cdsapirc`. Both the Zarr (Bearer) and Retrieve API (`PRIVATE-TOKEN`)
     paths take their token from here.
   - `dewpoint_to_vapor_pressure()`, `vapor_pressure_to_relative_humidity()`,
     `relative_humidity_to_vapor_pressure()` — Magnus/Buck formulas
   - Shared by **all three** forcing-adjustment functions (elevation, temperature/precip,
     glacier), which otherwise duplicate the same boilerplate four ways:
     - `_rebuild_forcing(stack, meta; changed_layers...)` — rebuild a forcing `DimStack`
       replacing only the named layers. `rebuild` preserves per-layer metadata, so
       broadcast results come back annotated without threading it by hand; unchanged
       variables are never enumerated at the call site.
     - `_adjust_longwave(LW, e, T, e′, T′)` — Konzelmann ε_cs recomputation preserving the
       diagnosed cloud increment Δε. Exact identity when `(e′,T′) == (e,T)`.
     - `_wrap_longitude(lon)` — → −180…180°E, for the geoid grid and the glacier table
       (both −180…180, unlike the 0–360 ERA5-Land grids). Handles any multiple wrap.

5. **`src/elevation_adjustment.jl`** - Elevation downscaling (workflow 2)
   - `climate_adjust_for_elevation(stack, Δz; lapse_rate, precip_scaling_method)` — per-variable
     correction of a forcing `DimStack` to a target elevation (Glover 1999 scheme). `Δz = 0`
     reproduces the input; physical-range validation re-runs on the result.
   - `empirical_lapse_rate()` fits a local rate from neighbouring grid cells (RACMO-style).
   - Exported climatological monthly tables: `GREENLAND_LAPSE_RATE`, `ARCTIC_LAPSE_RATE`,
     `ANTARCTICA_LAPSE_RATE` (K/km, Jan→Dec). See README for the full per-variable table and references.
   - Sibling **`src/climate_adjustment.jl`** holds the *direct* (non-elevation) perturbations,
     reusing this file's `saturation_vapor_pressure` / `konzelmann_clear_sky_emissivity` /
     `_SIGMA_SB` plus `utils.jl`'s `_rebuild_forcing` / `_adjust_longwave` rather than
     redefining any of them:
     - `temperature_adjust(stack, ΔT)` — uniform offset in K, propagated into `vapor_pressure`
       (constant RH) and `longwave_downward` (Konzelmann ε_cs, cloud increment Δε preserved).
       Deliberately leaves `pressure_air` alone (set by elevation/synoptics, not by warming).
     - `precipitation_adjust(stack, scaling)` — fractional multiplier, non-negative; nothing else
       depends on the precipitation rate.
     - Both are scalar-argument only (no monthly/per-step vectors, unlike `lapse_rate`), record
       **cumulative** metadata so repeated calls compose — `"temperature_offset"` additive,
       `"precipitation_scaling"` multiplicative — and re-run `validate_climate_forcing_units`.
       The identity perturbation (`ΔT=0`, `scaling=1`) is exact. The two are independent by
       design: no implicit Clausius–Clapeyron precipitation response to ΔT.
   - Sibling **`src/glacier_decoupling.jl`** + **`src/glacier_adjustment.jl`** hold the
     *ambient → on-glacier* correction (Shaw et al. 2025). Reanalysis 2 m T is an ambient
     temperature: no glacier boundary layer, so a SEB model fed it overestimates melt.
     - `glacier_decoupling_table()` / `glacier_decoupling(rgi_id)` / `glacier_decoupling(lat, lon)`
       read the **vendored** `data/shaw2025_glacier_decoupling.csv.gz` (186,792 RGI v6 glaciers,
       gzipped CSV, no network). Regenerate with `data/make_shaw2025_decoupling.jl` from the
       Zenodo `.mat` — that script needs HDF5.jl, which is why **HDF5 is deliberately not a
       package dependency**. The paper's 5-predictor regression is *not* re-implemented:
       Table S2 has sign errors in `a4`/`a5` and a corrected refit still only reaches R²=0.54
       against the authors' own `k`. `k_lower` is the upstream `K_CI`, an absolute **lower
       bound** (not a half-width) and unclamped, so it can be negative.
     - `climate_adjust_for_glacier(stack, k)` / `(stack; rgi_id)` applies
       `T′ = T + (k−1)·max(T − T_ref, 0)` — increment form so `k=1` is bit-exact, **gated at
       the melting point** by default because a bare `k` multiplier *warms* sub-freezing air.
       Recomputes `longwave_downward` only. **Leaves `vapor_pressure` untouched on purpose** —
       do *not* "fix" this by reusing `temperature_adjust`'s constant-RH propagation: SM10
       Eq. 4 pivots `e` about 6.11 hPa, the opposite sign. Wind is likewise untouched.
     - **Order matters**: `climate_adjust_for_elevation` first, then this. **RGI regions 05 and
       19 are absent** from the table, so lookups fail there — pass `k` explicitly.

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

8. **`src/geoid.jl`** - Geopotential → height, and model surface elevation
   - `geopotential2height(Φ, lat, lon; height_reference)` — full conversion, *not* the crude
     `z / 9.80665`: latitude-dependent gravity + effective Earth radius (List 1968, as in NCL's
     `gp2gmh`) gives orthometric height, then the geoid undulation `N` is added for `:wgs84`
     (default) ellipsoidal height. `height_reference=:orthometric` skips the geoid entirely, so
     it needs no network. A `Raster` method broadcasts over the grid via `_axis_vector`, which
     reshapes lon/lat per their own dimension — correct for either axis order.
   - `geoid_undulation` / `GEOID_MODELS` — EGM96 (15′, default) and EGM2008 (2.5′) grids
     **streamed** from the PROJ CDN (`cdn.proj.org`) as COGs over `/vsicurl/`, sharing
     `_configure_gdal_http()` with the DEM loader. Nothing is bundled. Longitude −180…180°E.
   - `surface_elevation(model, lat, lon)` — elevation of the model grid cell nearest a point,
     from the `:z` invariant. Defaults to `:orthometric` (comparable to sea-level DEMs).
     Accepts either symbol spelling: `:era5land` (dataset) or `:era5_land` (invariant registry),
     normalized by `_invariant_model`. Maps negative longitude into the invariant's native
     0–360°E before the `Near` lookup.

9. **`src/cds_retrieve.jl`** - Generic CDS **Retrieve API** (job) client
   - `cds_retrieve(dataset, inputs, dest; token, ...)` — submit → poll → download. Plus
     `cds_estimate_cost` (`/costing`) and `cds_valid_options` (`/constraints`), both cheap
     and unqueued, so the live catalogue can be consulted instead of a hardcoded table.
   - **Authenticates with a `PRIVATE-TOKEN` header**, *not* the `Authorization: Bearer`
     scheme `AuthenticatedHTTPStore` uses for the ARCO Zarr stores. Same key from
     `get_cds_api_key()`, different header — do not conflate the two.
   - Dataset-agnostic on purpose, so other archived CDS products can reuse it. CDSAPI.jl was
     deliberately not adopted (credential-name mismatch, HTTP v2 requirement, collapsed 4xx
     handling); it is the fallback if the protocol changes.
   - **Safe to call concurrently.** CDS enforces a *per-dataset* queued-request limit and
     answers an over-limit job with HTTP 200 + `status: "rejected"` — the reason appears only
     on `/results`, not in the job body. `_cds_rejection_detail` fetches it and
     `CDSQueueLimitError` marks it retryable, so `cds_retrieve` resubmits with backoff
     (30 s → 300 s) rather than failing. Note `"rejected"` is deliberately *not* in
     `_CDS_PENDING_STATUSES`: it must terminate the poll loop so the retry decision is made
     once. `timeout` bounds retries plus waiting, so a saturated queue cannot spin forever.

10. **`src/datasets/copernicus_albedo.jl`** - C3S satellite surface albedo (workflow 5)
   - `satellite_albedo(; time_range, extent, variable, ...)` → lazy `RasterSeries` over `Ti`.
   - **This is not an invariant and not a lazy remote read.** The product has no ARCO Zarr
     copy, no COG bucket and no OPeNDAP endpoint, so data is *ordered* via async CDS jobs
     (minutes of latency), cached per timestep on disk, then opened `lazy=true`. Do not route
     it through `climate_model_invariant`: `_open_invariant_raster` deliberately drops `Ti`.
   - **Cost = variables × |years| × |months| × |nominal_days|, limit 20**, and `area` does
     *not* reduce it. `_albedo_chunk_timesteps` splits requests accordingly; the cross-product
     (not date-count) formula is verified against the live `/costing` endpoint — under-counting
     a month-straddling group gets the request server-rejected.
   - **One ordered *variable* is a multi-layer NetCDF**, not a single grid: `AL_{DH|BH}_{BB|NI|VI}`
     plus `_ERR` uncertainties, `QFLAG`, and a scalar `crs`. An unnamed `Raster(path; lazy=true)`
     silently picks `crs` (a 1-element `Char`, no X/Y), so `_albedo_resolve_layer` /
     `_albedo_open_layer` always name the layer; `_ALBEDO_DEFAULT_LAYERS` holds the broadband
     default and the `layer` kwarg overrides it. Multi-variable requests build a `RasterStack`
     keyed by *variable* symbol, not by internal layer name.
   - Timesteps are day 10 / day 20 / end-of-month (verified via `/constraints`). Sentinel-3
     `v3_1` 300 m era only (2018–2024); AVHRR/VGT/PROBA eras use different grids and are out
     of scope. Longitude is −180…180°E (like the DEM, unlike the ERA5-Land invariants).
   - Requires one-time licence acceptance in the CDS web UI or every request 403s.
   - **Jobs run concurrently, and that is the whole performance story.** Latency is queue
     time, not compute — tens of minutes per job, and only loosely related to job size (a
     3-timestep order has finished in ~30 min while another ran past 2 h). So
     `_run_concurrent_jobs` fans the chunks out over `@async` tasks (I/O-bound, no threads
     needed): `max_concurrent_jobs=6` in flight, submissions spaced `submit_stagger=15` s
     apart because a simultaneous burst trips the queue limiter even at counts that succeed
     when spread out. Measured: 6 concurrent one-month orders — 3 done at 27–32 min, the
     rest at ~79 min, all 6 in 83 min, vs ~3 h serially. **The gain saturates**: CDS runs
     only ~3 of an account's jobs at a time, so this is ~2×, not 6×, and raising the cap
     mostly courts rejections. **Do not "optimize" by splitting a request into per-month
     calls** — that re-serialises the ordering and is slower; pass the full range and let
     the chunker do it.
   - **A large `area` subset fails server-side; omitting `area` succeeds.** This is
     counter-intuitive and cost us hours, so measure before "optimizing" here. All for one
     timestep (2019-06-10, `albb_dh`): global `area` **failed** after 28 min, hemisphere band
     `[90,-180,30,180]` after 10 min, quarter-sphere `[90,-180,0,-90]` after 39 min, and even
     a Greenland box `[84,-73,59,-11]` after 37 min — each returning `status:"failed"` with an
     **empty traceback** on `/results`. The *same* request with **no `area` key at all`**
     succeeded in 29 min. Tiny boxes (~0.4°×0.8°) work, so the limit sits between those.
     **For continental-or-larger coverage, order globally and subset locally**; do not tile
     with `area`, and do not retry a failed wide box.
   - **A job's `timeout` is not what makes big orders fail.** `_ALBEDO_JOB_TIMEOUT` was raised
     3 h → **24 h** because real orders routinely run past 3 h server-side, and giving up
     happens *before* anything is cached, so the abandoned queue time is pure loss. But
     waiting longer does not rescue a too-large `area` — those die in well under an hour.
   - **Killing a run abandons its server-side jobs permanently.** Nothing persists job IDs, so
     an interrupted call leaves `running` jobs consuming the account's ~3 concurrent slots
     until they are dismissed by hand (`DELETE /jobs/{id}`; `GET /jobs?limit=N` lists them).
     A follow-on run then queues behind ghosts of the one that was killed.
   - Sibling **`src/glacier_ice_albedo.jl`** reduces that record to GEMB's **bare-ice albedo**:
     `compute_glacier_ice_albedo(years; extent, ...)` → `RasterStack` over (`X`, `Y`, `Ti`) with
     `:glacier_ice_albedo` (Float32, NaN where unresolved) and `:n_valid_observations`.
     Per pixel-year it averages the darkest `percentile` (default 5 %) of valid broadband
     retrievals — on a glacier the annual albedo *minimum* is the bare-ice state, and averaging
     a low percentile rather than taking the single minimum stops one bad retrieval from setting
     the answer (~36 obs/yr, so 5 % is the darkest 1–2).
   - Loaded **after** `datasets/copernicus_albedo.jl` in `GEMB_ClimateForcing.jl` — it calls
     `satellite_albedo`, aliases its `_ALBEDO_YEARS`, and reads the backing product file out of
     each layer's `"file"` metadata key to find the QFLAG legend. It does **not** re-derive the
     private cache layout; if you need the file behind a layer, read that key.
   - **One `satellite_albedo` call per year, with `layer = [layer, :QFLAG, err_layer]`.** All
     three live in the same product NetCDF, so the extra layers cost no extra CDS job, and one
     call keeps the timesteps aligned by construction instead of by matching positional indices
     across three separate series. Every loader knob (`timeout`, `max_concurrent_jobs`,
     `poll_interval`, `force_download`) is forwarded via `albedo_kwargs...` rather than restated
     with its own default.
   - **The annual reduction streams** (`_LowPercentileTopK` / `_accumulate!` / `_finalize`):
     only the `ceil(percentile · n)` darkest values per pixel are retained, so a year holds ~3
     grids of state instead of all ~36 masked grids — 49 → 5.5 MiB at 600×600, and the
     saving grows with the extent, which is the point (a Greenland-scale extent is where this
     otherwise OOMs). Kernel *time* is unchanged (~45 ms either way); this is a memory
     change, not a speed one. Bit-identical to sorting every valid value — the same `k` values
     reach the same `mean` in the same order — and `_low_percentile_mean` remains as the batch
     wrapper the tests compare against.
   - **Four keywords exist only for grids far larger than a glacier basin**, all off or no-ops
     by default so regional behaviour is untouched. Measured on a real global product file:
     the grid is **120960 × 47040** (lon −180…180, lat **80°N…−60°S**, spacing 1/336°) =
     5.69 Gpx at **1.9 B/px on disk**, so **10.9 GB/timestep and 392 GB/year**. (A tiny
     area-subset file measures 9.2 B/px — 5× worse, because compression works far better on
     a global field. Do not size a global run from a subset.) One global timestep folds in
     **198 s at 21.5 GB peak RSS**, so a year is ~2 h of compute plus CDS queue time:
     - `block_rows=512` reads each timestep in Y-slabs instead of whole (a global timestep's
       three layers are ~110 GB read at once). **Bit-identical at any block size** — verified
       against the whole-grid path on real files at block_rows ∈ {1,7,13,17,30,512}.
     - `scratch_dir` / `scratch_threshold_pixels=200_000_000`: past that pixel count the
       accumulator (`kmax·npx·4 + 2·npx·4` bytes — ~117 GB globally, vs 103 GB of RAM) is
       `Mmap`-backed on disk. Also bit-identical to the in-memory path.
     - `batch_timesteps` + `discard_after_fold`: order a year in batches and **delete each
       batch's product files once folded**, bounding peak disk to one batch — required
       globally, where a year's 392 GB does not fit alongside 91 GB of scratch.
       `batch_timesteps` is a *deliberate pessimisation* — batches order serially, so
       wall-clock grows with the batch count, against the concurrent single call the default
       uses. Only reach for it when the year's files genuinely do not fit.
       `discard_after_fold` is destructive: it gives up the cache, so a re-run reorders from
       CDS. It must `GC.gc()` before `rm`, because NCDatasets holds each file open until its
       handles are collected.
     - `n_expected` comes from `_albedo_timesteps` over the *whole year*, not from the batch,
       so `kmax` is identical however the year is split. Batched folding is verified
       bit-identical to the single call on real cached files (1, 2 and 4 batches).
   - Three deliberate QC choices, each a trap to re-check before "fixing":
     - **`albedo_range`'s 0.3 floor is glaciological, not physical.** Exposed ice rarely goes
       below ~0.3 broadband, so darker pixels are usually rock/water/shadow/failed inversion.
       Heavily dust- or algae-darkened ablation zones need it lowered (~0.15) or the default
       clips the signal being measured.
     - **`snow_presence` is deliberately kept, not rejected.** A bright snow-covered timestep is
       discarded by the low percentile anyway; rejecting it up front only biases the sample count
       against `min_samples`. Do not add it to `GLACIER_ICE_ALBEDO_QFLAG_REJECT`.
     - **`brdf_warning_*` rejection is off by default** (`qflag_reject_brdf_warning=false`) — a
       warning, not a failure, and over snow/ice it fires often enough to thin the sample
       materially.
   - QFLAG bits are read from each file's own `flag_masks`/`flag_meanings` (`_qflag_table`),
     never hardcoded — the convention is version-specific — and an unreadable legend degrades to
     range-only QC with a warning rather than erroring.
   - Was formerly `examples/albedo_annual_low_percentile.jl`; that CLI script is gone, don't
     restore it. `examples/glacier_ice_albedo_example.jl` is its replacement: a thin driver
     over the package function (summary + NetCDF write), runnable on `include` with no
     `ARGS`/`PROGRAM_FILE` handling, keeping no analysis logic of its own.

11. **`src/simulate/simulate_climate_forcing.jl`** - Synthetic forcing (workflow 4)
   - `simulate_climate_forcing(set_id, time_step_hours=0)` — generates a full stochastic forcing
     `DimStack` from a named parameter set (`simulation_parameter_sets`, e.g. `"test_1"`), seeded
     RNG (`MersenneTwister`) for MATLAB-matching reproducibility.
   - Time↔decimal-year helpers `datetime2decyear` / `decyear2datetime` live here.
   - Climatological averaging (`forcing_climatology`) lives in **GEMB.jl**, not here — it
     operates on GEMB's `ClimateForcing` type (drops leap day 366 and partial years, averaging
     complete years into a one-year cycle). Convert a `DimStack` with `GEMB.ClimateForcing(ds)`
     first.

12. **`src/fit_climate/`** - Parameter fitting (workflow 4, inverse of simulate)
   - `fit_air_temperature`, `fit_precipitation`, `fit_longwave_irradiance_delta`,
     `fit_seasonal_daily_noise` — estimate `simulate_*` coefficients from observed series.
   - `simulate_coeffs_disp` prints a fitted coefficient NamedTuple as copy-pasteable Julia.
   - `varname2longname` maps variable symbols to human-readable names.

### MATLAB parity

Workflow-4 functions are direct ports of GEMB's MATLAB `simulate_*` / `fit_*` code and are
expected to reproduce its numerics. When editing them, preserve the RNG seeding order and the
exact coefficient math — divergence from MATLAB is a regression, not a cleanup opportunity.

### Deliberate performance choices

These look like things to tidy up. They are not — each was measured, and the "cleaner" form is
slower or breaks bit-exactness. Comments at each site record the same thing.

- **`konzelmann_clear_sky_emissivity` uses `sqrt(sqrt(sqrt(x)))`, not `x^(1/8)`.** Algebraically
  identical, ~9× faster (three hardware `sqrt`s vs a generic `pow`), differs by ≤1 ulp. This is
  the hot kernel of every longwave adjustment, called twice per time step.
- **`validate_climate_forcing_units` uses paired `minimum`/`maximum`, not `extrema`.** Base's
  `extrema` does not SIMD-vectorize, so two vectorized passes beat one scalar pass by ~23×
  (8.46 ms → 375 µs over a 32-yr hourly record). NaN propagation is identical. This function runs
  at the end of *every* adjustment, so it dominated their cost.
- **Do not fuse the adjustment broadcasts.** Both `_adjust_longwave` and
  `simulate_shortwave_irradiance` were tried as single fused kernels and both got *slower*
  (10.5 → 12.2 ms and 8.0 → 11.2 ms) — staged per-array broadcasts vectorize better than one long
  kernel, even though fusing cuts allocations ~4×.
- **`fit_longwave_irradiance_delta`'s M-step keeps its allocating broadcasts.** A single-pass
  accumulation loop is faster but reassociates the additions away from `sum`'s pairwise order,
  shifting the fitted coefficients ~1e-14 — a MATLAB-parity regression. `sum(f, indices)` is not
  a workaround: reducing lazily over indices blocks differently from reducing over an `Array`.
  The E-step is where the time goes and *is* optimized (the responsibility matrix is hoisted out
  of the iteration loop and row slicing removed: 112M → 1.6k allocations). Its `P[k] * (pdf)`
  grouping is also load-bearing for bit-exactness.
- **`PrecompileTools` workload in `src/GEMB_ClimateForcing.jl` is deliberately network-free** —
  no Zarr, CDS, or `/vsicurl/` — so precompilation stays offline and CI-safe. `climate_forcing`
  is therefore not covered. It costs ~1.4 s of precompile time and buys ~2.6 s of first-call
  latency.

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

Use `src/datasets/era5_land.jl` as a template. Note that `CONTRIBUTING.md` still describes
step 3 as returning a `GEMB.ClimateForcing` struct — that is stale; loaders return a
`DimStack` and the struct conversion lives in GEMB.jl.

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
  required by `climate_adjust_for_elevation` (and by GEMB.jl's `forcing_climatology`).
- **Rasters.jl** — invariant fields and the DEM are returned as lazy `Raster`/`RasterStack`;
  satellite albedo as a lazy `RasterSeries` over `Ti`
  (backends: `RastersNCDatasetsExt` via NCDatasets, `RastersArchGDALExt` via ArchGDAL).
  Note `RasterSeries` has **no metadata field** — `rebuild(series; metadata=…)` is accepted
  and silently discarded, so provenance must be attached to the layers instead.
- **GDAL** (via ArchGDAL) — `/vsicurl/` remote COG reads for the Copernicus DEM and the
  EGM96/EGM2008 geoid grids.
- Coordinate reference system: WGS84 (EPSG:4326) throughout.

### Grid-orientation gotcha
Both the ERA5-Land ARCO grid and the invariant NetCDFs use **0–359.9°E longitude** and
**descending latitude (90→−90°N)**. When cropping a `Raster`, `X = a .. b` / `Y = a .. b`
takes `min .. max` regardless of storage order — but verify axis order after any reprojection
and never assume ascending latitude.
