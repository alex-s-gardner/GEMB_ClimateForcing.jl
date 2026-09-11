# GEMB_ClimateForcing.jl

Load climate forcing data as a `DimStack` of climate variables. Converts seamlessly to [GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl) `ClimateForcing` via a package extension.

Pure Julia — no Python. Reanalysis data is read from Analysis-Ready, Cloud-Optimized (ARCO) Zarr stores over authenticated HTTPS, downloading only the requested time range and location (lazy loading).

## Capabilities

- **Reanalysis loading** — point time-series extraction from ERA5-Land ARCO Zarr stores (`climate_forcing`).
- **Elevation downscaling** — physically-based per-variable correction of a forcing stack to a target elevation for snow/ice surfaces (`climate_adjust_for_elevation`).
- **Climate perturbations** — uniform temperature offsets and fractional precipitation scaling for sensitivity/scenario experiments (`temperature_adjust`, `precipitation_adjust`).
- **On-glacier correction** — ambient → on-glacier air temperature via the Shaw et al. (2025) per-glacier decoupling factor (`climate_adjust_for_glacier`, `glacier_decoupling`).
- **Invariant fields** — lazy `Raster`s for land–sea mask, geopotential/orography, vegetation/soil/lake, and a 30 m global DEM (`climate_model_invariant`).
- **Chunk mapping** — visualize Zarr download locality before batch queries (`climate_chunk_map`).
- **Satellite albedo** — 10-daily C3S surface albedo (Sentinel-3, 300 m) as a lazy `RasterSeries`, ordered from the CDS Retrieve API (`satellite_albedo`).
- **Glacier bare-ice albedo** — observed bare-ice albedo as the mean of the darkest few percent of albedo retrievals: per pixel-year from C3S (`compute_glacier_ice_albedo`), the same statistic at a point list from MODIS MCD43A3 500 m in black-sky and white-sky forms (`compute_glacier_ice_albedo_modis`), or **one value per cell pooled over the whole 2000–2025 record** (`pool_ice_albedo_from_cache`), sampled at any GeoInterface point, line or polygon with `bare_ice_albedo`.
- **Global glacier cell list** — RGI 7.0 outlines rasterized onto the native MCD43A3 500 m grid: 3.36 M cells over 274,531 glaciers and 103 MODIS tiles, vendored and offline (`rgi7_modis_cells`).

## Installation

Requires **Julia 1.11 or newer**.

```julia
using Pkg
Pkg.develop(path="/path/to/GEMB_ClimateForcing.jl")
Pkg.instantiate()   # resolves EarthData.jl from git — see the note below
```

> **Note — `EarthData.jl` is pinned to `main`.** NASA CMR granule discovery is delegated to
> [`EarthData.jl`](https://github.com/evetion/EarthData.jl) (JuliaGeo) rather than
> reimplemented here, but the registered v0.1.0 predates the current CMR response schema and
> cannot parse a live query. `Project.toml` therefore pins the development branch through a
> `[sources]` block, which is what requires Julia 1.11 rather than the 1.10 LTS. The pin also
> means this package cannot be registered in General until upstream tags a release; missing
> functionality is being contributed upstream rather than kept local. `Pkg.instantiate()`
> fetches the pinned revision automatically — nothing extra to install.

## Download caches

Every product that downloads caches under one root, so a re-run reuses whatever is already on
disk. The root is `ENV["GEMB_CACHE_PATH"]` if set, otherwise this repository's own `data/`
directory, and each product appends its own gitignored subdirectory:

```
data/MCD43A3.061/         # MODIS granules and per-date samples
data/satellite_albedo/    # C3S ordered timesteps
data/invariant/<model>/   # ERA5-Land NetCDFs, Copernicus DEM tiles
```

```bash
export GEMB_CACHE_PATH=/big/volume/gemb_cache   # to put them somewhere else
```

**These caches reach hundreds of gigabytes** — the MCD43A3 per-date sample cache alone is 184 GB
— so set `GEMB_CACHE_PATH` if the volume holding the checkout is not sized for the product you are
building. Avoid a home directory, which is commonly a small SSD with a quota. Nothing here writes
to `$HOME`; the only paths derived from it are credential *reads* (`~/.cdsapirc`, `~/.netrc`,
`~/.edl_token`). Every function also takes an explicit `cache_path` that overrides the root, and a
`Pkg.add`-installed copy (read-only source tree) requires `GEMB_CACHE_PATH` and says so.

Losing the MCD43A3 sample cache is the expensive case: it turns a ~27 minute re-fold into a ~1 TB
re-download.

## Quick Start

ERA5-Land requires a free [CDS API key](https://cds.climate.copernicus.eu/). Register, copy your key, and set it in the environment:

```bash
export CDS_API_KEY="your-token-here"
```

```julia
using GEMB_ClimateForcing
using GEMB  # extension provides DimStack → ClimateForcing conversion
using Dates

# Load ERA5-Land for Summit, Greenland (returns a DimStack)
forcing_data = climate_forcing(
    :era5land, 72.58, -38.46;
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"],
)

# Convert to GEMB.ClimateForcing and run GEMB
cf = GEMB.ClimateForcing(forcing_data)
mp = GEMB.ModelParameters(output_frequency=:daily)
profile = GEMB.initialize_profile(mp, cf)
output = GEMB.gemb(profile, cf, mp)
```

See `examples/era5_land_example.jl` for a complete workflow.

## API Reference

### `climate_forcing(dataset, lat, lon; kwargs...)`

Load climate forcing and return a `DimStack`.

- `dataset::Symbol` — `:era5land`
- `lat::Real` — latitude [-90, 90]
- `lon::Real` — longitude [-180, 180] or [0, 360]

Keyword arguments:

- `time_range::Tuple{DateTime,DateTime}` — required
- `token::Union{String,Nothing}` — CDS API key (required for ERA5-Land)
- `chunk_strategy::Symbol=:geo` — `:geo` for point time-series (recommended), `:time` for spatial maps
- `cache_path::Union{String,Nothing}=nothing` — persistent disk cache (`Zarr.CachingStore`)

Returns a `DimStack` with `temperature_air`, `pressure_air`, `vapor_pressure`, `wind_speed`, `precipitation`, `shortwave_downward`, and `longwave_downward`, plus location and observation-height metadata.

### `climate_adjust_for_elevation(stack, delta_elevation; kwargs...)`

Downscale a forcing `DimStack` (must have a `Ti` dimension) for the elevation difference between the reanalysis grid cell and a target point, using physically-based per-variable corrections for **snow/ice (glacier and ice-sheet) surfaces**.

- `delta_elevation::Real` — `z_target − z_reanalysis` in metres (positive = target above the grid cell)
- `lapse_rate=6.5` — near-surface temperature lapse rate in **K/km**. Accepts a scalar, a length-12 monthly vector (Jan→Dec), or a per-time-step vector. Region-specific monthly tables are exported: `GREENLAND_LAPSE_RATE` (Fausto 2009), `ARCTIC_LAPSE_RATE` (Gardner 2009), `ANTARCTICA_LAPSE_RATE` (Fortuin & Oerlemans 1990). Use `empirical_lapse_rate` to fit the rate from neighbouring grid cells.
- `precip_scaling_method=nothing` — `nothing` leaves precipitation unchanged (RACMO practice); `:clausius_clapeyron` scales by `eₛ(T′)/eₛ(T)` (Glover 1999, elevation-desert effect).

| Variable | Adjustment | Key reference |
|----------|-----------|---------------|
| `temperature_air` | lapse `T − (Γ/1000)·Δz` | Glover 1999; Fausto 2009; Gardner 2009 |
| `pressure_air` | hydrostatic `P·exp(−g·Δz/(R_d·T̄))` | Glover 1999; Noël 2018 |
| `vapor_pressure` | constant relative humidity, recomputed at `T′` | Glover 1999; Curry & Webster 1999 |
| `longwave_downward` | Konzelmann (1994) clear-sky emissivity, preserving cloud increment Δε | Konzelmann 1994; Fiddes & Gruber 2014 |
| `shortwave_downward` | unchanged | — |
| `precipitation` | unchanged, or Clausius–Clapeyron `×eₛ(T′)/eₛ(T)` if requested | Glover 1999 |
| `wind_speed` | unchanged | — |

`Δz = 0` reproduces the input exactly. Physical-range validation re-runs on the result.

```julia
stack = climate_forcing(:era5land, 72.58, -38.46;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

# Downscale 250 m above the grid cell with Greenland monthly lapse rates
adjusted = climate_adjust_for_elevation(stack, 250.0; lapse_rate=GREENLAND_LAPSE_RATE)

# Fit the local gradient from neighbouring cells (RACMO-style)
Γ = empirical_lapse_rate(neighbour_T2m, neighbour_elevations)   # K/km
adjusted = climate_adjust_for_elevation(stack, 250.0; lapse_rate=Γ)
```

> The scheme follows Glover (1999, *J. Climate* 12, 551–563, Eqs. 15–20) and is the surface-field analogue of TopoSCALE (Fiddes & Gruber 2014); it is consistent with the RACMO2.3p2 studies of Noël et al. (2018, 2019). Near-surface lapse rates over melting ice are markedly shallower than the 6.5 K/km free-air default, so prefer a region-specific table or a locally-fitted rate.

### `temperature_adjust(stack, delta_temperature)`

Apply a **uniform temperature offset** (K) to a forcing `DimStack` and propagate it through the variables that depend on air temperature. Use for warming/cooling sensitivity experiments or to bias-correct a reanalysis against an observed temperature record. For a temperature change that arises from an elevation difference, use `climate_adjust_for_elevation` instead — it derives ΔT from a lapse rate and also corrects surface pressure.

| Variable | Adjustment |
|----------|-----------|
| `temperature_air` | `T + ΔT` |
| `vapor_pressure` | constant relative humidity, recomputed at `T′` (over-ice curve below 0 °C) |
| `longwave_downward` | Konzelmann (1994) clear-sky emissivity at `(e′, T′)`, preserving the cloud increment Δε |
| `pressure_air`, `wind_speed`, `precipitation`, `shortwave_downward` | unchanged |

`ΔT = 0` reproduces the input exactly. Metadata records `delta_temperature` and a **cumulative** `temperature_offset`, so repeated calls compose. Physical-range validation re-runs on the result — a large negative ΔT that drives temperature below 180 K (or longwave below 50 W/m²) raises an `ArgumentError` by design.

### `precipitation_adjust(stack, scaling)`

Rescale precipitation by a dimensionless fractional factor (`0.85` = 15 % drier, `1.15` = 15 % wetter), leaving every other variable unchanged. The scaling is uniform in time: it changes totals and event amplitude, not timing or intermittency. Precipitation phase is determined downstream from air temperature, so the scaling applies to whichever phase the temperature implies. `scaling` must be non-negative; `1.0` reproduces the input exactly.

Metadata records a **cumulative, multiplicative** `precipitation_scaling` (two calls of `1.1` record `1.21`) and a recomputed annual `precipitation_mean`.

The two perturbations are independent and commute — a temperature offset applies no implicit Clausius–Clapeyron precipitation response, so the accumulation change stays an explicit choice:

```julia
stack = climate_forcing(:era5land, 72.58, -38.46;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

warmer = temperature_adjust(stack, 2.0)              # +2 K scenario
wetter = precipitation_adjust(stack, 1.15)           # 15% more precipitation

# Warm and dry, on top of an elevation downscaling
scenario = precipitation_adjust(
    temperature_adjust(climate_adjust_for_elevation(stack, 250.0), 3.0), 0.85)
```

### `climate_adjust_for_glacier(stack, k)` / `climate_adjust_for_glacier(stack; rgi_id)`

Correct an **ambient** (off-glacier) forcing to **on-glacier** conditions. A melting surface pinned at 0 °C cools the air above it, and the resulting stable layer suppresses turbulent mixing and drives a katabatic wind. Reanalysis 2 m temperature carries none of this — ERA5-Land grid cells dwarf a valley glacier — so feeding it straight to a surface energy balance model overestimates melt (22 % of the mass-balance change per +1 °C; Greuell & Böhm 1998).

Applies the Shaw et al. (2025) decoupling factor `k` from the published **per-glacier lookup table** (186,792 RGI v6 glaciers, vendored in `data/`, no network needed). The paper's five-predictor regression is deliberately not re-implemented: its Table S2 has sign errors in `a4`/`a5`, and even a corrected refit reaches only R² = 0.54 against the authors' own published `k`. See [`docs/on_glacier_temperature_correction.md`](docs/on_glacier_temperature_correction.md).

| Variable | Adjustment |
|----------|-----------|
| `temperature_air` | `T + (k−1)·max(T − T_ref, 0)`, i.e. cooling proportional to how far ambient sits above melting |
| `longwave_downward` | Konzelmann (1994) clear-sky emissivity at the cooled `T′`, unchanged `e`, preserving the cloud increment Δε |
| `vapor_pressure`, `pressure_air`, `wind_speed`, `precipitation`, `shortwave_downward` | unchanged |

**Run `climate_adjust_for_elevation` first.** `k` multiplies an ambient temperature *at the glacier's elevation*; the reverse order is wrong.

`k` must be in `(0, 1]`; `1.0` reproduces the input exactly. Metadata records a **cumulative, multiplicative** `glacier_decoupling_factor` plus the lookup provenance (`glacier_decoupling_rgi_id`, `k_lower`, match distance).

```julia
stack = climate_forcing(:era5land, 45.97, 7.53;
                        time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
                        token=ENV["CDS_API_KEY"])

# Elevation correction to the glacier surface first, glacier correction second.
at_glacier = climate_adjust_for_elevation(stack, 2900.0 - metadata(stack)["elevation"])

on_glacier = climate_adjust_for_glacier(at_glacier; rgi_id="RGI60-11.02810")  # by RGI id
on_glacier = climate_adjust_for_glacier(at_glacier)                            # nearest centroid
on_glacier = climate_adjust_for_glacier(at_glacier, 0.83)                      # k directly

# Per-glacier k and its lower CI bound, for a sensitivity range.
row = glacier_decoupling("RGI60-11.02810")   # (; k, k_lower, lat, lon, ele, …)
```

Two limits to know:

- **`vapor_pressure` is left untouched by design.** This is *not* `temperature_adjust`'s constant-RH propagation. Shea & Moore (2010, Eq. 4) show `e_gla` pivots about 6.11 hPa (saturation over ice at 0 °C), so the boundary layer *adds* moisture when ambient air is drier than that — the opposite sign to constant-RH scaling. The correct scheme needs a flowpath-length raster.
- **RGI regions 05 (Greenland periphery) and 19 (Antarctic) are absent** from the Shaw table, so the lookup forms fail there; pass `k` explicitly. Cooling is gated to steps above melting by default (`apply_below_freezing=false`), since a bare `k` multiplier would *warm* sub-freezing air, and the regression is ablation-season only.

### `climate_model_invariant(; model, parameter, ...)`

Load a climate model's **time-invariant** fields as **lazy** `Raster`s — not present in the time-series Zarr stores. ECMWF distributes ERA5-Land invariants as global NetCDF files; they are downloaded once, cached, and opened lazily (no data read until indexed/cropped/collected).

```julia
using GEMB_ClimateForcing, Rasters

lsm = climate_model_invariant(parameter=:lsm)          # land-sea mask (0–1)
iceland = read(lsm[X = 335 .. 347, Y = 63 .. 67])      # crop then read

z = climate_model_invariant(parameter=:z)              # geopotential (m² s⁻²)
orography = z ./ 9.80665                                # elevation in metres

inv = climate_model_invariant()                         # all params as a lazy RasterStack
```

Available ERA5-Land parameters (GRIB shortName): `:lsm`, `:z`, `:cl` (lake cover), `:dl` (lake depth), `:cvl`/`:cvh` (low/high vegetation cover), `:tvl`/`:tvh` (low/high vegetation type), `:slt` (soil type), `:glm` (glacier mask). See `ERA5_LAND_INVARIANT_PARAMETERS`.

The 30 m global Copernicus DEM is available via `model=:copernicus_dem_30m`, served from Cloud-Optimized GeoTIFFs with byte-range reads, so only the bytes a crop needs are fetched:

```julia
dem = climate_model_invariant(model=:copernicus_dem_30m,
                              extent=Extents.Extent(X=(-38.5, -38.0), Y=(72.5, 72.8)))
dem[X(Near(-38.2)), Y(Near(72.6))]      # point elevation, metres above the EGM2008 geoid
```

Pass `cache_tiles=true` to download each covering tile (19–40 MB) once and read locally
afterwards. GDAL has **no on-disk cache for `/vsicurl/`** — its block cache lives in the process —
so without this, repeated lookups across sessions re-request the same bytes. Tiles are fetched
concurrently (`max_concurrent_downloads=4`), and the call refuses to run if the resolved cache is
a temp directory. Leave it off for a one-off window or a continental extent.

> [!NOTE]
> `surface_elevation(:era5land, lat, lon)` is the other elevation entry point, but it returns the
> **reanalysis grid cell's** elevation (~9 km) from the geopotential invariant, not a point
> elevation. In steep terrain the two differ by hundreds of metres — which is why
> `climate_adjust_for_elevation` exists.

> **Grid convention.** ERA5-Land invariants use **0–359.9°E** longitude and **descending** latitude (90→−90°N). The `X = a .. b` / `Y = a .. b` selector takes `min .. max` regardless of axis order.

### `satellite_albedo(; time_range, extent, variable, ...)`

Observed surface albedo from the C3S [Surface albedo 10-daily gridded data](https://cds.climate.copernicus.eu/datasets/satellite-albedo) product (Sentinel-3 OLCI+SLSTR, 300 m, `v3_1`, 2018–2024), returned as a **lazy `RasterSeries`** over `Ti`.

```julia
using GEMB_ClimateForcing, Rasters, Dates

alb = satellite_albedo(;
    time_range = (DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
    extent = Extents.Extent(X = (-48.0, -47.5), Y = (66.5, 67.0)),
    variable = :albb_dh,
)

lookup(alb, Ti)      # 2019-06-10, 2019-06-20, 2019-06-30
read(alb[1])         # materialise the first timestep (albedo fraction, 0–1)
```

Each ordered variable arrives as a NetCDF holding **several layers**: the broadband (`_BB`), near-infrared (`_NI`) and visible (`_VI`) albedos, a `_ERR` uncertainty for each, and a `QFLAG` quality mask. The full-spectrum broadband layer (`AL_DH_BB` / `AL_BH_BB`) is read by default; select another with `layer=`, and list what a cached file holds with `satellite_albedo_layers(path)`:

```julia
vis = satellite_albedo(; time_range = (DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
                       extent = Extents.Extent(X = (-48.0, -47.5), Y = (66.5, 67.0)),
                       layer = :AL_DH_VI)
```

Variables (`SATELLITE_ALBEDO_VARIABLES`) combine two axes — `albb` **broadband** vs `alsp` **spectral** (per-band), and `_dh` **directional-hemispherical** (black-sky, direct illumination) vs `_bh` **bi-hemispherical** (white-sky, fully diffuse). For a surface energy-balance model such as GEMB, `:albb_dh` (the default) is usually what you want; true albedo lies between the black- and white-sky values according to the diffuse fraction.

Timesteps follow the 10-daily ("decadal") convention: **day 10, day 20, and the last day of each month** (28/29/30/31). Only timesteps inside `time_range` are returned.

> **This is an order-and-cache pipeline, not a lazy remote read.** Unlike every other source here, this product has no ARCO Zarr copy, no COG bucket, and no OPeNDAP endpoint — the CDS catalogue exposes it only through a job-based Retrieve API. Each call **submits jobs that queue server-side for minutes**, caches the returned NetCDFs on disk, and opens them lazily. Repeat calls for cached timesteps submit no job at all.

> **Requests are size-limited, so long series mean many jobs.** CDS charges one unit per (variable × year × month × nominal-day) combination against a hard limit of 20, and narrowing `extent` reduces the download volume but **not** the cost. Requests are split automatically: one variable for a year is ~3 jobs, five years ~15. Ordering a decade will take a long time.

> **Accept the licence first.** Visit [the dataset's download tab](https://cds.climate.copernicus.eu/datasets/satellite-albedo?tab=download#manage-licences) once and accept the product terms, or every request fails with HTTP 403. This is the most common first-run failure.

> **Longitude is −180…180°E** here (as for the Copernicus DEM), *not* the 0–360°E convention of the ERA5-Land invariants. Omitting `extent` requests global 300 m data — 120960 × 47040 pixels (lat 80°N…−60°S), ~10.9 GB per variable per timestep — so a *small* extent saves a great deal. A *large* extent, however, is worse than none: CDS fails `area` subsets above roughly a Greenland-sized box, while the same request with no `area` succeeds, so continental-scale work should order globally and subset locally.

### `compute_glacier_ice_albedo(years; extent, ...)`

Observed **glacier bare-ice albedo**, reduced from the `satellite_albedo` record: for each pixel and calendar year, the mean of the darkest 5 % of that year's valid broadband retrievals. Returns a `RasterStack` over (`X`, `Y`, `Ti`).

```julia
using GEMB_ClimateForcing, Rasters, Dates, Statistics

ice = compute_glacier_ice_albedo(2019:2020;
    extent = Extents.Extent(X = (-48.0, -47.5), Y = (66.5, 67.0)))

ice[:glacier_ice_albedo]        # (X, Y, Ti) Float32, NaN where too few observations
ice[:n_valid_observations]      # how many retrievals each pixel-year drew on

# One bare-ice albedo per pixel, for GEMB:
albedo_ice = map(eachslice(ice[:glacier_ice_albedo]; dims=(X, Y))) do px
    v = filter(!isnan, collect(px))
    isempty(v) ? NaN32 : mean(v)
end
```

On a glacier the annual albedo minimum *is* the bare-ice state — seasonal snow has ablated, exposing ice at its most darkened by dust, black carbon and algae. This is the observational substitute for GEMB's bare-ice albedo, which is otherwise a tuned regional constant (~0.35–0.55). Averaging a low percentile rather than taking the single annual minimum keeps one bad retrieval from setting the answer: the 10-daily product gives ~36 observations a year, so 5 % averages the darkest 1–2.

Quality control runs per observation before any statistic is formed — missing/`_FillValue` pixels (which is what removes sea, cloud and shadow), albedo outside `albedo_range`, rejected `QFLAG` classes (read from each file's own `flag_masks`/`flag_meanings`, not hardcoded), and retrievals whose `_ERR` uncertainty exceeds `max_error`. Pixel-years below `min_samples` surviving observations are left `NaN` rather than computed from a handful of cloudy scenes.

> **The `albedo_range` floor of 0.3 is glaciological, not physical.** Exposed glacier ice rarely falls below ~0.3 broadband, so darker retrievals over a glacier pixel are usually rock, water, shadow or a failed inversion. For heavily dust- or algae-darkened ablation zones, lower it (toward ~0.15) or the default will clip the very signal being measured.

> **`snow_presence` is deliberately *not* rejected** by the QFLAG filter — a snow-covered timestep is a bright observation that the low percentile discards on its own, and rejecting it up front would bias the sample count instead. The v3.1 QFLAG legend is not a cloud mask at all; see `GLACIER_ICE_ALBEDO_QFLAG_REJECT`.

> **Budget 30–90 minutes per cold year** (~3 CDS jobs, ordered concurrently); cached years are nearly free. Caches default under this repository's `data/` directory (set `GEMB_CACHE_PATH` to move them); a lost cache means reordering everything. Request the full year range in one call — looping over months by hand re-serialises the ordering and is much slower.

See `examples/glacier_ice_albedo_example.jl` for a runnable workflow — per-year summary, the multi-year mean GEMB consumes, and a NetCDF write. It runs on `include` and leaves `run_example(years; kwargs...)` callable for other settings.

The albedo, `QFLAG` and `_ERR` layers all live in the same product file, so all three are read in a single `satellite_albedo` call per year at no extra CDS cost. The annual reduction streams the timesteps, keeping only the darkest few values per pixel, so peak memory is set by the percentile rather than by the number of observations in the year. Loader keywords (`timeout`, `max_concurrent_jobs`, `force_download`, …) are forwarded to `satellite_albedo`.

### `compute_glacier_ice_albedo_modis(lat, lon, years; ...)`

The **same statistic from MODIS MCD43A3** (500 m, daily, 2000-02-16 → present), evaluated at a **list of points** rather than over a gridded extent, and returning **black-sky** and **white-sky** albedo as separate layers.

Requires a free [NASA Earthdata Login](https://urs.earthdata.nasa.gov/users/new) and a bearer token:

```bash
export EARTHDATA_TOKEN="your-token-here"     # or a single line in ~/.edl_token
```

Tokens last 60 days and **you may hold at most two** — a third request returns HTTP 403. `earthdata_token_from_netrc()` will mint one from `~/.netrc` credentials; call it yourself once rather than from a loop.

```julia
using GEMB_ClimateForcing, DimensionalData, Statistics

# Points on the Russell Glacier ablation zone. The first two are ~20 m apart, i.e. inside
# one 500 m MODIS cell.
lat = [67.0900, 67.0902, 67.0950]
lon = [-50.0500, -50.0500, -49.9000]

ice = compute_glacier_ice_albedo_modis(lat, lon, 2019:2020; doy_range = :melt_season)

ice[:albedo_bsa]                  # (point, Ti) Float32 — black-sky, NaN where unresolved
ice[:albedo_wsa]                  # white-sky
ice[:n_valid_observations_bsa]    # retrievals each point-year drew on
ice[:cell_id]                     # "h16v02_r0699_c0124" — equal ⇒ identical albedo
ice[:latitude], ice[:longitude]   # centre of the cell actually sampled

# One value per point, for GEMB:
albedo_ice = [mean(filter(!isnan, collect(ice[:albedo_bsa][p, :]))) for p in 1:length(lat)]
```

**Points are deduplicated to unique MODIS cells**, so several points inside one 500 m cell cost one sample and are populated from a single derivation — `cell_id` makes that auditable, and points sharing it necessarily share their albedo bit-for-bit. A vector of `(lat, lon)` tuples is accepted as an alternative first argument.

> **Download volume is the cost and it is not reducible.** MODIS offers no server-side subsetting: the unit of transfer is a whole ~70 MB granule, of which ~85 % is layers never read. Budget `n_tiles × n_dates × 70 MB` — a Greenland point list spanning 8 tiles over a full year is ~200 GB. **`doy_range` is the fix**: at high latitude polar night yields no usable retrieval, so a melt-season window removes ~2.5× of the download at approximately zero cost in surviving samples. `stride` is the blunt fallback.
>
> **Pass `doy_range = :melt_season` rather than a hardcoded window.** The melt season is a different half of the year in each hemisphere, so a northern window like `(180, 220)` samples austral *midwinter* in Patagonia or Antarctica and resolves nothing. `:melt_season` reads the sign of the points' latitudes and picks `(120, 290)` or `(300, 110)`; a list straddling the equator falls back to the whole year, so split it by hemisphere and call twice. An explicit tuple with `first > last` wraps New Year, which is how the southern window is expressed — that pools the tail of one melt season with the start of the next inside a calendar year, which is acceptable for a darkest-percentile statistic. Granules are deleted after each date is folded (`keep_granules = false`, the default — *inverted* from the CDS path, where a re-order costs hours of queue latency rather than bandwidth-bound minutes), so peak disk stays at one date's tiles. Per-date sampled cell values are cached, so a re-run with the same points never re-downloads.

Quality control mirrors the C3S path except in one respect: MCD43A3's quality band is a small-integer **class**, not a bitmask, so `qa_keep` is a **whitelist**. There is no per-pixel uncertainty layer and hence no `max_error`. The glaciological 0.3 `albedo_range` floor and the decision to keep bright snowy observations both transfer unchanged.

The class decomposes: its **low bit is inversion quality** (even = full BRDF inversion, odd = magnitude inversion) and its **upper bits are Band 5/6 detector health** (`÷2` → 0 = both fine, 1 = Band 6 fill, 2 = Band 5 fill, 3 = both), with `255` for fill.

| | inversion | Band 5 | Band 6 | | | inversion | Band 5 | Band 6 |
|---|---|---|---|---|---|---|---|---|
| `0` | full | ok | ok | | `4` | full | fill | ok |
| `1` | magnitude | ok | ok | | `5` | magnitude | fill | ok |
| `2` | full | ok | fill | | `6` | full | fill | fill |
| `3` | magnitude | ok | fill | | `7` | magnitude | fill | fill |

The default `qa_keep = [0, 2, 4, 6]` is therefore **every full BRDF inversion**: classes 2, 4 and 6 are full-quality retrievals that merely lack a shortwave-infrared *spectral* band, which does not matter for the broadband `shortwave` albedo reduced here — and Aqua's Band 6 has 15 of 20 detectors non-functional, so they are common.

> [!TIP]
> **The QA whitelist is usually the binding filter, not data quality.** Over 200 real Iceland cells in late June, class `1` (magnitude inversion) was 64–72 % of pixels against 13–14 % for class `0`. If `min_samples` is unreachable, adding the odd classes — `qa_keep = [0, 1, 2, 4, 6]`, or all of `0:7` — multiplies the usable sample several-fold. The trade is that a magnitude inversion scales an *a priori* BRDF shape to the observed magnitude, so the level stays observation-driven but the black-sky/white-sky split rests on an assumed shape.

Requesting more layers costs **no extra download** — the spectral bands live in the same granules:

```julia
ice = compute_glacier_ice_albedo_modis(lat, lon, 2019;
    layers = (:Albedo_BSA_shortwave, :Albedo_WSA_shortwave, :Albedo_BSA_vis))
ice[:albedo_bsa_vis]
```

See `examples/glacier_ice_albedo_modis_example.jl` for a runnable workflow.

#### Choosing a source

| | `compute_glacier_ice_albedo` (C3S) | `compute_glacier_ice_albedo_modis` (MCD43A3) |
|---|---|---|
| Resolution | 300 m | 500 m |
| Cadence | 10-daily (~36/yr) | daily (~365/yr) |
| Record | 2018–2024 (Sentinel-3 era) | 2000-02-16 → present |
| Interface | gridded `extent` | point list |
| Output | one broadband albedo | **black-sky + white-sky** |
| Quality band | `QFLAG` bitmask + `_ERR` uncertainty | quality *class* whitelist, no uncertainty |
| Access | CDS async jobs — hours of queue latency | HTTPS granule download — bandwidth-bound |
| Credentials | CDS API key + licence acceptance | Earthdata Login token |

Prefer MODIS for long records, per-point work, or when black-sky/white-sky are needed separately; prefer C3S for finer spatial detail over a contiguous area. The reduction statistic is identical, so the two are directly comparable.

### `rgi7_modis_cells()` — every glacier on Earth, on the MODIS grid

RGI 7.0 glacier outlines rasterized onto the **native** MCD43A3 500 m sinusoidal grid, vendored
in `data/` and read with no network access. This is the point list the global bare-ice albedo
product is evaluated at.

```julia
cells = rgi7_modis_cells()
# RGI7ModisCells(3375423 rows over 274531 RGI 7.0 glaciers, 103 tiles, 103530 forced)

north, south = rgi7_hemisphere_split(cells)          # exact: v ≤ 8 is the northern hemisphere
pts = rgi7_modis_unique_cells(cells; index = north)  # 2,584,232 distinct cells, 63 tiles

ice = compute_glacier_ice_albedo_modis(pts, 2019; doy_range = (120, 290), cell_id = false)

rgi7_glacier_cells(cells, "RGI2000-v7.0-G-06-00241")  # that glacier's cell indices
lat, lon = rgi7_modis_cell_points(cells)              # cell centres, recomputed not stored
```

| | |
|---|---|
| Glaciers / area | 274,531 / 706,744 km² (RGI 7.0 published totals) |
| MODIS tiles | **103** — 63 northern, 40 southern, disjoint |
| Distinct cells | **3,360,716** (99.4 % area closure before the forced-cell rule) |
| Rows | 3,375,423 `(cell, glacier)` pairs |
| Vendored size | 12 MB across two gzipped CSVs |

Rasterizing *into* the product's own grid rather than warping the product onto a lat/lon grid is
the point: resampling a single day mixes neighbouring pixels' albedo *before* the
darkest-percentile reduction, biasing the very tail being measured.

> [!IMPORTANT]
> **Filter on `forced` before computing any aggregate.** 37.7 % of RGI 7.0 by count (but only
> ~2 % by area) is smaller than one 0.2146 km² cell, so those glaciers are represented by the
> single cell containing their representative point. That cell is majority *not* glacier — rock,
> moraine and water are all darker than ice — so its bare-ice albedo is **biased low**. The flag
> is carried through into the output product for exactly this reason.

Two further consequences worth knowing:

- **Rows are unique as `(cell, glacier)` pairs, not as cells.** A cell carries more than one
  glacier only where a sub-cell glacier would otherwise have none (14,707 rows, 0.44 %) — two
  tiny glaciers in one 463 m pixel genuinely share that pixel's albedo. Use
  `rgi7_modis_unique_cells` for anything that *samples*, or the download is over-counted.
- **A contested cell goes to the smaller glacier.** Otherwise a one-cell glacier would be erased
  by a large neighbour, and since the reducer is `min` over an area rank it is commutative, so
  the vendored file does not depend on the order the 19 regional shapefiles were read in.

### `rgi7_ice_albedo_climatology(years; reduction = median, ...)`

Reduce the per-year bare-ice albedo files over many years, per MODIS cell. Each of those files
already holds one cell's *annual* darkest-percentile mean, so a climatology is a pure read — no
download, no re-fold, no GeoParquet dependency.

```julia
clim = rgi7_ice_albedo_climatology(2001:2024)                          # median annual albedo
clim = rgi7_ice_albedo_climatology(2001:2024; reduction = mean, min_years = 10)
clim = rgi7_ice_albedo_climatology(2001:2024; reduction = x -> quantile(x, 0.1))

clim[:albedo_bsa]     # reduced black-sky albedo per cell
clim[:n_years_bsa]    # years each cell actually contributed
```

`reduction` is a **function handle** applied per cell to that cell's valid annual values. It
never sees a `NaN` — unresolved cell-years are dropped first, which is what makes `:n_years`
meaningful — and a cell with fewer than `min_years` valid years is reported `NaN` rather than
reduced over almost nothing. Returns a `DimStack` over `Dim{:point}` aligned with
`rgi7_modis_unique_cells`.

> [!NOTE]
> **This reduces annual statistics; it does not pool observations.** `median` gives the median
> *of the annual bare-ice albedos*, weighting every year equally. That is normally what
> "bare-ice albedo climatology" means, and it is a **different quantity** from one darkest-5 %
> taken over all years pooled, which would be dominated by whichever years were darkest. The
> pooled version is what `pool_ice_albedo_from_cache` computes (see below), from the per-date raw
> samples that stay on disk under `<cache_path>/samples/` — a CPU re-fold, not a re-download.

### `pool_ice_albedo_from_cache(cells; ...)` — one value per cell, over the whole record

Year-to-year MCD43A3 coverage is far too variable for an annual bare-ice albedo to be comparable:
the fraction of glacier cells resolving in a *single* year runs from 21.6 % (2001) down to 6.5 %
(2025) in the north, and 0.6–6.2 % in the south. Pooling every valid retrieval in the record into
one darkest-5 % mean per cell fixes that:

| | per-year | pooled |
|---|---|---|
| cells resolved | 6.5–21.6 % N, 0.6–6.2 % S | **95.0 % N, 72.6 % S** |
| glaciers with a value | 30.8 % | **83.7 %** |
| RGI 7.0 area with a value | 84.2 % | **97.5 %** |
| area with `n_valid ≥ 30` | 61.3 % | **93.7 %** |

Area-weighted global bare-ice albedo is **0.451**; the whole product is 34.7 MB for 3.36 M cells.

```julia
cells = rgi7_modis_unique_cells(rgi7_modis_cells())
north = filter(c -> c[2] <= 8, cells)          # v <= 8 is exactly the northern hemisphere
ice = pool_ice_albedo_from_cache(north)

ice[:albedo_bsa]                  # one value per cell over 2000–2025
ice[:n_valid_observations_bsa]    # retrievals behind it
ice[:k_used_bsa]                  # how many of the darkest were averaged
```

**It downloads nothing.** The per-date sample cache stores raw digital numbers and the QA class
rather than post-QC values, so `percentile`, `albedo_range` and `qa_keep` are all free to change by
re-folding — only the date list, cell list and layer set are fixed at download time. A requested
date that is not cached is an **error**, not a gap: a pooled statistic over an unknown subset of
the record is not interpretable.

`POOLED_ICE_ALBEDO_RANGE` is `(0.25, 1.0)`, not the per-year path's `(0.3, 1.0)` — every per-year
file has `minimum == exactly 0.3000`, meaning that floor was clipping the dark tail it was supposed
to bound. **No sample-count threshold is applied at all**; judge a cell by `n_valid`.

> [!WARNING]
> Region 19 (Antarctic periphery) averages **0.721** against 0.31–0.46 elsewhere: those cells
> rarely expose bare ice even in their darkest 5 %, so that number is snow. Filtering on albedo
> alone silently mixes the two populations.

### `bare_ice_albedo(geometry, reducer = nothing; ...)` — sample it at a geometry

```julia
using GEMB_ClimateForcing, Statistics
import GeoInterface as GI

a = bare_ice_albedo(GI.Point(-29.0341, 69.9979))
a[:albedo_bsa][1], a[:n_valid_bsa][1], a[:k_used_bsa][1]     # 0.253, 469, 24

poly = GI.Polygon([GI.LinearRing([(-29.15, 69.94), (-28.91, 69.94),
                                  (-28.91, 70.06), (-29.15, 70.06), (-29.15, 69.94)])])
bare_ice_albedo(poly)              # DimStack over Dim{:cell}, one row per selected cell
bare_ice_albedo(poly, median)      # NamedTuple of scalars: albedo_bsa 0.253, over 42 of 550 cells
```

Any GeoInterface geometry in **(lon, lat) degrees** — point, multipoint, line, ring, polygon or
multipolygon. Selection is nearest-cell throughout: a point takes the cell containing it, a line
every cell it passes through, a polygon every cell whose centre it covers (or `boundary=:touches`
/ `:inside`). The geometry is reprojected into the MODIS sinusoidal grid and burned there, so the
value reported for a cell is the value computed for that cell — nothing is interpolated or
resampled.

Layers are `:albedo_bsa` / `:albedo_wsa` with their `:n_valid_*` and `:k_used_*` counts, plus
`:latitude`, `:longitude`, `:cell_id` and `:in_product`.

> [!IMPORTANT]
> `:in_product = false` means the cell is **off-glacier**, not merely cloudy — both read `NaN`.
> RGI 7.0's glacier product **excludes the ice sheets** (region 05 is Greenland *periphery*), so
> an interior ice-sheet point is legitimately absent rather than missing data.

The first call parses a hemisphere table (~1.7 s) and memoizes it; later queries are ~0.03 ms at
99 MiB resident per hemisphere. Albedo is **not clamped to 1** — MCD43A3's valid range reaches
32766, and the QC range is a glaciological filter, not a physical clamp. Values are unsuitable
as-is for GEMB's `albedo_ice`, which asserts `0.2 ≤ albedo_ice ≤ 0.6`.

Regenerate the tables with `data/make_rgi7_modis_cells.jl` (offline apart from a one-off ~422 MB
shapefile fetch); build the global albedo product with `data/run_rgi7_ice_albedo_modis.jl`
(~1 TB of download per year — measured 745 GB north + 261 GB south — resumable), pool it with
`data/run_rgi7_pooled_albedo.jl` (no download, ~27 min), and write GeoParquet with
`data/make_rgi7_ice_albedo_parquet.jl`. See `CLAUDE.md` items 10c and 10d for the full cost model
and the reasoning behind each burn rule.

## ERA5-Land Details

Data is read from ECMWF's ARCO Zarr stores at `arco.datastores.ecmwf.int`. Variables are loaded from four store groups and converted for GEMB:

| Store group | Raw variables | Derived output |
|---|---|---|
| `sfc-2m-temperature` | `t2m`, `d2m` | `temperature_air` (K), `vapor_pressure` (Pa via dewpoint) |
| `sfc-pressure-precipitation` | `sp`, `tp` | `pressure_air` (Pa), `precipitation` (kg/m² = `tp`×1000) |
| `sfc-wind` | `u10`, `v10` | `wind_speed` (m/s magnitude) |
| `sfc-radiation-heat` | `ssrd`, `strd` | `shortwave_downward`, `longwave_downward` (W/m² = J/m² ÷ 3600) |

Grid: 0.1° (~9 km), 1801 × 3600 (lat × lon), hourly, 1950–present.

**Chunk layout** (storage order `time × lat × lon`):

| Strategy | time chunk | lat chunk | lon chunk | Optimized for |
|---|---|---|---|---|
| `:geo` (default) | 33,792 | 4 | 8 | Point time-series (~3.85 years/chunk, tiny spatial footprint) |
| `:time` | 1 | 1,024 | 1,024 | Spatial maps (one timestep, ~103°×103° tile) |

`:geo` is strongly preferred for single-point simulations: a multi-decade series at one location reads only a handful of chunks, whereas `:time` would load 1024×1024 tiles to recover a single point.

**Performance.** First load ~10–25 s for a year of hourly data; parallel loading of the four groups gives a 1.5–2× speedup. With `cache_path`, subsequent loads (including in later Julia sessions) skip the network entirely.

## Citation

If you use ERA5-Land data, please cite:

> Muñoz Sabater, J. (2019): ERA5-Land hourly data from 1950 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: [10.24381/cds.e2161bac](https://doi.org/10.24381/cds.e2161bac)

If you use the satellite surface albedo data, please cite:

> Copernicus Climate Change Service, Climate Data Store (2019): Surface albedo 10-daily gridded data from 1981 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: [10.24381/cds.ea87ed30](https://doi.org/10.24381/cds.ea87ed30)

## License

MIT License — see LICENSE file for details.

## Related Projects

- [GEMB.jl](https://github.com/alex-s-gardner/GEMB.jl) — Glacier Energy and Mass Balance model
- [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land) — dataset information
