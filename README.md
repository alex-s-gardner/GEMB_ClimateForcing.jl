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
- **Glacier bare-ice albedo** — observed bare-ice albedo per pixel-year, as the mean of each year's darkest few percent of albedo retrievals (`compute_glacier_ice_albedo`).

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/GEMB_ClimateForcing.jl")
```

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

The 30 m global Copernicus DEM is available via `model=:copernicus_dem_30m`, served from Cloud-Optimized GeoTIFFs with byte-range reads (tiles are never downloaded in full):

```julia
dem = climate_model_invariant(model=:copernicus_dem_30m,
                              extent=Extents.Extent(X=(-38.5, -38.0), Y=(72.5, 72.8)))
```

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

> **Longitude is −180…180°E** here (as for the Copernicus DEM), *not* the 0–360°E convention of the ERA5-Land invariants. Omitting `extent` requests global 300 m data — ~120960 × 47040 pixels per variable per timestep — so an extent is strongly recommended.

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

> **Budget 30–90 minutes per cold year** (~3 CDS jobs, ordered concurrently); cached years are nearly free. Pass a durable `cache_path`, since the default is under `tempdir()` and a lost cache means reordering everything. Request the full year range in one call — looping over months by hand re-serialises the ordering and is much slower.

The albedo, `QFLAG` and `_ERR` layers all live in the same product file, so all three are read in a single `satellite_albedo` call per year at no extra CDS cost. The annual reduction streams the timesteps, keeping only the darkest few values per pixel, so peak memory is set by the percentile rather than by the number of observations in the year. Loader keywords (`timeout`, `max_concurrent_jobs`, `force_download`, …) are forwarded to `satellite_albedo`.

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
