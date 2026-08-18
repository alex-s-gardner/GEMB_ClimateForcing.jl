# Example: derive glacier bare-ice albedo from the C3S 300 m satellite albedo record
#
# `compute_glacier_ice_albedo` reduces the 10-daily albedo record to, per pixel and
# calendar year, the mean of the darkest few percent of that year's valid retrievals.
# On a glacier the annual albedo minimum *is* the bare-ice state, so this is the
# observational substitute for GEMB's tuned bare-ice albedo constant.
#
# Requirements:
# - CDS API key from https://cds.climate.copernicus.eu/
#     export CDS_API_KEY="your-token-here"      (or a `key:` line in ~/.cdsapirc)
# - One-time licence acceptance for the dataset, or every request returns HTTP 403:
#     https://cds.climate.copernicus.eu/datasets/satellite-albedo?tab=download#manage-licences
#
# Run it:
#     julia --project=. examples/glacier_ice_albedo_example.jl
#
# or from the REPL, where you can re-run it with different settings without paying
# for the CDS orders again (they are cached on disk):
#     include("examples/glacier_ice_albedo_example.jl")
#     ice = run_example(2019:2020; extent = Extent(X = (-50.0, -48.0), Y = (66.0, 67.0)))
#
# Budget 30–90 minutes per *cold* year: the data is ordered from CDS via async jobs,
# and the wall time is queue latency rather than compute. Cached years are nearly free.

using GEMB_ClimateForcing
using Rasters
using Rasters.Extents
using Dates
using Statistics
using Printf

## ------------------------------------------------------------------ configuration

# Russell Glacier area, west Greenland — a small extent with a real ablation zone, so
# the run finishes in one CDS job round. The global grid is 120960 × 47040 pixels and
# ~10.9 GB *per timestep*, so a small extent is much cheaper.
#
# For global coverage use `extent = nothing`, NOT a −180…180 / −90…90 `Extent`: the
# latter becomes a CDS `area` subset, and large `area` subsets fail server-side. See the
# large-grid keywords in `compute_glacier_ice_albedo`.
const EXTENT = Extent(X = (-50.3, -49.5), Y = (66.9, 67.3))

# Sentinel-3 300 m era only (2018–2024 — `GLACIER_ICE_ALBEDO_YEARS`). Two years is
# enough to see interannual variability; add more for a climatology.
const YEARS = 2019:2020

# Durable cache for the ordered product files. The package default lives under
# `tempdir()`, which the OS may clear — and a lost cache means reordering hours of CDS
# jobs, so point this somewhere permanent.
const CACHE_PATH = joinpath(homedir(), "data", "CDS")

# Where the NetCDF result is written (gitignored).
const OUTPUT_DIR = joinpath(@__DIR__, "output")

## ------------------------------------------------------------------------- the run

"""
    run_example(years = YEARS; extent = EXTENT, cache_path = CACHE_PATH, kwargs...)

Compute bare-ice albedo for `years`, print a per-year summary, write a NetCDF, and
return the `RasterStack`. Extra keywords go to `compute_glacier_ice_albedo`.
"""
function run_example(years = YEARS; extent = EXTENT, cache_path = CACHE_PATH,
                     output_dir = OUTPUT_DIR, kwargs...)
    # Fail before any ordering if there is no usable token. `compute_glacier_ice_albedo`
    # would resolve this itself; checking here turns a mid-run failure into an
    # immediate, actionable one.
    try
        get_cds_api_key()
    catch
        error("""
        No CDS API key found. Get one from https://cds.climate.copernicus.eu/ and either
            export CDS_API_KEY="your-token-here"
        or put a `key:` line in ~/.cdsapirc. The satellite-albedo licence must also be
        accepted once in the CDS web UI, or every request returns HTTP 403.
        """)
    end

    println("="^64)
    println("Glacier bare-ice albedo from C3S 300 m satellite albedo")
    println("="^64)
    @printf("  years:  %d–%d\n", first(years), last(years))
    @printf("  extent: X %.2f…%.2f   Y %.2f…%.2f\n",
            extent.X[1], extent.X[2], extent.Y[1], extent.Y[2])
    println("  cache:  ", cache_path)
    println("\nOrdering from CDS if not already cached — a cold year takes 30–90 min.\n")

    ice = compute_glacier_ice_albedo(years; extent, cache_path, kwargs...)

    summarize(ice)
    path = write_result(ice, years, output_dir)
    println("\nWrote: ", path)
    println("  layers: glacier_ice_albedo, n_valid_observations")

    return ice
end

"""
    summarize(ice)

Print per-year statistics over the resolved pixels of a `compute_glacier_ice_albedo`
result, plus the multi-year mean that GEMB would actually consume.
"""
function summarize(ice)
    albedo = ice[:glacier_ice_albedo]
    counts = ice[:n_valid_observations]
    npixel = length(view(albedo, Ti(1)))

    println("\n", "="^64)
    println("Bare-ice albedo by year")
    println("="^64)
    for (i, t) in enumerate(lookup(albedo, Ti))
        # NaN marks a pixel-year that failed QC or fell below `min_samples`, so it must
        # be filtered out rather than fed to `mean` — see the `n_valid_observations`
        # layer for *why* a pixel is missing.
        vals = filter(!isnan, vec(view(albedo, Ti(i))))
        obs = filter(>(0), vec(view(counts, Ti(i))))
        if isempty(vals)
            @printf("%d: no pixel resolved (of %d) — try a lower `albedo_range` floor or `min_samples`\n",
                    year(t), npixel)
        else
            @printf("%d: mean %.3f  min %.3f  max %.3f   %d/%d px resolved, median %d obs/px\n",
                    year(t), mean(vals), minimum(vals), maximum(vals),
                    length(vals), npixel, isempty(obs) ? 0 : round(Int, median(obs)))
        end
    end

    # One value per pixel, averaged over years — this is the form GEMB wants for its
    # bare-ice albedo field.
    multiyear = map(eachslice(albedo; dims=(X, Y))) do px
        v = filter(!isnan, collect(px))
        isempty(v) ? NaN32 : mean(v)
    end
    resolved = filter(!isnan, vec(multiyear))
    if isempty(resolved)
        println("\nMulti-year mean: no resolved pixel.")
    else
        @printf("\nMulti-year mean bare-ice albedo: %.3f  (%d of %d px, range %.3f–%.3f)\n",
                mean(resolved), length(resolved), npixel,
                minimum(resolved), maximum(resolved))
        println("Use this as GEMB's bare-ice albedo instead of a tuned regional constant.")
    end

    return multiyear
end

"""
    write_result(ice, years, output_dir) -> String

Write the stack to NetCDF and return the path.

NetCDF rather than GeoTIFF on purpose: GDAL cannot write a third dimension with a
`Date` lookup ("no valid permutation of dimensions"), whereas NCDatasets stores the
`Ti` dimension natively.
"""
function write_result(ice, years, output_dir)
    mkpath(output_dir)
    path = joinpath(output_dir,
                    "glacier_ice_albedo_$(first(years))_$(last(years)).nc")
    Rasters.write(path, ice; force=true)
    return path
end

## ------------------------------------------------------------------------- execute

# Runs on `include`, and `run_example` stays callable afterwards with other settings.
ice = run_example()

## ---------------------------------------------------------------------- variations
#
# Heavily dust- or algae-darkened ablation zones: the default 0.3 floor is a
# *glaciological* bound, not a physical one, and will clip the very signal being
# measured. Lower it and the darkest retrievals survive:
#
#     ice = run_example(2019; albedo_range = (0.15, 1.0))
#
# A sparser record (high latitude, persistent cloud) resolves more pixels if the
# minimum sample count is relaxed — at the cost of noisier values:
#
#     ice = run_example(2019; min_samples = 5, percentile = 0.10)
#
# Loader knobs pass straight through to `satellite_albedo`:
#
#     ice = run_example(2019; max_concurrent_jobs = 3, force_download = true)
#
# Inspect why pixels are missing — `n_valid_observations` is reported even where the
# albedo is NaN, precisely so this is answerable after the fact:
#
#     using Statistics
#     n = ice[:n_valid_observations]
#     count(iszero, n), count(>(0), n)
