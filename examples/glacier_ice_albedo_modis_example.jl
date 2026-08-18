# Example: derive black-sky and white-sky glacier bare-ice albedo at a POINT LIST from the
# MODIS MCD43A3 500 m daily albedo product.
#
# `compute_glacier_ice_albedo_modis` applies the same statistic as its C3S sibling — per
# cell and calendar year, the mean of the darkest few percent of that year's valid
# retrievals — but from MODIS, which is longer (2000→present), denser (daily), and
# distributes black-sky (BSA, directional-hemispherical) and white-sky (WSA,
# bihemispherical) albedo as separate layers, which GEMB's radiation scheme can use
# directly.
#
# Requirements:
# - A NASA Earthdata Login account: https://urs.earthdata.nasa.gov/users/new
# - An EDL bearer token in the environment:
#     export EARTHDATA_TOKEN="your-token-here"     (or a single line in ~/.edl_token)
#   Mint one at https://urs.earthdata.nasa.gov/profile/edit/user_tokens. Tokens last 60
#   days and you may hold at most TWO — a third request returns HTTP 403. If you keep
#   credentials in ~/.netrc, `earthdata_token_from_netrc()` will fetch one for you; call it
#   yourself once, not in a loop, for the same two-token reason.
#
# Run it:
#     julia --project=. examples/glacier_ice_albedo_modis_example.jl
#
# or from the REPL, where a re-run with different QC settings is nearly free because the
# per-date sampled cell values are cached:
#     include("examples/glacier_ice_albedo_modis_example.jl")
#     ice = run_example(2019:2020; albedo_range = (0.15, 1.0))
#
# ⚠️  DOWNLOAD VOLUME IS THE COST, AND IT IS NOT REDUCIBLE.
# There is no server-side subsetting: the unit of transfer is a whole ~70 MB granule, of
# which ~85 % is layers never read. Budget roughly
#     n_tiles × n_dates × 70 MB
# so this example (1 tile × 41 days) is ~3 GB, while a Greenland point list spanning 8
# tiles over a full year is ~200 GB. `doy_range` is the real fix — at high latitude polar
# night yields no usable retrieval, so restricting to the melt season removes ~2.5× of the
# download at approximately zero cost in surviving samples. Granules are deleted after each
# date is folded (`keep_granules=false`, the default), so peak disk stays at one date.

using GEMB_ClimateForcing
using DimensionalData
using Dates
using Statistics
using Printf
# NCDatasets directly, not `Rasters.write`: the result is a point list, so it has no X/Y
# dimensions for the Rasters writer to key on.
using NCDatasets

## ------------------------------------------------------------------ configuration

# Points on and around the Russell Glacier ablation zone, west Greenland. The first two are
# deliberately ~20 m apart, i.e. inside a single 500 m MODIS cell: they demonstrate the
# deduplication — the value is derived ONCE for that cell and both points report it.
const POINTS = [
    (67.0900, -50.0500),   # ablation zone
    (67.0902, -50.0500),   # ~20 m north of the above — SAME MODIS cell
    (67.0950, -49.9000),   # a few km east
    (67.1500, -50.2000),   # up-glacier, closer to the equilibrium line
    (66.9500, -50.1000),   # south margin
]

# 2000-02-16 → present (`MCD43A3_YEARS`).
const YEARS = 2019:2020

# Melt season only. THE volume knob — see the warning above.
#
# `:melt_season` derives the window from the sign of the points' latitudes, which is the only
# safe default for a package that gets pointed at both ice sheets: the southern melt season
# runs Nov→Apr, so a hardcoded northern window like `(180, 220)` samples the *coldest* part
# of a Patagonian or Antarctic year and resolves nothing at all. An explicit tuple still
# works, and `first > last` wraps New Year:
#
#     doy_range = (180, 220)    # NH high summer — narrow, ~3 GB per tile here
#     doy_range = (300, 110)    # SH melt season, wrapping
#     doy_range = nothing       # whole year
const DOY_RANGE = :melt_season

# Durable cache for the sampled cell values (a few kB per date). The package default lives
# under `tempdir()`, which the OS may clear; a lost cache means re-downloading granules.
const CACHE_PATH = joinpath(homedir(), "data", "MODIS")

# Where the NetCDF result is written (gitignored).
const OUTPUT_DIR = joinpath(@__DIR__, "output")

## ------------------------------------------------------------------------- the run

"""
    run_example(years = YEARS; points = POINTS, doy_range = DOY_RANGE, kwargs...)

Compute black-sky and white-sky bare-ice albedo at `points` for `years`, print a summary,
write a NetCDF, and return the `DimStack`. Extra keywords go to
`compute_glacier_ice_albedo_modis`.
"""
function run_example(years = YEARS; points = POINTS, doy_range = DOY_RANGE,
                     cache_path = CACHE_PATH, output_dir = OUTPUT_DIR, kwargs...)
    # Fail before any download if there is no usable token. The driver would resolve this
    # itself; checking here turns a mid-run failure into an immediate, actionable one.
    try
        get_earthdata_token()
    catch
        error("""
        No Earthdata token found. Create an account at
            https://urs.earthdata.nasa.gov/users/new
        mint a token at
            https://urs.earthdata.nasa.gov/profile/edit/user_tokens
        and either `export EARTHDATA_TOKEN="..."` or put it on a line in ~/.edl_token.
        If you have ~/.netrc credentials, `earthdata_token_from_netrc()` fetches one.
        Note the limit of two concurrent tokens — a third request returns HTTP 403.
        """)
    end

    # Ask the package for the resolved window and its dates rather than differencing the
    # tuple: `:melt_season` has no bounds until the latitudes are known, and a wrapping
    # southern window has `last < first`, so `last - first + 1` would be negative.
    resolved = GEMB_ClimateForcing._resolve_doy_range(doy_range, first.(points))
    ndates = sum(length(GEMB_ClimateForcing._modis_dates(y, resolved, 1)) for y in years)

    println("="^70)
    println("Glacier bare-ice albedo (black-sky + white-sky) from MODIS MCD43A3 500 m")
    println("="^70)
    @printf("  points:    %d\n", length(points))
    @printf("  years:     %d–%d\n", first(years), last(years))
    @printf("  doy_range: %-12s (%s, %d dates)\n",
            isnothing(resolved) ? "all" : string(resolved[1], "–", resolved[2]),
            doy_range === :melt_season ?
                (all(>=(0), first.(points)) ? "northern melt season" :
                 all(<(0), first.(points)) ? "southern melt season" :
                 "points straddle the equator — whole year") : "explicit",
            ndates)
    println("  cache:     ", cache_path)
    @printf("\n~%.1f GB will be downloaded per tile touched and then discarded.\n", ndates * 70 / 1024)
    println("Downloads are bandwidth-bound, not queued — no CDS-style job latency here.\n")

    ice = compute_glacier_ice_albedo_modis(points, years; doy_range, cache_path, kwargs...)

    summarize(ice, points)
    path = write_result(ice, years, output_dir)
    println("\nWrote: ", path)

    return ice
end

"""
    summarize(ice, points)

Print the per-point, per-year black-sky and white-sky albedo, then the multi-year mean per
point that GEMB would actually consume.
"""
function summarize(ice, points)
    bsa = ice[:albedo_bsa]
    wsa = ice[:albedo_wsa]
    nobs = ice[:n_valid_observations_bsa]
    years = year.(lookup(ice, Ti))

    println("\n", "="^70)
    println("Albedo by point and year")
    println("="^70)
    println("  point    requested lat/lon       sampled cell         year   BSA    WSA   obs")
    for p in eachindex(points)
        for (i, y) in enumerate(years)
            # NaN marks a point-year that failed QC or fell below `min_samples`. The
            # observation count is reported either way, precisely so the two are
            # distinguishable: zero obs means QC rejected everything, a small nonzero count
            # means `min_samples` bit.
            b, w = bsa[p, i], wsa[p, i]
            @printf("  %5d    %7.4f, %8.4f   %-19s  %4d  %s  %s  %4d\n",
                    p, points[p][1], points[p][2], ice[:cell_id][p], y,
                    isnan(b) ? "  -  " : @sprintf("%.3f", b),
                    isnan(w) ? "  -  " : @sprintf("%.3f", w),
                    nobs[p, i])
        end
    end

    # The dedup, made visible: points sharing a cell_id necessarily share their values,
    # because the reduction ran once for that cell.
    shared = [(a, b) for a in eachindex(points), b in eachindex(points)
              if a < b && ice[:cell_id][a] == ice[:cell_id][b]]
    if isempty(shared)
        println("\nNo two points fell in the same 500 m MODIS cell.")
    else
        println("\nPoints sharing a MODIS cell (value derived once, reused):")
        for (a, b) in shared
            @printf("  points %d and %d → %s   identical: %s\n", a, b, ice[:cell_id][a],
                    all(i -> isequal(bsa[a, i], bsa[b, i]), eachindex(years)))
        end
    end

    # One value per point, averaged over years — the form GEMB wants for bare-ice albedo.
    println("\n", "="^70)
    println("Multi-year mean (feed this to GEMB instead of a tuned regional constant)")
    println("="^70)
    for p in eachindex(points)
        b = filter(!isnan, collect(bsa[p, :]))
        w = filter(!isnan, collect(wsa[p, :]))
        if isempty(b)
            @printf("  point %d: unresolved — try a lower `albedo_range` floor, a smaller \
                    `min_samples`, or `qa_keep=[0, 1]`\n", p)
        else
            @printf("  point %d: BSA %.3f   WSA %.3f\n", p, mean(b),
                    isempty(w) ? NaN : mean(w))
        end
    end

    # Black-sky is generally the darker of the two over snow and ice at typical solar
    # geometry. A sanity check on layer identity, not a physical law.
    both = [(bsa[i], wsa[i]) for i in eachindex(bsa) if !isnan(bsa[i]) && !isnan(wsa[i])]
    if !isempty(both)
        @printf("\nBSA ≤ WSA at %d of %d resolved point-years (expected for most).\n",
                count(p -> p[1] <= p[2], both), length(both))
    end

    return nothing
end

"""
    write_result(ice, years, output_dir) -> String

Write the point stack to NetCDF and return the path.

The `:cell_id` layer is dropped: it is a `String` vector, and the numeric
`latitude`/`longitude` of the sampled cell centre already identify the cell unambiguously —
two points with equal coordinates there shared a cell.
"""
function write_result(ice, years, output_dir)
    mkpath(output_dir)
    path = joinpath(output_dir,
                    "glacier_ice_albedo_modis_$(first(years))_$(last(years)).nc")
    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        defDim(ds, "point", size(ice[:albedo_bsa], 1))
        defDim(ds, "time", size(ice[:albedo_bsa], 2))
        defVar(ds, "time", [year(t) for t in lookup(ice, Ti)], ("time",))
        for k in keys(ice)
            k === :cell_id && continue
            a = ice[k]
            dnames = ndims(a) == 2 ? ("point", "time") : ("point",)
            v = defVar(ds, String(k), parent(a), dnames)
            for (mk, mv) in metadata(a)
                v.attrib[mk] = mv
            end
        end
        # Vector-valued metadata (`albedo_range`, `qa_keep`, `years`) is stringified: NetCDF
        # attributes are scalars or homogeneous arrays, and these are provenance, not data.
        for (mk, mv) in metadata(ice)
            ds.attrib[mk] = mv isa AbstractVector ? string(mv) : mv
        end
    end
    return path
end

## ------------------------------------------------------------------------- execute

# Runs on `include`, and `run_example` stays callable afterwards with other settings.
ice = run_example()

## ---------------------------------------------------------------------- variations
#
# Heavily dust- or algae-darkened ablation zones: the default 0.3 floor is a
# *glaciological* bound, not a physical one, and will clip the very signal being measured.
# Lower it and the darkest retrievals survive:
#
#     ice = run_example(2019; albedo_range = (0.15, 1.0))
#
# High latitude or persistent cloud: admitting magnitude inversions (QA class 1) roughly
# doubles the sample count at the cost of weaker retrievals, and relaxing `min_samples`
# resolves more points:
#
#     ice = run_example(2019; qa_keep = [0, 1], min_samples = 10)
#
# The spectral bands are in the same granules, so requesting more layers costs NO extra
# download — only the reduction runs again:
#
#     ice = run_example(2019; layers = (:Albedo_BSA_shortwave, :Albedo_WSA_shortwave,
#                                       :Albedo_BSA_vis, :Albedo_BSA_nir))
#     ice[:albedo_bsa_vis]
#
# Keep the granules if you plan to re-run with different QC (costs ~70 MB per tile-date of
# disk, but makes the second pass free):
#
#     ice = run_example(2019; keep_granules = true)
#
# `stride` is the blunt volume fallback — unlike `doy_range` it thins the melt season too:
#
#     ice = run_example(2019; doy_range = nothing, stride = 4)
#
# Southern hemisphere. `:melt_season` handles this automatically, but note what a *calendar*
# year means there: the wrapping window pools the tail of the 2018/19 season with the start of
# 2019/20. For a darkest-5 % statistic that is fine — both are bare-ice states of the same
# glacier. A point list straddling the equator falls back to the whole year, so split it:
#
#     patagonia = [(-49.30, -73.20), (-50.50, -73.50)]   # Southern Patagonian Icefield
#     ice = run_example(2019; points = patagonia)                     # → (300, 110)
#     ice = run_example(2019; points = patagonia, doy_range = (330, 60))  # narrower, wrapping
#
# Inspect why points are missing — the observation count is reported even where the albedo
# is NaN, precisely so this is answerable after the fact:
#
#     n = ice[:n_valid_observations_bsa]
#     count(iszero, n), count(>(0), n)
