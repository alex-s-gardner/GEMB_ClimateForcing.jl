#!/usr/bin/env julia
# Global glacier bare-ice albedo from MODIS MCD43A3, at every RGI 7.0 glacier cell.
#
# Reads the vendored cell list (data/rgi7_modis_cells.csv.gz, built by
# data/make_rgi7_modis_cells.jl) and evaluates the darkest-percentile annual albedo at all
# 3,360,716 distinct MCD43A3 cells, one calendar year at a time. Writes a gzipped-CSV
# intermediate per (year, hemisphere); data/make_rgi7_ice_albedo_parquet.jl turns those into
# the GeoParquet product.
#
#   export EARTHDATA_TOKEN="..."
#   julia --project=. -t auto data/run_rgi7_ice_albedo_modis.jl 2019 [2020 ...]
#
# Needs no scratch environment — everything here is a package dependency. Only the Parquet
# writer needs GeoParquet.jl, which is why it is a separate script.
#
# !! THIS DOWNLOADS ~1 TB PER YEAR. !!
#
# The cost is granule volume and is not reducible: MCD43A3 has no server-side subsetting, so
# the unit of transfer is a whole granule of which ~85 % is layers never read. Measured from
# real CMR granule sizes (three probe dates per hemisphere, each INSIDE its own melt season):
#
#   northern  63 tiles x 171 dates (doy 120-290)  4.36 GB/date  =  745 GB
#   southern  39 tiles x 176 dates (doy 300-110)  1.48 GB/date  =  261 GB
#                                                        TOTAL  =  0.98 TB/yr
#
# Note the package's own `@info "Download volume estimate"` assumes 70 MB/granule and so
# reads ~1.2 TB; the figures above are measured and are the ones to plan against. Southern
# granules are far smaller (median 6 MB) because those tiles are mostly ocean, and one
# southern tile has no granule at all — a permanent gap, handled rather than fatal.
#
# Sampling the south in July would have suggested 0.79 GB/date; in its OWN season it is 1.48.
# Measure each hemisphere inside its own window or the estimate is off by 2x.
#
# Peak disk stays near one date's tiles (4.4 GB north, 1.5 GB south) because
# keep_granules=false deletes each date's granules once folded. Re-running is cheap: the
# per-date sample cache means a killed run re-downloads nothing.
#
# WHY TWO CALLS PER YEAR, NOT ONE
#
# `doy_range=:melt_season` returns `nothing` — the WHOLE YEAR — for any cell list spanning
# both hemispheres, and a global list always does. That is 365 dates instead of ~171, i.e.
# ~2.4 TB instead of ~1.16 TB, for no benefit. So each hemisphere is run separately with an
# EXPLICIT window. The split is exact and free: _MODIS_UL_Y == 9 * _MODIS_TILE_SPAN_M to the
# bit, so `v <= 8` is precisely the northern hemisphere and no MODIS tile straddles the
# equator. The two tile sets are therefore disjoint, so nothing is downloaded twice.
#
# Skipping the southern half is not an option: it is 23.2 % of global glacier area, region 19
# (Subantarctic and Antarctic Islands) alone being 133,432 km².
#
# NOTE on the southern window: (300, 110) wraps New Year, so a calendar year pools the tail of
# one melt season with the start of the next. That is deliberate and documented in
# src/glacier_ice_albedo_modis.jl — both are bare-ice states of the same glacier — and the
# window actually used is recorded per row in the output.

using GEMB_ClimateForcing
using CodecZlib
using Dates
using DimensionalData
using Printf

const G = GEMB_ClimateForcing

## ------------------------------------------------------------------------------- configuration

const YEARS = [2019]

# Granule cache and sample cache. Wants to live on a filesystem with room for one date's
# tiles plus ~2.8 GB of sample cache per year — NOT under tempdir(), which the OS may clear
# out from under a multi-day run.
const CACHE_PATH = get(ENV, "RGI7_CACHE_PATH", G._default_modis_cache())
# Overridable to keep a test run's intermediates out of data/, matching the Parquet writer.
const OUTPUT_DIR = get(ENV, "RGI7_ALBEDO_DIR", joinpath(@__DIR__))

# The statistic. Deliberately the package defaults; see compute_glacier_ice_albedo_modis.
const PERCENTILE = 0.05
const MIN_SAMPLES = 30
const ALBEDO_RANGE = (0.3, 1.0)
# Every full BRDF inversion — see MCD43A3_QA_KEEP. Add the odd classes to admit
# magnitude inversions if min_samples turns out unreachable at high latitude.
const QA_KEEP = MCD43A3_QA_KEEP

# Bandwidth-bound with no queue limiter, so higher than the package default of 4. Measure
# before pushing past ~16: LP DAAC throttles, and a stalled connection costs a retry cycle.
const MAX_CONCURRENT_DOWNLOADS = 8

# Measured mean granule volume per date, per hemisphere, in GB — from real CMR sizes over
# three probe dates inside each melt season (north 4.20/4.27/4.60, south 1.46/1.48/1.51).
# Used for the plan printout instead of the package's 70 MB/granule assumption, which
# overestimates by ~20 % globally and by ~2.5x in the south.
const MEASURED_GB_PER_DATE = Dict(:north => 4.36, :south => 1.48)

const HEMISPHERES = (:north, :south)

## ------------------------------------------------------------------------------------ helpers

"""
    say(args...)

`println` that flushes. Julia block-buffers a redirected `stdout`, so a `nohup`-ed run would
otherwise produce an empty log for hours and look hung — measured at 0 bytes after 3.5 minutes
while 966 MB were happily downloading. Every progress line in this script goes through here or
through [`@printf_flush`](@ref).
"""

"""
    assert_bulk_cache(path, why)

Refuse to run with the cache under `tempdir()`.

These runs are hundreds of gigabytes and the cache is the whole point of being able to resume
or re-fold, so a temp-directory cache is never what was intended: it is reaped by the OS, and
`/tmp` is usually on the root filesystem. Failing here costs a second; discovering it later
costs a re-download.
"""
function assert_bulk_cache(path::AbstractString, why::AbstractString)
    startswith(abspath(path), abspath(tempdir())) || return nothing
    error("""
        The MCD43A3 cache resolves to $(path), which is under tempdir().
        $(why)
        Set a bulk-storage location on a large volume (NOT a home directory, which is
        commonly a small SSD with a quota):
            export GEMB_CACHE_PATH=/big/volume/gemb_cache   # all products
            export RGI7_CACHE_PATH=/big/volume/MCD43A3.061  # this cache only
        """)
end

function say(args...)
    println(args...)
    flush(stdout)
    return nothing
end

"""
    @printf_flush(fmt, args...)

`@printf` to `stdout` that flushes afterwards — same reason as [`say`](@ref). Deliberately not
applied to the `@printf(io, ...)` calls that write the intermediates: those go to a file that is
closed on completion, so buffering there is correct.
"""
macro printf_flush(fmt, args...)
    # `fmt` is interpolated rather than escaped: `@printf` requires a *literal* format string
    # at macro-expansion time and rejects an escaped expression with "First argument to
    # `@printf` after `io` must be a format string".
    return quote
        @printf($fmt, $(map(esc, args)...))
        flush(stdout)
    end
end

"""
    hemisphere_config(cells, hemi) -> (; index, cells, doy_range, label)

Cell indices, the distinct cells, and the explicit melt-season window for one hemisphere.

`doy_range` is passed explicitly rather than as `:melt_season` so the window appears in the
log and in the output, and so a future change to the resolver cannot silently double the
download.
"""
function hemisphere_config(t::G.RGI7ModisCells, hemi::Symbol)
    north, south = rgi7_hemisphere_split(t)
    if hemi === :north
        return (; index = north, cells = rgi7_modis_unique_cells(t; index = north),
                  doy_range = G._MODIS_MELT_SEASON_NORTH, label = "north")
    elseif hemi === :south
        return (; index = south, cells = rgi7_modis_unique_cells(t; index = south),
                  doy_range = G._MODIS_MELT_SEASON_SOUTH, label = "south")
    end
    throw(ArgumentError("hemisphere must be :north or :south, got $(hemi)"))
end

intermediate_path(dir::AbstractString, year::Integer, hemi::Symbol) =
    joinpath(dir, @sprintf("rgi7_ice_albedo_modis_%d_%s.csv.gz", year, hemi))

"""
    already_done(path, ncell) -> Bool

Whether an intermediate is present and has the expected number of data rows.

The coarse half of resumability: a killed run is restarted by re-invoking the script, and a
completed (year, hemisphere) is skipped without touching CMR or the network. The fine half is
the per-date sample cache inside `compute_glacier_ice_albedo_modis`, which makes a crash
*within* a hemisphere re-download nothing.
"""
function already_done(path::AbstractString, ncell::Integer)
    isfile(path) || return false
    n = 0
    try
        open(GzipDecompressorStream, path) do io
            for line in eachline(io)
                (isempty(line) || startswith(line, '#') || startswith(line, "h,")) && continue
                n += 1
            end
        end
    catch e
        @warn "unreadable intermediate, will recompute" path exception = e
        return false
    end
    n == ncell && return true
    @warn "intermediate has the wrong row count, will recompute" path rows = n expected = ncell
    return false
end

"""
    write_intermediate(path, cells, ice, doy_range)

One `(year, hemisphere)` result as gzipped CSV, keyed by cell.

Plain text on purpose: these files are the handoff to a script running in a *different*
environment (the Parquet writer needs GeoParquet.jl, which is deliberately not a package
dependency), and they outlive package versions. Keyed by `(h, v, row, col)` and written in the
cells' own canonical order, so the Parquet step is a merge, not a hash join.
"""
function write_intermediate(path::AbstractString, cells::Vector{NTuple{4,Int}}, ice,
                            doy_range::Tuple{Int,Int})
    bsa = parent(ice[:albedo_bsa])
    wsa = parent(ice[:albedo_wsa])
    nb = parent(ice[:n_valid_observations_bsa])
    nw = parent(ice[:n_valid_observations_wsa])
    tmp = path * ".part"
    open(GzipCompressorStream, tmp, "w") do io
        println(io, "# Glacier bare-ice albedo at RGI 7.0 MCD43A3 cells, one hemisphere-year.")
        println(io, "# Written by data/run_rgi7_ice_albedo_modis.jl; intermediate for")
        println(io, "# data/make_rgi7_ice_albedo_parquet.jl. Not the published product.")
        println(io, "#")
        @printf(io, "# doy_range     : %d,%d%s\n", doy_range[1], doy_range[2],
                doy_range[1] > doy_range[2] ? "  (wraps New Year)" : "")
        @printf(io, "# percentile    : %g\n", PERCENTILE)
        @printf(io, "# min_samples   : %d\n", MIN_SAMPLES)
        @printf(io, "# albedo_range  : %g,%g\n", ALBEDO_RANGE[1], ALBEDO_RANGE[2])
        @printf(io, "# qa_keep       : %s\n", join(QA_KEEP, ","))
        println(io, "#")
        println(io, "# albedo_bsa/wsa: black-sky / white-sky broadband shortwave, unitless.")
        println(io, "#                 NaN where unresolved. NO upper clamp — MCD43A3 albedo")
        println(io, "#                 can legitimately exceed 1.0.")
        println(io, "# n_valid_*     : retrievals contributing to the percentile mean.")
        println(io, "h,v,row,col,albedo_bsa,albedo_wsa,n_valid_bsa,n_valid_wsa")
        for (i, c) in enumerate(cells)
            @printf(io, "%d,%d,%d,%d,%s,%s,%d,%d\n", c[1], c[2], c[3], c[4],
                    isnan(bsa[i, 1]) ? "NaN" : @sprintf("%.4f", bsa[i, 1]),
                    isnan(wsa[i, 1]) ? "NaN" : @sprintf("%.4f", wsa[i, 1]),
                    nb[i, 1], nw[i, 1])
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    summarize(ice, label)

Print the resolved fraction and albedo distribution for one hemisphere-year.

`NaN` is filtered explicitly: it means the cell failed QC or fell below `min_samples`, which
is expected at high latitude and for `forced` cells, and averaging over it would silently
propagate.
"""
function summarize(ice, label::AbstractString)
    for key in (:albedo_bsa, :albedo_wsa)
        a = vec(parent(ice[key]))
        ok = filter(!isnan, a)
        if isempty(ok)
            @printf_flush("  %-11s %-6s no cell resolved\n", label, key)
            continue
        end
        s = sort(ok)
        @printf_flush("  %-11s %-6s resolved %7d/%7d (%5.1f%%)  min %.3f  p05 %.3f  med %.3f  p95 %.3f  max %.3f\n",
                label, key, length(ok), length(a), 100 * length(ok) / length(a),
                s[1], s[max(1, round(Int, 0.05 * end))], s[cld(end, 2)],
                s[min(end, round(Int, 0.95 * end))], s[end])
    end
    return nothing
end

## ---------------------------------------------------------------------------------- the run

"""
    run_rgi7_albedo(years = YEARS; hemispheres = HEMISPHERES, kwargs...)

Compute and write the global product, one `(year, hemisphere)` at a time.

Idempotent: a completed pair is skipped, so re-invoking after a crash resumes. Keyword
arguments are forwarded to [`compute_glacier_ice_albedo_modis`](@ref).
"""
function run_rgi7_albedo(years = YEARS; hemispheres = HEMISPHERES, cache_path = CACHE_PATH,
                         output_dir = OUTPUT_DIR, dry_run::Bool = false, kwargs...)
    # Checked up front, so a missing token fails now rather than three hours in.
    token = try
        get_earthdata_token()
    catch e
        error("""
        No NASA Earthdata token available ($(sprint(showerror, e))).

        Mint one at https://urs.earthdata.nasa.gov/profile/edit/user_tokens (60-day life,
        MAXIMUM TWO concurrent — a third request 403s) and then:

            export EARTHDATA_TOKEN="your-token-here"

        Or write it to ~/.edl_token. `earthdata_token_from_netrc()` will mint one from
        ~/.netrc, but it is a deliberate, explicit call for exactly that reason.
        """)
    end

    year_list = vec(collect(years))
    t = rgi7_modis_cells()

    say("="^78)
    say("Global RGI 7.0 glacier bare-ice albedo from MODIS MCD43A3 v061")
    say("="^78)
    say(t)
    @printf_flush("years          : %s\n", join(year_list, ", "))
    @printf_flush("cache          : %s\n", cache_path)
    @printf_flush("output         : %s\n", output_dir)
    @printf_flush("statistic      : darkest %g%% mean, min_samples=%d, albedo_range=%s, qa_keep=%s\n",
            100 * PERCENTILE, MIN_SAMPLES, ALBEDO_RANGE, QA_KEEP)
    @printf_flush("token          : resolved (%d chars)\n", length(token))
    assert_bulk_cache(cache_path,
        "A full year is ~1 TB of granule download, and the per-date sample cache is what " *
        "makes the run resumable and re-foldable.")

    total_gb = 0.0
    plan = Dict{Symbol,Any}()
    for hemi in hemispheres
        cfg = hemisphere_config(t, hemi)
        ntile = length(unique(c -> (c[1], c[2]), cfg.cells))
        # Summed over every requested year, not the first year scaled up: 2000 is a partial
        # year (the record starts 2000-02-16), so a southern window loses its first 45 days and
        # extrapolating from it under-reports every later year by ~25 %.
        ndates = [length(G._modis_dates(y, cfg.doy_range, 1)) for y in year_list]
        gb = sum(ndates) * MEASURED_GB_PER_DATE[hemi]
        total_gb += gb
        plan[hemi] = cfg
        @printf_flush("  %-6s %9d cells  %3d tiles  doy %3d-%-3d  %3d dates%s  ~%6.0f GB total\n",
                cfg.label, length(cfg.cells), ntile, cfg.doy_range[1], cfg.doy_range[2],
                sum(ndates),
                allequal(ndates) ? "" : " ($(minimum(ndates))–$(maximum(ndates))/yr)", gb)
        # The sample cache is keyed on a hash of the WHOLE cell list, so regenerating the
        # vendored table orphans (does not delete) every cached date. Log it to make that
        # visible when a "resume" unexpectedly re-downloads.
        @printf_flush("         sample-cache key %s\n", G._modis_cells_key(cfg.cells))
    end
    @printf_flush("\nESTIMATED TOTAL DOWNLOAD: %.1f TB over %d year(s)\n", total_gb / 1024,
            length(year_list))
    say("="^78)

    if dry_run
        say("\ndry_run=true — nothing downloaded.")
        return nothing
    end

    mkpath(output_dir)
    for y in year_list, hemi in hemispheres
        cfg = plan[hemi]
        path = intermediate_path(output_dir, y, hemi)
        if already_done(path, length(cfg.cells))
            @printf_flush("\n[%d %s] already complete, skipping: %s\n", y, cfg.label, basename(path))
            continue
        end
        @printf_flush("\n[%d %s] %d cells, doy %d-%d — starting %s\n", y, cfg.label,
                length(cfg.cells), cfg.doy_range[1], cfg.doy_range[2], Dates.now())
        t0 = time()
        ice = compute_glacier_ice_albedo_modis(cfg.cells, y;
                                              percentile = PERCENTILE,
                                              min_samples = MIN_SAMPLES,
                                              albedo_range = ALBEDO_RANGE,
                                              qa_keep = QA_KEEP,
                                              doy_range = cfg.doy_range,
                                              token = token,
                                              cache_path = cache_path,
                                              keep_granules = false,
                                              cell_id = false,
                                              max_concurrent_downloads = MAX_CONCURRENT_DOWNLOADS,
                                              progress = true, verbose = true,
                                              kwargs...)
        @printf_flush("[%d %s] folded in %.1f h\n", y, cfg.label, (time() - t0) / 3600)
        summarize(ice, "$(y) $(cfg.label)")

        # A year is only written if every date got every granule that exists. `degraded_dates`
        # lists dates that lost granules to a FAILED DOWNLOAD rather than to a genuine record
        # gap; those dates are deliberately left out of the sample cache, so simply re-invoking
        # this script retries exactly them and nothing else.
        #
        # Refusing to write is the conservative choice and it is deliberate: `already_done`
        # would otherwise skip this year forever, silently baking a network outage into a
        # multi-decade product. Measured during the first 2000-2025 attempt: an LP DAAC DNS
        # outage cost all 63 northern tiles of 2002-09-09.
        nbad = get(DimensionalData.metadata(ice), "n_degraded_dates", 0)
        if nbad > 0
            bad = get(DimensionalData.metadata(ice), "degraded_dates", String[])
            say("[$(y) $(cfg.label)] NOT WRITING: $(nbad) date(s) lost granules to failed " *
                "downloads, not to a record gap.")
            say("[$(y) $(cfg.label)]   affected: " * join(first(bad, 8), ", ") *
                (nbad > 8 ? " … (+$(nbad - 8) more)" : ""))
            say("[$(y) $(cfg.label)]   those dates are uncached, so re-run this script to " *
                "retry just them once the service is back.")
            ice = nothing
            GC.gc()
            continue
        end

        write_intermediate(path, cfg.cells, ice, cfg.doy_range)
        @printf_flush("[%d %s] wrote %s (%.1f MB)\n", y, cfg.label, basename(path),
                filesize(path) / 1024^2)
        ice = nothing
        GC.gc()
    end
    say("\nDone. Reduce several years to a climatology with, e.g.:")
    say("  rgi7_ice_albedo_climatology(", first(year_list), ":", last(year_list),
            "; reduction = median)")
    say("\nOr build the GeoParquet product with:")
    say("  julia --project=/some/scratch/env data/make_rgi7_ice_albedo_parquet.jl ",
        join(year_list, " "))
    return nothing
end

## -------------------------------------------------------------------------------- execute

# Years from the command line, else YEARS. Deliberately not run on plain `include` — unlike
# examples/, this is a ~1.2 TB/year job and must be asked for explicitly.
if abspath(PROGRAM_FILE) == @__FILE__
    years = isempty(ARGS) ? YEARS : parse.(Int, ARGS)
    run_rgi7_albedo(years)
end

## ---- variations
#
# Plan only, no download — prints cells, tiles, dates and the volume estimate:
#   run_rgi7_albedo([2019]; dry_run = true)
#
# One hemisphere at a time, e.g. to run the two halves on different machines:
#   run_rgi7_albedo([2019]; hemispheres = (:north,))
#
# Keep granules, making a re-run with different QC free at ~1.2 TB of disk:
#   run_rgi7_albedo([2019]; keep_granules = true)
#
# A cheap smoke test before committing to a real year — three dates, one tile's worth of
# cells, min_samples relaxed so they resolve:
#   t = rgi7_modis_cells()
#   cells = filter(c -> (c[1], c[2]) == (17, 2), rgi7_modis_unique_cells(t))[1:200]
#   compute_glacier_ice_albedo_modis(cells, 2019; doy_range = (180, 182), min_samples = 1,
#                                    cell_id = false, cache_path = CACHE_PATH)
