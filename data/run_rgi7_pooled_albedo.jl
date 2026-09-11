#!/usr/bin/env julia
#
# Global RGI 7.0 bare-ice albedo, pooled over the whole MCD43A3 record.
#
# One value per glacierized 500 m cell: the mean of the darkest 5 % of *every* valid retrieval
# from 2000–2025, with the number of retrievals behind it. This is not the per-year product
# `run_rgi7_ice_albedo_modis.jl` writes, and not a reduction of it — see
# `src/pooled_ice_albedo.jl` for why coverage variability between years makes the pooled form
# the more useful statistic.
#
# THIS SCRIPT DOWNLOADS NOTHING. It re-folds the per-date sample cache that the per-year run
# left on disk (~184 GB across 8983 dates), which holds raw DNs and the QA class rather than
# post-QC values — so the percentile, the QA whitelist and the albedo range are all free to
# change here, while the date list, cell list and layer set are fixed at download time.
#
# Runs in the package environment; every dependency is already one.
#
#   julia --project=. data/run_rgi7_pooled_albedo.jl                     # both hemispheres
#   julia --project=. data/run_rgi7_pooled_albedo.jl north               # one hemisphere
#   RGI7_ALBEDO_FLOOR=0.20 julia --project=. data/run_rgi7_pooled_albedo.jl
#
# Two hemispheres are independent folds and can be run as two concurrent processes; each is
# single-threaded and holds its own accumulator (4.3 GiB north, 1.3 GiB south).

using GEMB_ClimateForcing
using CodecZlib
using Dates
using Printf
using DimensionalData

const G = GEMB_ClimateForcing

const CACHE_PATH = get(ENV, "RGI7_CACHE_PATH", G._default_modis_cache())
const OUTPUT_DIR = get(ENV, "RGI7_ALBEDO_DIR", @__DIR__)
const PERCENTILE = parse(Float64, get(ENV, "RGI7_ALBEDO_PERCENTILE", "0.05"))
# 0.25, not the per-year run's 0.30: every per-year file has `minimum == exactly 0.3000`, i.e.
# that floor was clipping the dark tail it was meant to bound.
const FLOOR = parse(Float64, get(ENV, "RGI7_ALBEDO_FLOOR", "0.25"))
const CEILING = parse(Float64, get(ENV, "RGI7_ALBEDO_CEILING", "1.0"))

say(args...) = (println(args...); flush(stdout))

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


pooled_path(dir, hemi) = joinpath(dir, "rgi7_ice_albedo_pooled_$(hemi).csv.gz")

"""
    write_pooled(path, cells, ice)

One hemisphere's pooled result as gzipped CSV, keyed by cell and in the cells' canonical order.

Plain text for the same reason the per-year intermediates are: the GeoParquet step runs in a
different environment, and these files outlive package versions.
"""
function write_pooled(path::AbstractString, cells::Vector{NTuple{4,Int}}, ice)
    bsa = parent(ice[:albedo_bsa])
    wsa = parent(ice[:albedo_wsa])
    nb = parent(ice[:n_valid_observations_bsa])
    nw = parent(ice[:n_valid_observations_wsa])
    kb = parent(ice[:k_used_bsa])
    kw = parent(ice[:k_used_wsa])
    md = DimensionalData.metadata(ice)
    tmp = path * ".part"
    open(GzipCompressorStream, tmp, "w") do io
        println(io, "# Global RGI 7.0 glacier bare-ice albedo, POOLED over the whole record.")
        println(io, "# Written by data/run_rgi7_pooled_albedo.jl.")
        println(io, "#")
        @printf(io, "# percentile   : %g\n", md["percentile"])
        @printf(io, "# albedo_range : %g,%g\n", md["albedo_range"]...)
        @printf(io, "# qa_keep      : %s\n", join(md["qa_keep"], ","))
        @printf(io, "# n_dates      : %d  (%s … %s)\n", md["n_dates"], md["date_range"]...)
        @printf(io, "# cells_key    : %s\n", md["cells_key"])
        println(io, "#")
        println(io, "# albedo_bsa/wsa: black-sky / white-sky broadband shortwave, unitless.")
        println(io, "#                 Mean of the darkest `percentile` of ALL valid")
        println(io, "#                 retrievals in the record. NaN where none survived QC.")
        println(io, "#                 NO upper clamp — MCD43A3 albedo can exceed 1.0.")
        println(io, "# n_valid_*     : retrievals passing QA and albedo_range, whole record.")
        println(io, "# k_used_*      : how many of the darkest were averaged,")
        println(io, "#                 ceil(percentile * n_valid).")
        println(io, "#")
        println(io, "# NOT a per-year value and NOT a reduction of per-year values. Judge a")
        println(io, "# cell by n_valid; no sample-count threshold is applied.")
        println(io, "h,v,row,col,albedo_bsa,albedo_wsa,n_valid_bsa,n_valid_wsa,k_used_bsa,k_used_wsa")
        for (i, c) in enumerate(cells)
            @printf(io, "%d,%d,%d,%d,%.4f,%.4f,%d,%d,%d,%d\n",
                    c[1], c[2], c[3], c[4], bsa[i], wsa[i], nb[i], nw[i], kb[i], kw[i])
        end
    end
    mv(tmp, path; force = true)
    return path
end

function report(ice, label)
    bsa = parent(ice[:albedo_bsa])
    nb = parent(ice[:n_valid_observations_bsa])
    res = count(!isnan, bsa)
    n = length(bsa)
    say(@sprintf("  %s: %d of %d cells resolved (%.1f %%)", label, res, n, 100 * res / n))
    if res > 0
        vals = filter(!isnan, bsa)
        say(@sprintf("    albedo_bsa  min %.3f  mean %.3f  max %.3f", minimum(vals),
                     sum(vals) / length(vals), maximum(vals)))
        counts = filter(>(0), nb)
        say(@sprintf("    n_valid     min %d  median %d  max %d", minimum(counts),
                     sort(counts)[cld(length(counts), 2)], maximum(counts)))
    end
    return nothing
end

function run_pooled(hemispheres = (:north, :south); cache_path = CACHE_PATH,
                    dir = OUTPUT_DIR, percentile = PERCENTILE,
                    albedo_range = (FLOOR, CEILING), force = false)
    say("=" ^ 74)
    say("Global RGI 7.0 pooled bare-ice albedo — re-fold from cache, no download")
    say("=" ^ 74)
    say(@sprintf("cache        : %s", cache_path))
    say(@sprintf("output       : %s", dir))
    say(@sprintf("percentile   : %g", percentile))
    say(@sprintf("albedo_range : %g … %g", albedo_range...))
    say("")
    assert_bulk_cache(cache_path,
        "This script only re-folds samples that are already cached; it downloads " *
        "nothing, so an empty cache means there is nothing to pool.")

    t = rgi7_modis_cells()
    north, south = rgi7_hemisphere_split(t)
    index_of = Dict(:north => north, :south => south)

    for hemi in hemispheres
        path = pooled_path(dir, hemi)
        if isfile(path) && !force
            say("$(hemi): $(basename(path)) already present, skipping (force=true to redo)")
            continue
        end
        cells = rgi7_modis_unique_cells(t; index = index_of[hemi])
        say(@sprintf("%s: %d cells", hemi, length(cells)))
        t0 = time()
        ice = pool_ice_albedo_from_cache(cells; cache_path = cache_path,
                                        percentile = percentile,
                                        albedo_range = albedo_range)
        report(ice, String(hemi))
        write_pooled(path, cells, ice)
        say(@sprintf("  wrote %s (%.1f MB) in %.1f min", basename(path),
                     filesize(path) / 1e6, (time() - t0) / 60))
        # The accumulator is GiB-scale; drop it before the next hemisphere.
        ice = nothing
        GC.gc()
        say("")
    end
    say("done.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    hemis = isempty(ARGS) ? (:north, :south) : Tuple(Symbol.(ARGS))
    run_pooled(hemis)
end
