# Annual "darkest-5%" broadband surface albedo from the C3S 300 m product.
#
# For every pixel and every calendar year in the requested range this computes the
# **mean of the lowest 5 % of valid broadband albedo observations recorded that year**
# (the 10-daily product gives ~36 observations/year, so 5 % is the darkest 1–2 of them).
# That statistic is a proxy for the bare-ice / maximally-darkened surface state, which is
# why it is worth isolating from the annual mean.
#
# Quality control, applied per observation before any statistic is formed:
#   1. missing / `_FillValue` pixels (NCDatasets turns these into `missing` or `NaN`)
#   2. albedo outside the physical range 0–1
#   3. QFLAG classes that undermine the retrieval — by default the `no_obs_in_last_decade_*`
#      bits, i.e. no actual observation backed that band in the 10-day window. The bit
#      table is read from each file's own `flag_masks`/`flag_meanings` attributes rather
#      than hardcoded, because the QFLAG convention is version-specific. Note the v3.1
#      legend is *not* a land/cloud mask (see `QFLAG_REJECT_PATTERNS`); sea and cloud
#      pixels come through as fill values and are caught by step 1.
#   4. optionally, observations whose companion `_ERR` uncertainty exceeds `MAX_ERROR`
#   5. pixel-years with fewer than `MIN_SAMPLES` surviving observations are left missing,
#      so a "5 % minimum" is never computed from two cloudy scenes
#
# Data is *ordered* from CDS (async jobs) and cached on disk, so the first run of a
# multi-year request is slow and reruns are nearly free. One variable-year costs ~3 CDS
# jobs; see `satellite_albedo`'s docstring on the cost limit.
#
# **Server-side latency is long and highly variable, and not simply a function of job size**:
# observed on this dataset — a 3-timestep order finishing in ~30 min, another 3-timestep
# order still `running` after 2 h, and a 12-timestep order still `running` after 3 h.
#
# Since that latency is queue time rather than compute, the fix is concurrency, not smaller
# jobs: `satellite_albedo` submits a year's jobs together (`MAX_CONCURRENT_JOBS`) and polls
# them in parallel, which roughly halves a cold year's wall-clock. Splitting the request into
# months by hand would *undo* this. Still budget ~30–90 min per cold year, and note
# `JOB_TIMEOUT` is generous because a timeout kills the process before the finished result is
# cached, forcing a reorder.
#
# Requirements:
#   export CDS_API_KEY="your-token-here"
#   plus one-time licence acceptance at
#   https://cds.climate.copernicus.eu/datasets/satellite-albedo?tab=download
#
# Usage:
#   julia --project=. examples/albedo_annual_low_percentile.jl            # YEARS below
#   julia --project=. examples/albedo_annual_low_percentile.jl 2019 2020  # override years

using GEMB_ClimateForcing
import GEMB_ClimateForcing as GCF
using Rasters
using NCDatasets
using DimensionalData
using Dates
using Statistics
using Printf

## ---------------------------------------------------------------- configuration

# Southwest Greenland ablation zone (Kangerlussuaq transect) — small enough to run in a
# few minutes of ordering per year. Longitude is −180…180 for this product.
const EXTENT = Extent(X=(-50.0, -47.0), Y=(66.5, 67.5))

# Sentinel-3 300 m era only (2018–2024).
const YEARS = 2019:2019

const VARIABLE = :albb_dh          # broadband directional-hemispherical (black-sky)
const LAYER = :AL_DH_BB            # full-spectrum broadband albedo layer in that file
const ERR_LAYER = :AL_DH_BB_ERR    # companion uncertainty layer (set to `nothing` to skip)

const PERCENTILE = 0.05            # fraction of darkest observations to average
const MIN_SAMPLES = 10             # min valid observations per pixel-year
const MAX_ERROR = 0.20             # reject obs whose _ERR exceeds this (absolute albedo)
const ALBEDO_RANGE = (0.0, 1.0)

# QFLAG classes to reject, matched case-insensitively as substrings of the file's own
# `flag_meanings`.
#
# The v3.1 legend (19 bitmask classes, read from a delivered file) is **not** a land/cloud
# mask — there is no sea, water, cloud or shadow class in it at all. It is:
#
#   bit 0            snow_presence
#   bits 1–9         no_obs_in_last_decade_for_<band>   (Oa03 Oa04 Oa07 Oa17 Oa21 S1 S2 S5 S6)
#   bits 10–18       brdf_warning_for_<band>            (same nine bands)
#
# So: sea, cloud and shadow pixels arrive as `_FillValue`/`missing` instead, and are
# already dropped by the missing-value check — QFLAG's job here is narrower.
#
# `no_obs_in_last_decade_*` marks a band with no observation backing the retrieval in that
# 10-day window, which is exactly the cloud-persistent case that produces spurious low
# albedo, so it is rejected by default. `snow_presence` is deliberately *kept* — it is the
# signal, not an error.
const QFLAG_REJECT_PATTERNS = [
    "no_obs_in_last_decade",
    # Generic names from other product versions/eras, harmless when absent from the legend.
    "sea", "water", "ocean", "cloud", "shadow",
    "invalid", "fail", "error", "no_data", "nodata", "missing", "unusable", "unfilled",
]

# Also reject `brdf_warning_*` (the BRDF inversion for a band was flagged as poorly
# constrained). Off by default: it is a warning rather than a failure, and on a snow/ice
# surface it fires often enough to thin the sample materially.
const QFLAG_REJECT_BRDF_WARNING = false

# Per-job give-up time, in seconds. Well above `satellite_albedo`'s 1800 s default: a
# 12-timestep albedo order was measured still `running` after 1840 s, so the default
# aborts a job that would have succeeded.
const JOB_TIMEOUT = 10800

# CDS jobs to keep in flight. Jobs are queued and processed server-side, so wall-clock is
# set by how many wait in parallel, not by CPU. Measured with six one-month orders: three
# finished in 27–32 min, all six in 83 min, versus ~3 h in series. The gain saturates —
# CDS runs only ~3 of an account's jobs at once — and it rejects submissions past a
# per-dataset queue limit, so raising this much above 6 buys nothing.
const MAX_CONCURRENT_JOBS = 6

const OUTPUT_DIR = joinpath(@__DIR__, "output")

## ------------------------------------------------------------------- QFLAG table

"""
    qflag_table(path) -> Vector{@NamedTuple{meaning::String, mask::Int, is_value::Bool}}

The QFLAG legend as declared by the file itself: CF `flag_meanings` paired with either
`flag_masks` (bitwise) or `flag_values` (exact value). Returns an empty vector when the
attributes are absent, so callers can degrade to range-only QC instead of guessing bits.
"""
function qflag_table(path::AbstractString)
    NCDataset(String(path), "r") do ds
        haskey(ds, "QFLAG") || return NamedTuple[]
        att = ds["QFLAG"].attrib
        meanings = get(att, "flag_meanings", nothing)
        meanings === nothing && return NamedTuple[]
        names = split(strip(String(meanings)))
        masks = get(att, "flag_masks", nothing)
        is_value = masks === nothing
        codes = is_value ? get(att, "flag_values", nothing) : masks
        codes === nothing && return NamedTuple[]
        codes = Int.(collect(codes))
        n = min(length(names), length(codes))
        return [(meaning=names[i], mask=codes[i], is_value=is_value) for i in 1:n]
    end
end

"""
    qflag_rejects(table; patterns=QFLAG_REJECT_PATTERNS, brdf_warning=…) -> (bitmask, values, matched)

Split the legend into a single OR-ed `bitmask` (for `flag_masks` files) and a list of
exact `values` (for `flag_values` files), covering every class whose meaning matches one
of `patterns`. `matched` names the classes for logging — QC that silently drops data is
worse than none, and this product's legend does not contain the class names one would
expect (see [`QFLAG_REJECT_PATTERNS`](@ref)).

Set `brdf_warning=true` to additionally reject the nine `brdf_warning_*` classes.
"""
function qflag_rejects(table; patterns=QFLAG_REJECT_PATTERNS,
                       brdf_warning::Bool=QFLAG_REJECT_BRDF_WARNING)
    pats = brdf_warning ? vcat(patterns, "brdf_warning") : patterns
    bitmask = 0
    values = Int[]
    matched = String[]
    for entry in table
        m = lowercase(entry.meaning)
        any(p -> occursin(p, m), pats) || continue
        push!(matched, entry.meaning)
        entry.is_value ? push!(values, entry.mask) : (bitmask |= entry.mask)
    end
    return bitmask, values, matched
end

"""
    cached_albedo_file(extent) -> Union{Nothing,String}

Path of any cached product NetCDF for `extent`, used only to read the QFLAG legend.
Reaches into the cache layout of `satellite_albedo`; returns `nothing` if that layout
changes, which downgrades QC rather than erroring.
"""
function cached_albedo_file(extent)
    try
        area = GCF._albedo_area(extent)
        store = joinpath(GCF._default_albedo_cache(), "files",
                         GCF._albedo_era_tag(), GCF._albedo_area_tag(area))
        files = GCF._albedo_cached_files(store)
        isempty(files) && return nothing
        return first(sort(collect(values(files))))
    catch err
        @debug "Could not locate a cached albedo file for the QFLAG legend" exception = err
        return nothing
    end
end

## ------------------------------------------------------------------ per-obs mask

_isbad(x::Missing) = true
_isbad(x::Number) = !isfinite(x)

"""
    valid_albedo(alb, qflag, err; bitmask, values, max_error) -> Matrix{Float32}

One timestep's albedo as a plain `Float32` matrix with every rejected observation set to
`NaN`. `qflag` and `err` may be `nothing` when those layers are unavailable.
"""
function valid_albedo(alb::AbstractMatrix, qflag, err;
                      bitmask::Integer=0, values::AbstractVector{<:Integer}=Int[],
                      max_error::Union{Nothing,Real}=nothing)
    lo, hi = ALBEDO_RANGE
    out = Matrix{Float32}(undef, size(alb))
    @inbounds for i in eachindex(out, alb)
        a = alb[i]
        if _isbad(a) || a < lo || a > hi
            out[i] = NaN32
            continue
        end
        if qflag !== nothing
            q = qflag[i]
            if _isbad(q)
                out[i] = NaN32
                continue
            end
            qi = Int(q)
            if (bitmask != 0 && (qi & bitmask) != 0) || (!isempty(values) && qi in values)
                out[i] = NaN32
                continue
            end
        end
        if err !== nothing && max_error !== nothing
            e = err[i]
            # A missing uncertainty is not itself disqualifying; an oversized one is.
            if !_isbad(e) && e > max_error
                out[i] = NaN32
                continue
            end
        end
        out[i] = Float32(a)
    end
    return out
end

## --------------------------------------------------------------- the statistic

"""
    low_percentile_mean(obs, percentile, min_samples) -> (mean, count)

Per-pixel mean of the darkest `percentile` fraction of valid observations in `obs`
(a vector of same-sized matrices, invalid entries `NaN`).

The count of averaged observations is `max(1, ceil(percentile * n_valid))`, so a pixel
with 36 valid scenes averages its darkest 2 at 5 %. Pixels with fewer than `min_samples`
valid observations yield `NaN` — the statistic is meaningless on a handful of scenes.
Returns the statistic and the number of valid observations it was drawn from.
"""
function low_percentile_mean(obs::Vector{<:AbstractMatrix{Float32}},
                             percentile::Real, min_samples::Integer)
    isempty(obs) && throw(ArgumentError("no observations supplied"))
    sz = size(first(obs))
    all(o -> size(o) == sz, obs) || throw(ArgumentError(
        "observation grids differ in size: $(unique(size.(obs)))"))

    result = fill(NaN32, sz)
    counts = zeros(Int32, sz)
    buf = Vector{Float32}(undef, length(obs))

    @inbounds for i in eachindex(result)
        n = 0
        for o in obs
            v = o[i]
            isnan(v) && continue
            n += 1
            buf[n] = v
        end
        counts[i] = n
        n >= min_samples || continue
        k = max(1, ceil(Int, percentile * n))
        vals = view(buf, 1:n)
        sort!(vals)
        result[i] = mean(view(vals, 1:k))
    end
    return result, counts
end

## ------------------------------------------------------------------------ driver

"""
    annual_low_albedo(years; extent, ...) -> (albedo::Raster, counts::Raster)

Order (or reuse) the 10-daily albedo for each year in `years` and reduce it to the
darkest-`percentile` annual mean. Both outputs carry a `Ti` dimension of `Date(year)`.
"""
function annual_low_albedo(years=YEARS; extent=EXTENT, variable=VARIABLE, layer=LAYER,
                           err_layer=ERR_LAYER, percentile=PERCENTILE,
                           min_samples=MIN_SAMPLES, max_error=MAX_ERROR,
                           timeout=JOB_TIMEOUT, verbose=true)
    year_list = collect(years)
    isempty(year_list) && throw(ArgumentError("no years requested"))

    annual = Matrix{Float32}[]
    annual_counts = Matrix{Int32}[]
    template = nothing

    for y in year_list
        verbose && @info "Loading $(y) albedo (ordering from CDS if not cached)"
        time_range = (DateTime(y, 1, 1), DateTime(y, 12, 31, 23, 59, 59))
        # One call for the whole year: `satellite_albedo` splits it into the ~3 jobs the
        # cost limit allows and runs them **concurrently**, so the year costs about one
        # job's queue time rather than three in series. Do not "help" by looping over
        # months — that serialises the ordering and is much slower.
        alb = satellite_albedo(; time_range, extent, variable, layer, timeout, verbose,
                               max_concurrent_jobs=MAX_CONCURRENT_JOBS)
        # These two reuse the *same* cached files — no extra CDS jobs are submitted.
        qf = try
            satellite_albedo(; time_range, extent, variable, layer=:QFLAG, timeout,
                             verbose=false)
        catch err
            @warn "QFLAG layer unavailable; falling back to range-only QC" exception = err
            nothing
        end
        er = if err_layer === nothing
            nothing
        else
            try
                satellite_albedo(; time_range, extent, variable, layer=err_layer,
                                 timeout, verbose=false)
            catch err
                @warn "Uncertainty layer $(err_layer) unavailable; skipping error QC" exception = err
                nothing
            end
        end

        bitmask, values, matched = if qf === nothing
            0, Int[], String[]
        else
            path = cached_albedo_file(extent)
            table = path === nothing ? NamedTuple[] : qflag_table(path)
            if isempty(table)
                @warn """
                No CF QFLAG legend (`flag_masks`/`flag_meanings`) found in the product \
                files, so QFLAG-based rejection is skipped and only missing/out-of-range \
                observations are dropped.
                """
                0, Int[], String[]
            else
                qflag_rejects(table)
            end
        end
        verbose && !isempty(matched) &&
            @info "Rejecting QFLAG classes" year = y classes = matched bitmask = bitmask values = values

        dates = collect(lookup(alb, Ti))
        obs = Matrix{Float32}[]
        for (i, d) in enumerate(dates)
            a = Array(read(alb[i]))
            q = qf === nothing ? nothing : Array(read(qf[i]))
            e = er === nothing ? nothing : Array(read(er[i]))
            m = valid_albedo(a, q, e; bitmask, values, max_error)
            push!(obs, m)
            verbose && @info @sprintf("  %s  valid %5.1f%%  median %.3f", d,
                                      100 * count(!isnan, m) / length(m),
                                      let v = filter(!isnan, vec(m))
                                          isempty(v) ? NaN : median(v)
                                      end)
        end
        isempty(obs) && error("No albedo timesteps returned for $(y)")

        template === nothing && (template = read(alb[1]))
        stat, counts = low_percentile_mean(obs, percentile, min_samples)
        push!(annual, stat)
        push!(annual_counts, counts)
        verbose && @info @sprintf("%d: %d/%d pixels resolved, mean of darkest %.0f%% = %.3f",
                                  y, count(!isnan, stat), length(stat), 100 * percentile,
                                  let v = filter(!isnan, vec(stat))
                                      isempty(v) ? NaN : mean(v)
                                  end)
    end

    xy = dims(template, (X, Y))
    ti = Ti(Date.(year_list, 1, 1))
    meta = Dict{String,Any}(
        "statistic" => "mean of lowest $(round(100 * percentile; digits=2))% of valid " *
                       "observations per calendar year",
        "variable" => variable,
        "layer" => layer,
        "percentile" => percentile,
        "min_samples" => min_samples,
        "max_error" => max_error,
        "albedo_range" => ALBEDO_RANGE,
        "qflag_reject_patterns" => QFLAG_REJECT_PATTERNS,
        "units" => "1 (dimensionless reflectance fraction, 0–1)",
        "source" => "C3S satellite surface albedo, Sentinel-3 OLCI+SLSTR 300 m, v3.1",
    )

    albedo = Raster(cat(annual...; dims=3), (xy..., ti);
                    name=Symbol("albedo_low_$(Int(round(100 * percentile)))pct"),
                    missingval=NaN32, metadata=meta)
    counts = Raster(cat(annual_counts...; dims=3), (xy..., ti);
                    name=:n_valid_observations, missingval=Int32(0), metadata=meta)
    return albedo, counts
end

## --------------------------------------------------------------------- run

function main(args=ARGS)
    years = if isempty(args)
        YEARS
    else
        ys = parse.(Int, args)
        length(ys) == 1 ? (ys[1]:ys[1]) : minimum(ys):maximum(ys)
    end

    albedo, counts = annual_low_albedo(years)

    mkpath(OUTPUT_DIR)
    tag = "$(first(years))_$(last(years))_p$(Int(round(100 * PERCENTILE)))"
    # NetCDF, not GeoTIFF: GDAL cannot write a third dimension with a `Date` lookup
    # ("no valid permutation of dimensions"), whereas NCDatasets stores `Ti` natively.
    out_path = joinpath(OUTPUT_DIR, "albedo_annual_low_$(tag).nc")
    Rasters.write(out_path, RasterStack((; albedo_low=albedo, n_valid=counts)); force=true)

    println("\n", "="^64)
    println("Mean of the lowest $(round(100 * PERCENTILE; digits=1))% annual broadband albedo")
    println("="^64)
    for (i, t) in enumerate(lookup(albedo, Ti))
        slice = filter(!isnan, vec(view(albedo, Ti(i))))
        n = filter(>(0), vec(view(counts, Ti(i))))
        if isempty(slice)
            @printf("%d: no pixel met the %d-observation minimum\n", year(t), MIN_SAMPLES)
        else
            @printf("%d: mean %.3f  min %.3f  max %.3f  (%d px, median %d valid obs/px)\n",
                    year(t), mean(slice), minimum(slice), maximum(slice),
                    length(slice), isempty(n) ? 0 : round(Int, median(n)))
        end
    end
    println("\nWrote: $(out_path)  (layers: albedo_low, n_valid)")
    return albedo, counts
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
