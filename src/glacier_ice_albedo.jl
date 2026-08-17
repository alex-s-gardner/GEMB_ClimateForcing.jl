"""
Glacier bare-ice albedo from the C3S 300 m satellite albedo product.

`compute_glacier_ice_albedo` reduces the 10-daily albedo record to, per pixel and
calendar year, the **mean of the darkest few percent of that year's valid
observations**. On a glacier surface the annual minimum-albedo state is bare ice: snow
has ablated away and the exposed ice is at its most-darkened. Averaging a low percentile
rather than taking the single minimum keeps the statistic from being set by one bad
retrieval.

This is the observational counterpart to GEMB's prescribed bare-ice albedo, which is
otherwise a tuned constant (~0.35–0.55 depending on region and dust/algae loading).

Built on [`satellite_albedo`](@ref), so the same caveats apply: data is *ordered* from
CDS via async jobs rather than read from a cloud store, a cold year takes tens of minutes,
and results are cached on disk so reruns are nearly free.
"""

"""
Calendar years covered by default: the Sentinel-3 300 m era. Earlier AVHRR/VGT/PROBA
eras of this product use different grids and are out of scope, so this is exactly the
range [`satellite_albedo`](@ref) accepts — aliased, not restated, so the two cannot drift.
"""
const GLACIER_ICE_ALBEDO_YEARS = _ALBEDO_YEARS

# Broadband directional-hemispherical (black-sky) albedo, and the full-spectrum layer
# within that file plus its companion uncertainty. Black-sky is the right choice for a
# surface energy-balance model; see `SATELLITE_ALBEDO_VARIABLES`.
const _ICE_ALBEDO_VARIABLE = :albb_dh
const _ICE_ALBEDO_LAYER = :AL_DH_BB
const _ICE_ALBEDO_ERR_LAYER = :AL_DH_BB_ERR

const _ICE_ALBEDO_PERCENTILE = 0.05   # fraction of darkest observations to average
const _ICE_ALBEDO_MIN_SAMPLES = 10    # min valid observations per pixel-year
const _ICE_ALBEDO_MAX_ERROR = 0.20    # reject obs whose _ERR exceeds this (absolute albedo)

# Accept-range for an individual observation. The 0.3 floor is a *glaciological* bound,
# not the physical one: exposed glacier ice rarely falls below ~0.3 broadband, so darker
# retrievals over a glacier pixel are usually rock, water, shadow or a failed inversion.
# Lower it (toward ~0.15) for heavily dust- or algae-darkened ablation zones, where 0.3
# would clip the very signal being measured.
const _ICE_ALBEDO_RANGE = (0.3, 1.0)

"""
QFLAG classes to reject, matched case-insensitively as substrings of the file's own
`flag_meanings`.

The v3.1 legend (19 bitmask classes, read from a delivered file) is **not** a land/cloud
mask — there is no sea, water, cloud or shadow class in it at all. It is:

    bit 0            snow_presence
    bits 1–9         no_obs_in_last_decade_for_<band>   (Oa03 Oa04 Oa07 Oa17 Oa21 S1 S2 S5 S6)
    bits 10–18       brdf_warning_for_<band>            (same nine bands)

So sea, cloud and shadow pixels arrive as `_FillValue`/`missing` instead, and are already
dropped by the missing-value check — QFLAG's job here is narrower.

`no_obs_in_last_decade_*` marks a band with no observation backing the retrieval in that
10-day window, which is exactly the cloud-persistent case that produces spurious low
albedo, so it is rejected by default. `snow_presence` is deliberately *kept* — a
snow-covered timestep is a bright observation that the low percentile discards on its own,
and rejecting it here would bias the sample count instead.

The generic names beyond the first entry are harmless when absent from this legend; they
exist so the same filter still bites on other product versions/eras.
"""
const GLACIER_ICE_ALBEDO_QFLAG_REJECT = [
    "no_obs_in_last_decade",
    "sea", "water", "ocean", "cloud", "shadow",
    "invalid", "fail", "error", "no_data", "nodata", "missing", "unusable", "unfilled",
]

## ------------------------------------------------------------------- QFLAG table

"""
    _qflag_table(path) -> Vector{@NamedTuple{meaning::String, mask::Int, is_value::Bool}}

The QFLAG legend as declared by the file itself: CF `flag_meanings` paired with either
`flag_masks` (bitwise) or `flag_values` (exact value). Read from the file rather than
hardcoded because the convention is version-specific. Returns an empty vector when the
attributes are absent, so callers can degrade to range-only QC instead of guessing bits.
"""
function _qflag_table(path::AbstractString)
    NCDataset(String(path), "r") do ds
        haskey(ds, "QFLAG") || return NamedTuple[]
        att = ds["QFLAG"].attrib
        meanings = get(att, "flag_meanings", nothing)
        isnothing(meanings) && return NamedTuple[]
        names = split(strip(String(meanings)))
        masks = get(att, "flag_masks", nothing)
        is_value = isnothing(masks)
        codes = is_value ? get(att, "flag_values", nothing) : masks
        isnothing(codes) && return NamedTuple[]
        codes = Int.(collect(codes))
        n = min(length(names), length(codes))
        return [(meaning=names[i], mask=codes[i], is_value=is_value) for i in 1:n]
    end
end

"""
    _qflag_rejects(table; patterns, brdf_warning=false) -> (bitmask, values, matched)

Split the legend into a single OR-ed `bitmask` (for `flag_masks` files) and a list of
exact `values` (for `flag_values` files), covering every class whose meaning matches one
of `patterns`. `matched` names the classes for logging — QC that silently drops data is
worse than none, and this product's legend does not contain the class names one would
expect (see [`GLACIER_ICE_ALBEDO_QFLAG_REJECT`](@ref)).

Set `brdf_warning=true` to additionally reject the nine `brdf_warning_*` classes. Off by
default: that is a warning rather than a failure, and on a snow/ice surface it fires often
enough to thin the sample materially.
"""
function _qflag_rejects(table; patterns=GLACIER_ICE_ALBEDO_QFLAG_REJECT,
                        brdf_warning::Bool=false)
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
    _albedo_qflag_legend(path; patterns, brdf_warning=false) -> (bitmask, values, matched)

Resolve the reject bitmask from a cached product file. Returns all-empty (so QC degrades
to range-only) rather than erroring if the legend cannot be read.
"""
function _albedo_qflag_legend(path; patterns=GLACIER_ICE_ALBEDO_QFLAG_REJECT,
                              brdf_warning::Bool=false)
    table = isnothing(path) ? NamedTuple[] : _qflag_table(path)
    if isempty(table)
        @warn """
        No CF QFLAG legend (`flag_masks`/`flag_meanings`) found in the product files, so \
        QFLAG-based rejection is skipped and only missing/out-of-range observations are \
        dropped.
        """
        return 0, Int[], String[]
    end
    return _qflag_rejects(table; patterns=patterns, brdf_warning=brdf_warning)
end

## ------------------------------------------------------------------ per-obs mask

_albedo_isbad(::Missing) = true
_albedo_isbad(x::Number) = !isfinite(x)

"""
    _valid_albedo(alb, qflag, err; bitmask, values, max_error, albedo_range) -> Matrix{Float32}

One timestep's albedo as a plain `Float32` matrix with every rejected observation set to
`NaN`. `qflag` and `err` may be `nothing` when those layers are unavailable.

Rejects, in order: missing/`_FillValue` pixels (NCDatasets surfaces these as `missing` or
`NaN`), albedo outside `albedo_range`, flagged QFLAG classes, and observations whose
companion uncertainty exceeds `max_error`. A *missing* uncertainty is not itself
disqualifying; an oversized one is.
"""
function _valid_albedo(alb::AbstractMatrix, qflag, err;
                       bitmask::Integer=0, values::AbstractVector{<:Integer}=Int[],
                       max_error::Union{Nothing,Real}=nothing,
                       albedo_range::Tuple{Real,Real}=_ICE_ALBEDO_RANGE)
    lo, hi = albedo_range
    out = Matrix{Float32}(undef, size(alb))
    @inbounds for i in eachindex(out, alb)
        a = alb[i]
        if _albedo_isbad(a) || a < lo || a > hi
            out[i] = NaN32
            continue
        end
        if !isnothing(qflag)
            q = qflag[i]
            if _albedo_isbad(q)
                out[i] = NaN32
                continue
            end
            qi = Int(q)
            if (bitmask != 0 && (qi & bitmask) != 0) || (!isempty(values) && qi in values)
                out[i] = NaN32
                continue
            end
        end
        if !isnothing(err) && !isnothing(max_error)
            e = err[i]
            if !_albedo_isbad(e) && e > max_error
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
    _LowPercentileTopK(sz, n_expected, percentile, min_samples)

Streaming accumulator for the darkest-percentile mean: retains, per pixel, only the
`kmax = ceil(percentile * n_expected)` smallest values seen so far, plus a count of valid
observations.

Streaming rather than holding the year's masked grids is what keeps this affordable at
scale. The statistic needs the `k = ceil(percentile * n_valid)` smallest values, and
`k ≤ kmax` because `n_valid ≤ n_expected` — so the retained set always contains them. At
the default 5 % over ~36 timesteps that is 2 retained values per pixel instead of 36:
49 → 5.5 MiB at 600×600, and the gap widens with the extent. Kernel time is a wash
(~45 ms either way over that grid); this buys memory, not speed.

The result is bit-identical to sorting every valid value and averaging the first `k`: the
same `k` values, in the same ascending order, reach the same `mean`.
"""
mutable struct _LowPercentileTopK
    percentile::Float64
    min_samples::Int
    kmax::Int
    n_expected::Int
    n_seen::Int
    sz::Tuple{Int,Int}
    # Field types are abstract so the same struct can hold either plain in-memory arrays
    # (the default) or `Mmap.mmap` views backed by scratch files (`scratch_dir`), which is
    # what makes a global grid possible at all — see `_topk_scratch_arrays`. The kernels
    # below are the hot path, so they take these as typed arguments rather than reading the
    # fields repeatedly.
    vals::AbstractMatrix{Float32}    # kmax × npixel, the retained smallest values (unsorted)
    nfilled::AbstractVector{Int32}   # how many of the kmax slots are in use
    counts::AbstractVector{Int32}    # valid observations per pixel
end

"""
    _topk_scratch_arrays(dir, kmax, npixel) -> (vals, nfilled, counts)

Allocate the accumulator's three arrays as `Mmap.mmap` views over files in `dir`.

At global 300 m (120960 × 60480 = 7.32 Gpx) the in-memory form needs
`kmax·npixel·4 + 2·npixel·4` ≈ 117 GB, which does not fit in RAM. Paging the same bytes
through the filesystem trades that for ~95 GB of scratch disk and lets the OS keep only
the working set resident. Values stay `Float32` and counts `Int32`, so results are
bit-identical to the in-memory path — this is a *placement* change, not a numeric one.

`nfilled`/`counts` are zero-filled because `Mmap.mmap` guarantees a new file reads as
zeros, which is exactly the required initial state; `vals` is written before it is read
(guarded by `nfilled`), so it is left undefined.
"""
function _topk_scratch_arrays(dir::AbstractString, kmax::Int, npixel::Int)
    mkpath(dir)
    open_map(name, T, dims) = open(joinpath(dir, name), "w+") do io
        Mmap.mmap(io, Array{T,length(dims)}, dims)
    end
    vals = open_map("topk_vals.bin", Float32, (kmax, npixel))
    nfilled = open_map("topk_nfilled.bin", Int32, (npixel,))
    counts = open_map("topk_counts.bin", Int32, (npixel,))
    return vals, nfilled, counts
end

function _LowPercentileTopK(sz::Tuple{Int,Int}, n_expected::Integer,
                            percentile::Real, min_samples::Integer;
                            scratch_dir::Union{Nothing,AbstractString}=nothing)
    (0 < percentile <= 1) || throw(ArgumentError(
        "percentile must be in (0, 1]; got $(percentile)"))
    n_expected >= 1 || throw(ArgumentError("no observations supplied"))
    kmax = max(1, ceil(Int, percentile * n_expected))
    npixel = prod(sz)
    vals, nfilled, counts = if isnothing(scratch_dir)
        Matrix{Float32}(undef, kmax, npixel), zeros(Int32, npixel), zeros(Int32, npixel)
    else
        _topk_scratch_arrays(scratch_dir, kmax, npixel)
    end
    return _LowPercentileTopK(Float64(percentile), Int(min_samples), kmax, Int(n_expected),
                              0, sz, vals, nfilled, counts)
end

"""
    _accumulate!(acc, m)

Fold one timestep's masked albedo grid (invalid entries `NaN`) into `acc`.
"""
function _accumulate!(acc::_LowPercentileTopK, m::AbstractMatrix{Float32})
    size(m) == acc.sz || throw(ArgumentError(
        "observation grid $(size(m)) does not match $(acc.sz)"))
    acc.n_seen += 1
    acc.n_seen <= acc.n_expected || throw(ArgumentError(
        "more observations accumulated than the $(acc.n_expected) declared"))
    # Offset 0: a whole-grid fold writes pixel p of the accumulator from element p of `m`.
    _accumulate_block!(acc.vals, acc.nfilled, acc.counts, m, 0, acc.kmax)
    return acc
end

"""
    _accumulate_block!(vals, nfilled, counts, m, offset, kmax)

Fold `m` into the accumulator arrays starting at linear pixel `offset`, i.e. element `i`
of `m` updates pixel `offset + i`.

Split out of [`_accumulate!`](@ref) for two reasons. It takes the three arrays as
*concretely typed* parameters, so the loop specializes on whether they are `Array`s or
`Mmap` views instead of dynamically dispatching every element access through
`_LowPercentileTopK`'s abstract fields — that indirection would otherwise cost far more
than the mmap paging it enables. And the `offset` lets a caller fold one Y-block of a
grid at a time ([`_accumulate_blockwise!`](@ref)) without materialising the whole grid.
"""
function _accumulate_block!(vals::AbstractMatrix{Float32}, nfilled::AbstractVector{Int32},
                            counts::AbstractVector{Int32}, m::AbstractArray{Float32},
                            offset::Int, kmax::Int)
    @inbounds for i in eachindex(m)
        v = m[i]
        isnan(v) && continue
        p = offset + i
        counts[p] += Int32(1)
        nf = Int(nfilled[p])
        if nf < kmax
            vals[nf + 1, p] = v
            nfilled[p] = Int32(nf + 1)
        else
            # Evict the largest of the retained values, if this one is smaller. kmax is
            # tiny (2 at the defaults), so a linear scan beats any heap bookkeeping.
            worst, worst_v = 1, vals[1, p]
            for j in 2:kmax
                if vals[j, p] > worst_v
                    worst, worst_v = j, vals[j, p]
                end
            end
            v < worst_v && (vals[worst, p] = v)
        end
    end
    return nothing
end

"""
    _finalize(acc) -> (albedo, counts)

Reduce the accumulator to the per-pixel statistic and the valid-observation count.

The number of observations averaged is `max(1, ceil(percentile * n_valid))`, so a pixel
with 36 valid scenes averages its darkest 2 at 5 %. Pixels with fewer than `min_samples`
valid observations yield `NaN` — a "darkest 5 %" drawn from a handful of scenes is noise —
but their count is still reported, so a caller can see *why* a pixel is missing.
"""
function _finalize(acc::_LowPercentileTopK)
    result = fill(NaN32, acc.sz)
    _finalize_into!(result, acc.vals, acc.nfilled, acc.counts,
                    acc.kmax, acc.percentile, acc.min_samples)
    # `reshape` of an mmap view stays an mmap view, so a global-grid count layer is not
    # copied into RAM here; `Raster` construction later decides what to materialise.
    return result, reshape(acc.counts, acc.sz)
end

"""
    _finalize_into!(result, vals, nfilled, counts, kmax, percentile, min_samples)

Write the per-pixel statistic into `result`, whose length is the number of pixels.

Typed-argument inner kernel for the same reason as [`_accumulate_block!`](@ref): reading
`_LowPercentileTopK`'s abstract fields elementwise would dynamically dispatch on every
access.
"""
function _finalize_into!(result::AbstractArray{Float32}, vals::AbstractMatrix{Float32},
                         nfilled::AbstractVector{Int32}, counts::AbstractVector{Int32},
                         kmax::Int, percentile::Float64, min_samples::Int)
    buf = Vector{Float32}(undef, kmax)
    @inbounds for p in eachindex(result)
        n = Int(counts[p])
        n >= min_samples || continue
        nf = Int(nfilled[p])
        nf >= 1 || continue
        k = min(max(1, ceil(Int, percentile * n)), nf)
        for j in 1:nf
            buf[j] = vals[j, p]
        end
        kept = view(buf, 1:nf)
        sort!(kept)
        result[p] = mean(view(kept, 1:k))
    end
    return result
end

"""
    _accumulate_blockwise!(acc, alb, qflag, err; block_rows, qc...) -> n_valid

Read one timestep in horizontal blocks of `block_rows` grid rows and fold each block into
`acc`, returning the number of observations that survived QC.

**Why blocks.** `_accumulate!`'s caller normally does `parent(read(layer))`, materialising
the whole grid. At global 300 m one layer is 7.32 Gpx — 37 GB as
`Union{Missing,Float32}`, and three layers (albedo, `QFLAG`, `_ERR`) is ~110 GB, which
exceeds RAM before any statistic is formed. Reading a slab of rows at a time bounds that
by `block_rows / ny`: at the default 512 rows a global block is ~0.9 GB across all three
layers, and a small extent is read in a single block, so nothing changes for regional runs.

Blocks are taken along the **second** (Y) dimension because the product's chunk layout is
`(nx, ny, 1)` — full rows of X — so a Y-slab is contiguous in the file and reads without
random access. The linear pixel offset works because Julia is column-major: block `j0:j1`
of an `(nx, ny)` grid is exactly linear pixels `(j0-1)·nx + 1 … j1·nx`, contiguous and in
the same order the accumulator indexes.
"""
function _accumulate_blockwise!(acc::_LowPercentileTopK, alb, qflag, err;
                                block_rows::Integer=512,
                                bitmask::Integer=0,
                                values::AbstractVector{<:Integer}=Int[],
                                max_error::Union{Nothing,Real}=nothing,
                                albedo_range::Tuple{Real,Real}=_ICE_ALBEDO_RANGE)
    nx, ny = size(alb)
    (nx, ny) == acc.sz || throw(ArgumentError(
        "observation grid $((nx, ny)) does not match $(acc.sz)"))
    acc.n_seen += 1
    acc.n_seen <= acc.n_expected || throw(ArgumentError(
        "more observations accumulated than the $(acc.n_expected) declared"))

    n_valid = 0
    step = max(1, Int(block_rows))
    for j0 in 1:step:ny
        j1 = min(j0 + step - 1, ny)
        rows = j0:j1
        a = alb[:, rows]
        q = isnothing(qflag) ? nothing : qflag[:, rows]
        e = isnothing(err) ? nothing : err[:, rows]
        m = _valid_albedo(a, q, e; bitmask=bitmask, values=values,
                          max_error=max_error, albedo_range=albedo_range)
        n_valid += count(!isnan, m)
        _accumulate_block!(acc.vals, acc.nfilled, acc.counts, m, (j0 - 1) * nx, acc.kmax)
    end
    return n_valid
end

"""
    _low_percentile_mean(obs, percentile, min_samples) -> (albedo, counts)

Per-pixel mean of the darkest `percentile` fraction of valid observations in `obs`
(a vector of same-sized matrices, invalid entries `NaN`).

Non-streaming convenience wrapper over [`_LowPercentileTopK`](@ref) for when the whole
set of grids is already in hand; the driver streams instead.
"""
function _low_percentile_mean(obs::Vector{<:AbstractMatrix{Float32}},
                              percentile::Real, min_samples::Integer)
    isempty(obs) && throw(ArgumentError("no observations supplied"))
    acc = _LowPercentileTopK(size(first(obs)), length(obs), percentile, min_samples)
    for m in obs
        _accumulate!(acc, m)
    end
    return _finalize(acc)
end

## ------------------------------------------------------------------------ driver

"""
    compute_glacier_ice_albedo(years=GLACIER_ICE_ALBEDO_YEARS; extent=nothing, kwargs...)
        -> RasterStack

Bare-ice albedo per pixel and calendar year: the **mean of the darkest `percentile`
fraction of that year's valid broadband albedo observations** from the C3S 300 m product.

Returns a `RasterStack` over (`X`, `Y`, `Ti`) with two layers:
- `:glacier_ice_albedo` — the statistic (`Float32`, dimensionless 0–1, `NaN` where
  fewer than `min_samples` observations survived QC)
- `:n_valid_observations` — how many observations that pixel-year's statistic drew on
  (`Int32`), so a marginal pixel can be identified after the fact

The `Ti` lookup is `Date(year, 1, 1)` for each requested year.

# Why this statistic
On a glacier the annual albedo minimum *is* the bare-ice state: seasonal snow has
ablated, exposing ice at its most darkened by dust, black carbon and algae. The 10-daily
product gives ~36 observations a year, so 5 % averages the darkest 1–2 of them — low
enough to isolate the ice surface, but an average rather than a single scene, so one bad
retrieval cannot set the answer. GEMB otherwise takes bare-ice albedo as a tuned regional
constant; this is the observational substitute.

# Quality control
Applied per observation, before any statistic is formed:
1. missing / `_FillValue` pixels (this is what removes sea, cloud and shadow — see
   [`GLACIER_ICE_ALBEDO_QFLAG_REJECT`](@ref), the QFLAG legend is *not* a cloud mask)
2. albedo outside `albedo_range`
3. QFLAG classes matching `qflag_reject`, read from each file's own `flag_masks` /
   `flag_meanings` attributes rather than hardcoded
4. observations whose companion `_ERR` uncertainty exceeds `max_error`
5. pixel-years left `NaN` below `min_samples` surviving observations

# Arguments
- `years`: any iterable of integers, or a single `Integer` — `2019`, `2019:2021`,
  `[2018, 2022]`. Restricted to the Sentinel-3 300 m era (2018–2024).

# Keywords
- `extent`: an `Extents.Extent` or `(; X, Y)` NamedTuple in −180…180 longitude. Defaults
  to `nothing` (global), which is ~120960 × 47040 pixels *per timestep* — pass an extent.
- `percentile = $(_ICE_ALBEDO_PERCENTILE)`: fraction of darkest observations to average.
- `min_samples = $(_ICE_ALBEDO_MIN_SAMPLES)`: minimum valid observations per pixel-year.
- `max_error = $(_ICE_ALBEDO_MAX_ERROR)`: reject observations whose `_ERR` exceeds this;
  `nothing` skips the check.
- `albedo_range = $(_ICE_ALBEDO_RANGE)`: accept-range per observation. The 0.3 floor is
  glaciological, not physical — lower it for heavily darkened ablation zones.
- `variable`, `layer`, `err_layer`: which product variable and NetCDF layers to read.
  Defaults are broadband black-sky (`:albb_dh` / `:AL_DH_BB`), with `err_layer` derived
  from `layer` as `<layer>_ERR`; set `err_layer=nothing` to skip uncertainty QC.
- `qflag_reject`, `qflag_reject_brdf_warning`: QFLAG filtering, see above.
- `token`: CDS API key; defaults to [`get_cds_api_key()`](@ref).
- `cache_path`: product-file cache. Defaults to `satellite_albedo`'s per-user temp
  directory — pass a durable path, since a lost cache means reordering hours of CDS jobs.
- `timeout`, `max_concurrent_jobs`, `verbose`: passed through to
  [`satellite_albedo`](@ref), whose defaults apply.

## Very large grids
Four keywords exist only for grids far larger than a glacier basin. All are no-ops or
off by default, so a regional run behaves exactly as before:
- `block_rows = 512`: read each timestep in Y-slabs of this many rows rather than whole,
  so a global timestep (~67 GB across its three layers) never materialises. Results are
  bit-identical to a whole-grid read at any block size.
- `scratch_dir`, `scratch_threshold_pixels = 200_000_000`: past that pixel count the
  accumulator (`kmax·npx·4 + 2·npx·4` bytes — ~117 GB globally) is backed by `Mmap`
  scratch files instead of RAM. Pass `scratch_dir` to force it, and to control *where*.
- `batch_timesteps`: order the year in batches of this many timesteps instead of one
  call. **This is a deliberate pessimisation** — batches are ordered serially, so
  wall-clock grows with the batch count, whereas the default single call orders the
  year's ~3 jobs concurrently. Use it only with:
- `discard_after_fold = false`: **delete each batch's product files once folded.**
  Destructive: it gives up the cache, so a re-run reorders from CDS. Together these two
  bound peak disk to one batch, which is the only way a global year (~36 × 67 GB, far
  past any normal free space) fits at all.

# Cost
A cold year is ~3 CDS jobs, ordered and polled concurrently, so budget **30–90 minutes
per year**; cached years are nearly free. The whole year is requested in one
`satellite_albedo` call on purpose — looping over months by hand re-serialises the
ordering and is much slower.

# Prerequisite
A CDS API key, plus one-time licence acceptance at
<https://cds.climate.copernicus.eu/datasets/satellite-albedo?tab=download#manage-licences>.

# Example
```julia
using Dates, Rasters

ice = compute_glacier_ice_albedo(2019:2020;
    extent = Extent(X = (-50.0, -45.0), Y = (66.0, 68.0)))

ice[:glacier_ice_albedo]                      # (X, Y, Ti) Float32, NaN where unresolved
mean(filter(!isnan, vec(ice[:glacier_ice_albedo])))

# Multi-year bare-ice albedo for GEMB, one value per pixel:
using Statistics
albedo_ice = map(eachslice(ice[:glacier_ice_albedo]; dims=(X, Y))) do px
    v = filter(!isnan, collect(px))
    isempty(v) ? NaN32 : mean(v)
end
```

# See also
[`satellite_albedo`](@ref) for the underlying loader and its layer/variable options.
"""
compute_glacier_ice_albedo(year::Integer; kwargs...) =
    compute_glacier_ice_albedo(year:year; kwargs...)

function compute_glacier_ice_albedo(years=GLACIER_ICE_ALBEDO_YEARS;
                                    extent=nothing,
                                    variable::Symbol=_ICE_ALBEDO_VARIABLE,
                                    layer::Symbol=_ICE_ALBEDO_LAYER,
                                    err_layer::Union{Nothing,Symbol}=Symbol(layer, :_ERR),
                                    percentile::Real=_ICE_ALBEDO_PERCENTILE,
                                    min_samples::Integer=_ICE_ALBEDO_MIN_SAMPLES,
                                    max_error::Union{Nothing,Real}=_ICE_ALBEDO_MAX_ERROR,
                                    albedo_range::Tuple{Real,Real}=_ICE_ALBEDO_RANGE,
                                    qflag_reject=GLACIER_ICE_ALBEDO_QFLAG_REJECT,
                                    qflag_reject_brdf_warning::Bool=false,
                                    token::Union{Nothing,String}=nothing,
                                    cache_path::Union{Nothing,String}=nothing,
                                    verbose::Bool=true,
                                    # Large-grid controls. Both are no-ops at regional size:
                                    # a small grid is read in one block and kept in RAM, so
                                    # results and cost are unchanged from before.
                                    block_rows::Integer=512,
                                    scratch_dir::Union{Nothing,AbstractString}=nothing,
                                    scratch_threshold_pixels::Real=200_000_000,
                                    # Disk-constrained ordering. `nothing` = one call for
                                    # the whole year, the concurrent fast path.
                                    batch_timesteps::Union{Nothing,Integer}=nothing,
                                    discard_after_fold::Bool=false,
                                    # `timeout`, `max_concurrent_jobs`, `poll_interval`,
                                    # `force_download`, … go straight to the loader rather
                                    # than being restated here with defaults that could drift.
                                    albedo_kwargs...)
    year_list = collect(years)
    isempty(year_list) && throw(ArgumentError("no years requested"))
    (0 < percentile <= 1) || throw(ArgumentError(
        "percentile must be in (0, 1]; got $(percentile)"))
    min_samples >= 1 || throw(ArgumentError(
        "min_samples must be at least 1; got $(min_samples)"))
    first(albedo_range) < last(albedo_range) || throw(ArgumentError(
        "albedo_range must be (low, high) with low < high; got $(albedo_range)"))

    annual = Matrix{Float32}[]
    annual_counts = Matrix{Int32}[]
    xy = nothing
    # The QFLAG legend is a property of the product version, not of the year, so it is
    # read from the first year's files and reused.
    legend = nothing

    want = Symbol[layer, :QFLAG]
    isnothing(err_layer) || push!(want, err_layer)

    for y in year_list
        verbose && @info "Loading $(y) albedo (ordering from CDS if not cached)"

        # Every timestep the year *should* yield, known before any ordering, so the
        # accumulator's `n_expected` (hence `kmax`) is right even when the year is
        # ordered in several batches.
        year_dates = _albedo_timesteps((DateTime(y, 1, 1), DateTime(y, 12, 31, 23, 59, 59)))
        isempty(year_dates) && error("No albedo timesteps fall in $(y)")

        # One batch = one `satellite_albedo` call. The default is a single batch for the
        # whole year, which is the fast path: that call splits the year into the ~3 jobs
        # the cost limit allows and runs them **concurrently**, so the year costs about
        # one job's queue time rather than three in series.
        #
        # `batch_timesteps` exists only for grids whose *product files* do not fit on
        # disk — a global year is ~36 × 67 GB — and it is a deliberate pessimisation:
        # batches are ordered one after another, so wall-clock grows roughly with the
        # batch count. Do not set it to work around anything else.
        batches = isnothing(batch_timesteps) ? [year_dates] :
                  [year_dates[i:min(i + Int(batch_timesteps) - 1, end)]
                   for i in 1:Int(batch_timesteps):length(year_dates)]

        # Two independent scale problems, both only biting on very large grids, so both are
        # switched on by pixel count rather than changing behaviour for regional runs:
        #  * the accumulator itself is `kmax·npx·4 + 2·npx·4` bytes — 117 GB globally — so
        #    past `scratch_threshold_pixels` it is mmap-backed on disk instead of in RAM.
        #  * each timestep's three layers are ~110 GB globally if read whole, so they are
        #    read and folded in Y-blocks (`_accumulate_blockwise!`).
        acc = nothing
        scratch = nothing
        n_folded = 0

        for (bi, batch) in enumerate(batches)
            length(batches) > 1 && verbose &&
                @info "  batch $(bi)/$(length(batches)): $(first(batch)) … $(last(batch))" n_timesteps=length(batch)
            time_range = (DateTime(first(batch)),
                          DateTime(last(batch)) + Hour(23) + Minute(59) + Second(59))

            # All three layers in one call: they live in the same product files, so the
            # extra layers cost no extra CDS job, and one call keeps the timesteps
            # aligned by construction.
            series = try
                satellite_albedo(; time_range, extent, variable, layer=want, verbose,
                                 token, cache_path, albedo_kwargs...)
            catch err
                # ONLY a layer-availability failure may fall through to the retry. That
                # error is raised by `_albedo_resolve_layer` *after* the CDS jobs have
                # completed and their files are cached, so re-calling `satellite_albedo`
                # re-reads the cache and orders nothing — which is what makes the
                # fallback cheap.
                #
                # Every other failure must propagate. A job timeout, a queue-limit give-up
                # or an HTTP error all arrive here mid-order, with nothing cached;
                # retrying then re-submits the *same* request (`layer` is not part of the
                # request hash, so it reuses the same job dir, finds no archive and orders
                # again). That turned a 3 h timeout into an endless 3-hourly resubmit loop
                # that never cached a byte and steadily filled the account's job queue.
                # `_cds_wait` raises timeouts and rejections via `error(...)` →
                # `ErrorException`, while the layer error is an `ArgumentError`, so the
                # type alone separates them.
                err isa ArgumentError || rethrow()
                @warn """
                Could not open all of $(want) — falling back to the albedo layer alone, so \
                QC degrades to missing/out-of-range only.
                """ exception = err
                satellite_albedo(; time_range, extent, variable, layer, verbose,
                                 token, cache_path, albedo_kwargs...)
            end
            dates = lookup(series, Ti)
            isempty(dates) && error("No albedo timesteps returned for $(first(batch)) … $(last(batch))")

            # A single surviving layer comes back as a series of `Raster`s rather than
            # `RasterStack`s, so name the accessors once instead of branching per timestep.
            multi = first(series) isa AbstractRasterStack
            _layer_of(st, l) = multi ? st[l] : st
            has(l) = multi && haskey(first(series), l)
            qflag_key = has(:QFLAG) ? :QFLAG : nothing
            err_key = (!isnothing(err_layer) && has(err_layer)) ? err_layer : nothing

            if isnothing(legend)
                path = isnothing(qflag_key) ? nothing :
                       get(metadata(_layer_of(first(series), qflag_key)), "file", nothing)
                legend = isnothing(qflag_key) ? (0, Int[], String[]) :
                         _albedo_qflag_legend(path; patterns=qflag_reject,
                                              brdf_warning=qflag_reject_brdf_warning)
                verbose && !isempty(legend[3]) && @info "Rejecting QFLAG classes" classes =
                    legend[3] bitmask = legend[1] values = legend[2]
            end
            bitmask, values, _ = legend

            # Stream the timesteps into the accumulator rather than holding the batch's
            # masked grids: only the darkest `ceil(percentile * n)` values per pixel are
            # ever needed.
            spent = String[]
            for (i, d) in enumerate(dates)
                st = series[i]
                alb = _layer_of(st, layer)
                if isnothing(acc)
                    sz = size(alb)
                    npx = prod(sz)
                    use_scratch = !isnothing(scratch_dir) ||
                                  npx > scratch_threshold_pixels
                    if use_scratch
                        scratch = isnothing(scratch_dir) ?
                                  mktempdir(; prefix="gemb_ice_albedo_") : String(scratch_dir)
                        kmax = max(1, ceil(Int, percentile * length(year_dates)))
                        gib = npx * (kmax * 4 + 8) / 2^30
                        verbose && @info "Large grid: backing the accumulator with scratch files" pixels=npx scratch_gib=round(gib; digits=1) dir=scratch
                    end
                    acc = _LowPercentileTopK(sz, length(year_dates), percentile, min_samples;
                                             scratch_dir=use_scratch ? scratch : nothing)
                end
                # Lazy `Raster`s: index the block rather than `read`ing the whole layer, so
                # a global timestep never materialises in full.
                q = isnothing(qflag_key) ? nothing : st[qflag_key]
                e = isnothing(err_key) ? nothing : st[err_key]
                nvalid = _accumulate_blockwise!(acc, alb, q, e; block_rows=block_rows,
                                                bitmask, values, max_error, albedo_range)
                n_folded += 1
                if verbose
                    @info "  $(d)  valid $(round(100 * nvalid / prod(size(alb)); digits=1))%"
                end
                isnothing(xy) && (xy = dims(alb, (X, Y)))
                # The backing file is recorded per layer by the loader; collect it now,
                # while the layer is still open, and delete only after the whole batch has
                # been folded (all three layers share one file).
                discard_after_fold && for l in (layer, qflag_key, err_key)
                    isnothing(l) && continue
                    p = get(metadata(_layer_of(st, l)), "file", nothing)
                    isnothing(p) || push!(spent, String(p))
                end
            end

            # Free the batch's product files before ordering the next one. Destructive and
            # opt-in: it trades the cache (and so any cheap re-run) for peak disk, which is
            # the only way a global year fits at all.
            if discard_after_fold
                series = nothing
                GC.gc()   # NCDatasets holds the file open until the handles are collected
                freed = 0
                for p in unique(spent)
                    isfile(p) || continue
                    freed += filesize(p)
                    rm(p; force=true)
                end
                verbose && @info "  discarded batch product files" n=length(unique(spent)) freed_gib=round(freed / 2^30; digits=2)
            end
        end

        n_folded == length(year_dates) || @warn "Fewer timesteps folded than $(y) should have" folded=n_folded expected=length(year_dates)

        stat, counts = _finalize(acc)
        push!(annual, stat)
        push!(annual_counts, counts)
        if verbose
            v = filter(!isnan, vec(stat))
            @info "$(y): $(length(v))/$(length(stat)) pixels resolved, mean of darkest " *
                  "$(round(100 * percentile; digits=1))% = " *
                  "$(isempty(v) ? NaN : round(mean(v); digits=3))"
        end
    end

    ti = Ti(Date.(year_list, 1, 1))
    meta = Dict{String,Any}(
        "statistic" => "mean of lowest $(round(100 * percentile; digits=2))% of valid " *
                       "observations per calendar year",
        "interpretation" => "glacier bare-ice albedo (annual minimum-albedo surface state)",
        "variable" => variable,
        "layer" => layer,
        "err_layer" => err_layer,
        "percentile" => percentile,
        "min_samples" => min_samples,
        "max_error" => max_error,
        "albedo_range" => albedo_range,
        "qflag_reject" => collect(qflag_reject),
        "qflag_reject_brdf_warning" => qflag_reject_brdf_warning,
        "qflag_rejected_classes" => legend[3],
        "units" => "1 (dimensionless reflectance fraction, 0-1)",
        "source" => "C3S satellite surface albedo, Sentinel-3 OLCI+SLSTR 300 m, v3_1",
    )

    albedo = Raster(cat(annual...; dims=3), (xy..., ti);
                    name=:glacier_ice_albedo, missingval=NaN32, metadata=meta)
    # Same provenance, different units — a count is not a reflectance.
    count_meta = merge(meta, Dict{String,Any}(
        "units" => "1 (count of valid observations)",
        "statistic" => "number of observations surviving QC per calendar year"))
    counts = Raster(cat(annual_counts...; dims=3), (xy..., ti);
                    name=:n_valid_observations, missingval=Int32(0), metadata=count_meta)
    return RasterStack((; glacier_ice_albedo=albedo, n_valid_observations=counts);
                       metadata=meta)
end
