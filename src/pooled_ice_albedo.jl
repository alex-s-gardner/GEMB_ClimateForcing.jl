"""
Bare-ice albedo pooled over the whole record, re-folded from the per-date sample cache.

`pool_ice_albedo_from_cache` answers a different question from
[`compute_glacier_ice_albedo_modis`](@ref) and from [`rgi7_ice_albedo_climatology`](@ref):
**one** darkest-`percentile` mean over every valid retrieval from every date, rather than a
per-year statistic or a reduction of per-year statistics.

That distinction is the whole point of the file. Coverage varies enormously between years — the
fraction of cells resolving in a year runs from 21.6 % (2001) down to 6.5 % (2025) in the north
and 0.6–6.2 % in the south — so a per-year value is missing for most cell-years, and a mean or
median *of* those values weights a sparse year equally with a well-covered one. Pooling asks
only how dark the cell got across the record, and reports how many retrievals said so.

Nothing here downloads. The cache holds raw DNs and the QA class (not post-QC values), so
`percentile`, `albedo_range` and `qa_keep` can all be changed by a re-fold; only the date list,
cell list and layer set are fixed at download time. A missing date is an error rather than a
silent gap, because a pooled statistic computed over an unknown subset of the record is not
interpretable.
"""

# Lower QC bound of 0.25 rather than the 0.30 the per-year run used. The floor is
# glaciological, not physical: it exists to reject rock, water, shadow and failed inversions,
# and every per-year file written at 0.30 has `minimum == exactly 0.3000`, i.e. the floor was
# clipping the dark tail it was meant to bound. Dust- and algae-darkened ablation zones sit
# below 0.30.
const POOLED_ICE_ALBEDO_RANGE = (0.25, 1.0)

# Same 5 % as the per-year path. Over a pooled record it is far more robust, not less: 5 % of a
# cell's ~355 lifetime retrievals averages ~18 values, against the darkest 1–2 in one year.
const _POOLED_ICE_PERCENTILE = 0.05

"""
    _pooled_layer_names(layer) -> (albedo, n_valid, k_used)

Output layer names for one MCD43A3 albedo layer.

The first two come from [`_modis_ice_albedo_layer_names`](@ref) so a pooled stack is keyed the
same way `compute_glacier_ice_albedo_modis` keys its per-year one; the third mirrors the second.
"""
function _pooled_layer_names(layer::Symbol)
    alb, nval = _modis_ice_albedo_layer_names(layer)
    return (alb, nval, Symbol("k_used", chopprefix(String(nval), "n_valid_observations")))
end

"""
    _pooled_k_used(n, percentile) -> Int

How many retained values `_finalize_into!` averages for a cell with `n` valid retrievals.

Derived rather than returned from the accumulator, which is exact: that function uses
`k = min(max(1, ceil(percentile·n)), nfilled)` with `nfilled == min(n, kmax)` and
`kmax == ceil(percentile·n_expected)`. Since `n ≤ n_expected`, `ceil(percentile·n) ≤ kmax`, and
since `percentile ≤ 1`, `ceil(percentile·n) ≤ n` — so both `min` arms are inactive and
`k == max(1, ceil(percentile·n))` whenever `n ≥ 1`.
"""
_pooled_k_used(n::Integer, percentile::Real) =
    n < 1 ? 0 : max(1, ceil(Int, percentile * n))

"""
    _pooled_cached_dates(cache_path, cells_key) -> Vector{Date}

Every date cached for `cells_key`, ascending.

The date list is taken from the cache rather than re-derived from a `doy_range`, so the pooled
statistic covers exactly what is on disk and `n_expected` is exact. Re-deriving it would risk
declaring more dates than are cached, which inflates `kmax` and silently retains more values
than the percentile calls for.
"""
function _pooled_cached_dates(cache_path::AbstractString, cells_key::AbstractString)
    dir = joinpath(cache_path, "samples", cells_key)
    isdir(dir) || throw(ArgumentError(
        "no sample cache at $(dir). A pooled re-fold reads cached samples only; run " *
        "`compute_glacier_ice_albedo_modis` (or data/run_rgi7_ice_albedo_modis.jl) first to " *
        "populate it."))
    dates = Date[]
    for f in readdir(dir)
        endswith(f, ".tsv") || continue
        push!(dates, Date(chop(f; tail = 4)))
    end
    isempty(dates) && throw(ArgumentError("sample cache at $(dir) holds no .tsv files"))
    sort!(dates)
    return dates
end

"""
    pool_ice_albedo_from_cache(cells; kwargs...) -> DimStack

Bare-ice albedo per MODIS cell, pooled over every cached date.

For each cell, the mean of the darkest `percentile` of *all* valid retrievals in the record —
not a per-year value and not a reduction of per-year values. Returns a `DimStack` over
`Dim{:point}(1:length(cells))`, in the order `cells` was given.

# Keywords
- `cache_path = _default_modis_cache()`: the granule/sample cache root — the same default the
  per-year MODIS path uses, so the two cannot disagree about where samples live. Set
  `ENV["GEMB_CACHE_PATH"]` (see [`_gemb_cache_root`](@ref)) or pass this explicitly.
- `percentile = $(_POOLED_ICE_PERCENTILE)`: fraction of the darkest retrievals averaged.
- `albedo_range = POOLED_ICE_ALBEDO_RANGE`: `(0.25, 1.0)`. See that constant for why the floor
  differs from the per-year run's 0.30.
- `qa_keep = MCD43A3_QA_KEEP`: QA classes accepted — every full BRDF inversion.
- `layer = MCD43A3_ALBEDO_LAYERS`: which albedo layers to reduce.
- `dates = nothing`: restrict to a subset of the cached dates. `nothing` uses all of them,
  which is what "across all years" means.
- `progress = true`: report a plain-text bar per folded date.

# Returns
Layers `:albedo_bsa` / `:albedo_wsa` (`NaN32` where no retrieval survived QC), each with
`:n_valid_*` (retrievals passing QA and the range) and `:k_used_*` (how many of the darkest
were averaged), plus `:latitude`, `:longitude` and `:cell_id`.

# Example
```julia
cells = rgi7_modis_unique_cells(rgi7_modis_cells())
north = filter(c -> c[2] <= 8, cells)
ice = pool_ice_albedo_from_cache(north)
ice[:albedo_bsa]     # one value per cell, over the whole record
ice[:n_valid_bsa]    # how many retrievals stand behind it
```

!!! note "Memory scales with `percentile × n_dates`"
    The accumulator retains `kmax = ceil(percentile · n_dates)` values per cell — 223 for a
    26-year northern melt-season record at 5 %, against 9 for a single year. At 2.58 M cells
    and two layers that is 4.3 GiB, held in RAM. A larger `percentile` or a longer record
    scales it linearly.
"""
function pool_ice_albedo_from_cache(cells::AbstractVector{<:NTuple{4,Integer}};
                                    cache_path::AbstractString = _default_modis_cache(),
                                    percentile::Real = _POOLED_ICE_PERCENTILE,
                                    albedo_range::Tuple{Real,Real} = POOLED_ICE_ALBEDO_RANGE,
                                    qa_keep::AbstractVector{<:Integer} = MCD43A3_QA_KEEP,
                                    layer = MCD43A3_ALBEDO_LAYERS,
                                    dates::Union{Nothing,AbstractVector{Date}} = nothing,
                                    progress::Bool = true)
    isempty(cells) && throw(ArgumentError("no cells supplied"))
    (0 < percentile <= 1) || throw(ArgumentError(
        "percentile must be in (0, 1]; got $(percentile)"))
    albedo_range[1] < albedo_range[2] || throw(ArgumentError(
        "albedo_range must be (lo, hi) with lo < hi; got $(albedo_range)"))
    isempty(qa_keep) && throw(ArgumentError(
        "qa_keep is empty, so every observation would be rejected; pass at least " *
        "one QA class (see MCD43A3_QA_KEEP)."))
    layer_list = layer isa Symbol ? [layer] : collect(layer)
    for l in layer_list
        l in MCD43A3_ALBEDO_LAYERS || throw(ArgumentError(
            "$(l) is not an albedo layer; expected one of $(MCD43A3_ALBEDO_LAYERS)"))
    end

    cell_list = [(Int(c[1]), Int(c[2]), Int(c[3]), Int(c[4])) for c in cells]
    # The cache is keyed by a SHA of the *whole* cell list in the order given, so the list must
    # match the one the samples were downloaded for, bit for bit.
    cells_key = _modis_cells_key(cell_list)
    ncell = length(cell_list)

    all_dates = _pooled_cached_dates(cache_path, cells_key)
    date_list = isnothing(dates) ? all_dates : sort(collect(dates))
    isempty(date_list) && throw(ArgumentError("no dates to fold"))
    missing_dates = setdiff(date_list, all_dates)
    isempty(missing_dates) || throw(ArgumentError(
        "$(length(missing_dates)) requested dates are not cached (first: " *
        "$(first(missing_dates))). A pooled statistic over an unknown subset of the record " *
        "is not interpretable, so this is an error rather than a gap."))

    quality_of = Dict(l => _mcd43a3_quality_layer(l) for l in layer_list)
    cache_layers = unique(vcat(layer_list, collect(values(quality_of))))

    # One accumulator for the entire record — this is what makes the statistic pooled rather
    # than per-year, and why `kmax` is sized from every date at once.
    accs = Dict(l => _LowPercentileTopK((ncell, 1), length(date_list), percentile, 1)
                for l in layer_list)
    if progress
        @info "Pooled bare-ice albedo re-fold (no download)" cells = ncell dates = length(date_list) kmax = accs[first(layer_list)].kmax percentile albedo_range accumulator_GiB = round(length(layer_list) * accs[first(layer_list)].kmax * ncell * 4 / 2^30; digits = 2)
    end

    t0 = time()
    for (i, date) in enumerate(date_list)
        path = _modis_sample_cache_path(cache_path, cells_key, date)
        samples = _modis_sample_cache_read(path, cache_layers, ncell)
        isnothing(samples) && throw(ArgumentError(
            "the cached samples for $(date) at $(path) could not be used — the file is " *
            "absent, empty, or was written for a different cell list or layer set. Re-run " *
            "the per-year driver for that date before pooling."))
        for l in layer_list
            masked = _valid_albedo(reshape(samples[l], ncell, 1),
                                   reshape(samples[quality_of[l]], ncell, 1), nothing;
                                   keep_values = qa_keep, scale = _MCD43A3_SCALE,
                                   albedo_range = albedo_range)
            _accumulate!(accs[l], masked)
        end
        progress && _report_progress("MCD43A3 pooled", i, length(date_list), t0)
    end

    point_dim = Dim{:point}(1:ncell)
    out = Pair{Symbol,Any}[]
    for l in layer_list
        cell_albedo, cell_counts = _finalize(accs[l])
        alb = Vector{Float32}(undef, ncell)
        nval = Vector{Int32}(undef, ncell)
        kused = Vector{Int32}(undef, ncell)
        for p in 1:ncell
            alb[p] = cell_albedo[p, 1]
            n = cell_counts[p, 1]
            nval[p] = n
            kused[p] = _pooled_k_used(n, percentile)
        end
        alb_key, nval_key, kused_key = _pooled_layer_names(l)
        push!(out, alb_key => DimArray(alb, (point_dim,);
            metadata = Dict("units" => "1",
                            "long_name" => "bare-ice albedo, mean of the darkest " *
                                           "$(round(100 * percentile; digits = 2)) % of all " *
                                           "valid retrievals in the record ($(alb_key))",
                            "percentile" => percentile,
                            "albedo_range" => collect(albedo_range),
                            "n_dates" => length(date_list),
                            "date_range" => [string(first(date_list)), string(last(date_list))],
                            "note" => "Pooled over every date, NOT a per-year value and NOT " *
                                      "a reduction of per-year values. Not clamped to 1.")))
        push!(out, nval_key => DimArray(nval, (point_dim,);
            metadata = Dict("units" => "count",
                            "long_name" => "retrievals passing QA and albedo_range, " *
                                           "over the whole record")))
        push!(out, kused_key => DimArray(kused, (point_dim,);
            metadata = Dict("units" => "count",
                            "long_name" => "darkest retrievals averaged into the reported " *
                                           "albedo, ceil(percentile * n_valid)")))
    end

    centres = [_modis_cell_center(c...) for c in cell_list]
    push!(out, :latitude => DimArray([c[1] for c in centres], (point_dim,);
        metadata = Dict("units" => "degrees_north",
                        "long_name" => "latitude of the MODIS cell centre")))
    push!(out, :longitude => DimArray([c[2] for c in centres], (point_dim,);
        metadata = Dict("units" => "degrees_east",
                        "long_name" => "longitude of the MODIS cell centre")))
    push!(out, :cell_id => DimArray(_modis_cell_id.(cell_list), (point_dim,);
        metadata = Dict("long_name" => "MCD43A3 sinusoidal grid cell")))

    return DimStack(NamedTuple(out); metadata = Dict{String,Any}(
        "source" => "MODIS MCD43A3 v061, bare-ice albedo pooled over the whole record",
        "percentile" => percentile,
        "albedo_range" => collect(albedo_range),
        "qa_keep" => collect(qa_keep),
        "n_cells" => ncell,
        "n_dates" => length(date_list),
        "date_range" => [string(first(date_list)), string(last(date_list))],
        "cells_key" => cells_key,
        "note" => "One value per cell over all dates. Judge a cell by n_valid, not by the " *
                  "albedo alone; no sample-count threshold is applied.",
    ))
end
