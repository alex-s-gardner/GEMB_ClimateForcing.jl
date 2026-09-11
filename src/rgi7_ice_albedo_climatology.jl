"""
Multi-year reduction of the per-year glacier bare-ice albedo files written by
`data/run_rgi7_ice_albedo_modis.jl`.

Those files are the expensive part of the pipeline — ~1 TB of MCD43A3 download per year — and
each row is already one cell's **annual darkest-`percentile` mean**. So a climatology is a pure
read: no download, no re-fold, no GeoParquet dependency.

!!! note "This reduces annual statistics; it does not pool observations"
    `reduction` is applied to the per-year values, so `median` gives *the median of the annual
    bare-ice albedos*, weighting every year equally. That is usually what "bare-ice albedo
    climatology" means, and it is **not** the same quantity as one darkest-5 % taken over all
    years pooled — which would be dominated by whichever years happened to be darkest.

    Pooling across years is possible but is a different computation: it needs the per-date raw
    samples, which `compute_glacier_ice_albedo_modis` also leaves on disk under
    `<cache_path>/samples/`, and re-folding them is CPU minutes rather than a re-download.
    Nothing here reads those; the annual scalars are all this file needs.
"""

# Column order written by data/run_rgi7_ice_albedo_modis.jl's `write_intermediate`. Named here
# rather than assumed positionally at the parse site, and checked against the file's own header
# so a format change fails loudly instead of transposing albedo into a count.
const _RGI7_ALBEDO_COLUMNS = ("h", "v", "row", "col",
                              "albedo_bsa", "albedo_wsa", "n_valid_bsa", "n_valid_wsa")

# The two albedo layers, paired with the count column that qualifies each.
const _RGI7_ALBEDO_LAYERS = (:albedo_bsa => :n_valid_bsa, :albedo_wsa => :n_valid_wsa)

"""
    rgi7_ice_albedo_path(dir, year, hemisphere) -> String

Path of one per-year intermediate, matching what `data/run_rgi7_ice_albedo_modis.jl` writes.
"""
rgi7_ice_albedo_path(dir::AbstractString, year::Integer, hemisphere::Symbol) =
    joinpath(dir, "rgi7_ice_albedo_modis_$(year)_$(hemisphere).csv.gz")

"""
    _rgi7_read_albedo_year!(bsa, wsa, uniq, path) -> Int

Scatter one intermediate into per-unique-cell columns, returning the number of rows placed.

Rows are located by `searchsortedfirst` on the canonical cell list rather than by trusting the
file's order. The two *are* built from the same sort, but a north file and a south file
**interleave** in it (`h` is the outer key, so `h17v02` sorts before `h17v10`), so a positional
assumption would mis-attribute rows silently rather than erroring.
"""
function _rgi7_read_albedo_year!(bsa::AbstractVector{Float32}, wsa::AbstractVector{Float32},
                                 uniq::Vector{NTuple{4,Int}}, path::AbstractString)
    n = 0
    open(GzipDecompressorStream, path) do io
        header_seen = false
        for line in eachline(io)
            (isempty(line) || startswith(line, '#')) && continue
            if !header_seen && startswith(line, "h,")
                cols = Tuple(split(line, ','))
                cols == _RGI7_ALBEDO_COLUMNS || error(
                    "$(path) has columns $(cols); expected $(_RGI7_ALBEDO_COLUMNS). The file " *
                    "was written by a different version of data/run_rgi7_ice_albedo_modis.jl.")
                header_seen = true
                continue
            end
            f = split(line, ',')
            length(f) == 8 || error("malformed row in $(path): $line")
            cell = (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]), parse(Int, f[4]))
            k = searchsortedfirst(uniq, cell)
            (k <= length(uniq) && uniq[k] == cell) || error(
                "$(path) has cell $(cell), which is absent from the vendored cell list; the " *
                "run and data/rgi7_modis_cells.csv.gz are out of sync")
            bsa[k] = parse(Float32, f[5])
            wsa[k] = parse(Float32, f[6])
            n += 1
        end
        header_seen || error("$(path) has no header row")
    end
    return n
end

"""
    rgi7_ice_albedo_climatology(years; reduction=median, dir, hemispheres, min_years,
                                cells, verbose) -> DimStack

Reduce the per-year glacier bare-ice albedo files over `years`, per MODIS cell.

`reduction` is a **function handle** applied to each cell's vector of valid annual values —
`median` (default), `mean`, `minimum`, `maximum`, `std`, or anything else that maps a
`Vector{Float32}` to a number, e.g. `x -> quantile(x, 0.25)`. It never sees a `NaN`: unresolved
cell-years are dropped first, which is also what makes `:n_years` meaningful. A cell with fewer
than `min_years` valid years is reported as `NaN` rather than reduced over almost nothing.

Returns a `DimStack` over `Dim{:point}`, aligned with `rgi7_modis_unique_cells(cells)` — the
distinct cells in canonical order — carrying the reduced `:albedo_bsa` and `:albedo_wsa`, the
per-cell `:n_years_bsa` / `:n_years_wsa`, and the cell centre `:latitude` / `:longitude`.

# Keywords
- `reduction = median`: applied per cell over that cell's valid annual values.
- `dir = joinpath(@__DIR__, "..", "data")`: where the per-year files live.
- `hemispheres = (:north, :south)`: which halves to read. Both, unless only one was run.
- `min_years = 1`: minimum valid years for a cell to be reduced at all. Raise it (e.g. `10`)
  when the reduction is only meaningful over a decent sample.
- `cells = rgi7_modis_cells()`: pass a pre-loaded table to avoid re-reading it.
- `verbose = true`: report per-year coverage and the resolved fraction.

# Example
```julia
clim = rgi7_ice_albedo_climatology(2001:2024)                       # median annual albedo
clim = rgi7_ice_albedo_climatology(2001:2024; reduction = mean, min_years = 10)
clim = rgi7_ice_albedo_climatology(2001:2024; reduction = x -> quantile(x, 0.1))

clim[:albedo_bsa]        # reduced black-sky albedo per cell
clim[:n_years_bsa]       # how many years each cell actually contributed
```

!!! warning "Memory scales with cells × years"
    The annual values are materialised before reduction, because an arbitrary function handle
    cannot be folded incrementally. At the full 3.36M cells that is ~13 MB per layer-year, so
    ~670 MB for two layers over 25 years. Pass a shorter `years`, or one hemisphere at a time,
    if that is too much.
"""
function rgi7_ice_albedo_climatology(years;
                                     reduction = median,
                                     dir::AbstractString = joinpath(@__DIR__, "..", "data"),
                                     hemispheres = (:north, :south),
                                     min_years::Integer = 1,
                                     cells::Union{Nothing,RGI7ModisCells} = nothing,
                                     verbose::Bool = true)
    year_list = vec(collect(years))
    isempty(year_list) && throw(ArgumentError("no years requested"))
    min_years >= 1 || throw(ArgumentError("min_years must be at least 1, got $(min_years)"))
    isempty(hemispheres) && throw(ArgumentError("no hemispheres requested"))

    t = isnothing(cells) ? rgi7_modis_cells() : cells
    uniq = rgi7_modis_unique_cells(t)
    nu = length(uniq)
    ny = length(year_list)

    # nu x ny per layer. NaN is "no value", which is also what the files themselves carry for
    # a cell that failed QC or fell below min_samples, so no separate mask is needed.
    bsa = fill(NaN32, nu, ny)
    wsa = fill(NaN32, nu, ny)

    for (yi, y) in enumerate(year_list)
        got = 0
        for hemi in hemispheres
            path = rgi7_ice_albedo_path(dir, y, hemi)
            isfile(path) || error("""
                No bare-ice albedo file for $(y) $(hemi) at
                  $(path)
                Produce it with:
                  julia --project=. -t auto data/run_rgi7_ice_albedo_modis.jl $(y)
                or pass `hemispheres = (:$(first(hemispheres)),)` if only one half was run.
                """)
            got += _rgi7_read_albedo_year!(view(bsa, :, yi), view(wsa, :, yi), uniq, path)
        end
        if verbose
            res = count(!isnan, view(bsa, :, yi))
            @info "read" year = y rows = got coverage_pct = round(100 * got / nu; digits = 1) resolved_pct =
                round(100 * res / nu; digits = 1)
        end
        # Only a complaint when *both* halves were asked for: with `hemispheres` restricted the
        # other hemisphere's cells are supposed to stay NaN, and warning about it every year
        # would train the user to ignore the one message that means something.
        both = :north in hemispheres && :south in hemispheres
        (got == nu || !both) || @warn """
            year $(y) covers $(got) of $(nu) cells although both hemispheres were requested;
            the rest stay NaN for this year, so that run is incomplete.
            """ year = y rows = got expected = nu
    end

    out = Pair{Symbol,Any}[]
    point_dim = Dim{:point}(1:nu)
    buf = Vector{Float32}(undef, ny)
    for (layer, values) in ((:albedo_bsa, bsa), (:albedo_wsa, wsa))
        red = fill(NaN32, nu)
        nyr = zeros(Int32, nu)
        for i in 1:nu
            m = 0
            @inbounds for j in 1:ny
                v = values[i, j]
                isnan(v) && continue
                m += 1
                buf[m] = v
            end
            nyr[i] = m
            # `reduction` never sees a NaN, so `median`/`mean` need no nan-aware variant, and
            # a cell with too few years is reported unresolved rather than quietly reduced
            # over one value.
            m >= min_years && (red[i] = Float32(reduction(view(buf, 1:m))))
        end
        push!(out, layer => DimArray(red, (point_dim,);
            metadata = Dict("units" => "1",
                            "long_name" => "$(_rgi7_reduction_name(reduction)) of annual " *
                                           "glacier bare-ice albedo ($(String(layer)))",
                            "reduction" => _rgi7_reduction_name(reduction),
                            "years" => collect(year_list),
                            "min_years" => min_years)))
        push!(out, Symbol("n_years_", split(String(layer), '_')[2]) =>
            DimArray(nyr, (point_dim,);
                metadata = Dict("units" => "count",
                                "long_name" => "years with a resolved annual value")))
    end

    lat = Vector{Float64}(undef, nu)
    lon = Vector{Float64}(undef, nu)
    for (i, c) in enumerate(uniq)
        lat[i], lon[i] = _modis_cell_center(c...)
    end
    push!(out, :latitude => DimArray(lat, (point_dim,);
        metadata = Dict("units" => "degrees_north",
                        "long_name" => "latitude of the MODIS cell centre")))
    push!(out, :longitude => DimArray(lon, (point_dim,);
        metadata = Dict("units" => "degrees_east",
                        "long_name" => "longitude of the MODIS cell centre")))

    if verbose
        res = count(!isnan, out[1].second)
        @info "climatology" cells = nu years = ny reduction = _rgi7_reduction_name(reduction) resolved =
            res resolved_pct = round(100 * res / nu; digits = 1)
    end

    return DimStack(NamedTuple(out); metadata = Dict{String,Any}(
        "source" => "MODIS MCD43A3 v061 annual glacier bare-ice albedo, reduced over years",
        "reduction" => _rgi7_reduction_name(reduction),
        "years" => collect(year_list),
        "hemispheres" => String.(collect(hemispheres)),
        "min_years" => min_years,
        "n_cells" => nu,
        "note" => "Reduction of ANNUAL darkest-percentile means, weighting years equally — " *
                  "not a percentile pooled across years.",
    ))
end

# A printable name for the reduction, for metadata. Anonymous functions have gensym'd names
# like "#17#18", which are meaningless in a file, so those degrade to "custom".
function _rgi7_reduction_name(reduction)
    s = string(nameof(reduction))
    return (isempty(s) || startswith(s, '#')) ? "custom" : s
end
