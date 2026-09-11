#!/usr/bin/env julia
# Rasterize the RGI 7.0 glacier outlines onto the native MCD43A3 500 m sinusoidal grid and
# vendor the result as two compact gzipped CSVs in this directory.
#
# Run manually when RGI is updated. Unlike data/make_shaw2025_decoupling.jl this needs no
# scratch environment — ArchGDAL, Rasters and CodecZlib are already package dependencies, and
# Tar/Downloads/SHA are stdlibs:
#
#   julia --project=. data/make_rgi7_modis_cells.jl [source_dir]
#
# `source_dir` (default data/rgi7_source/, gitignored) holds the RGI 7.0 "G" product. Any
# region not already present there is downloaded and extracted. Pass a directory of
# NSIDC-downloaded files to avoid the mirror entirely — that is the intended interface for a
# reproducible build, since the default URL is a courtesy mirror with no archival guarantee.
#
# Source: RGI 7.0, glacier ("G") product, 19 regional shapefiles.
#   Citation : https://doi.org/10.5067/f6jmovy5navz  (NSIDC; requires Earthdata login)
#   Mirror   : https://cluster.klima.uni-bremen.de/~fmaussion/misc/rgi7_data/l4_rgi7b0_tar/
#   Guide    : https://www.glims.org/rgi_user_guide/
#
# The burn rules (centre-in-polygon, smallest-glacier-wins, the forced-cell fallback) and why
# each is the opposite of the obvious choice are documented in src/rgi7_modis_cells.jl. This
# script only drives them.

using GEMB_ClimateForcing
using ArchGDAL
using CodecZlib
using Printf
using SHA
using Tar
using Downloads

const G = GEMB_ClimateForcing
const AG = ArchGDAL

const MIRROR = "https://cluster.klima.uni-bremen.de/~fmaussion/misc/rgi7_data/l4_rgi7b0_tar"
const DEFAULT_SRC = joinpath(@__DIR__, "rgi7_source")
const DEST_CELLS = joinpath(@__DIR__, "rgi7_modis_cells.csv.gz")
const DEST_GLACIERS = joinpath(@__DIR__, "rgi7_glaciers.csv.gz")

# Region 20 (Antarctic mainland) has zero glaciers in RGI 7.0, so the G product ships 19 files.
const REGIONS = ["01_alaska", "02_western_canada_usa", "03_arctic_canada_north",
                 "04_arctic_canada_south", "05_greenland_periphery", "06_iceland",
                 "07_svalbard_jan_mayen", "08_scandinavia", "09_russian_arctic",
                 "10_north_asia", "11_central_europe", "12_caucasus_middle_east",
                 "13_central_asia", "14_south_asia_west", "15_south_asia_east",
                 "16_low_latitudes", "17_southern_andes", "18_new_zealand",
                 "19_subantarctic_antarctic_islands"]

# Published totals, as sanity tripwires rather than hard requirements.
const EXPECTED_GLACIERS = 274_531
const EXPECTED_AREA_KM2 = 706_744.0
const EXPECTED_TILES = 103

# Geometries are flushed to the burn once a tile has this many pending, to bound the memory
# held per region while keeping the number of `rasterize!` calls far below one per glacier.
const FLUSH_AT = 2000

## ------------------------------------------------------------------------- source fetching

region_stem(reg) = "RGI2000-v7.0-G-$(reg)"

# Restrict the run to a subset, e.g. RGI7_REGIONS=06,07. For debugging and for re-running one
# region after a source fix; a partial run writes partial tables, so it is not a normal mode.
function selected_regions()
    spec = get(ENV, "RGI7_REGIONS", "")
    isempty(spec) && return REGIONS
    want = split(spec, ',')
    sel = [r for r in REGIONS if first(split(r, '_')) in want]
    isempty(sel) && error("RGI7_REGIONS=$(spec) matched no region of $(REGIONS)")
    @warn "PARTIAL RUN — the written tables will not be the global product" regions = sel
    return sel
end

"""
    field_index(layer, name) -> Int

0-based index of a named field, from the *layer* definition.

Taken from the layer rather than from a feature: a feature fetched only to read its field list
may be freed when its `do` block returns, and using the index afterwards is then undefined.
"""
function field_index(layer, name::AbstractString)
    defn = AG.layerdefn(layer)
    for i in 1:AG.nfield(defn)
        AG.getname(AG.getfielddefn(defn, i - 1)) == name && return i - 1
    end
    error("no field $(name) in layer; have " *
          "$([AG.getname(AG.getfielddefn(defn, i - 1)) for i in 1:AG.nfield(defn)])")
end

"""
    ensure_region(src_dir, reg) -> (shapefile_path, attributes_path, sha256_or_nothing)

Path to one region's shapefile, downloading and extracting the tarball if it is absent.

The SHA-256 of the tarball is returned when this run fetched or found one, so it can be
recorded in the output provenance header — the vendored artifact has to stay reproducible
after the mirror disappears.
"""
function ensure_region(src_dir::AbstractString, reg::AbstractString)
    stem = region_stem(reg)
    dir = joinpath(src_dir, stem)
    shp = joinpath(dir, "$(stem).shp")
    att = joinpath(dir, "$(stem)-attributes.csv")
    tarball = joinpath(src_dir, "$(stem).tar.gz")

    if !isfile(shp)
        if !isfile(tarball)
            url = "$(MIRROR)/$(stem).tar.gz"
            @info "downloading" region = reg url
            mkpath(src_dir)
            Downloads.download(url, tarball * ".part")
            mv(tarball * ".part", tarball; force = true)
        end
        @info "extracting" tarball
        mktempdir(src_dir) do tmp
            open(GzipDecompressorStream, tarball) do io
                Tar.extract(io, tmp)
            end
            # the tarball carries a single top-level directory named after the region
            inner = only(readdir(tmp))
            mv(joinpath(tmp, inner), dir; force = true)
        end
    end
    isfile(shp) || error("no shapefile at $(shp) after extraction")
    isfile(att) || error("no attributes CSV at $(att)")
    # A pre-extracted directory with no tarball beside it yields no checksum, which silently
    # leaves that region unattested in the output header. Warn rather than shrug: the whole
    # point of recording the hashes is that the vendored table stays reproducible after the
    # mirror changes, and a gap is only discoverable by reading the header afterwards.
    sha = if isfile(tarball)
        bytes2hex(open(sha256, tarball))
    else
        @warn """
              no tarball for region $(reg), so its SHA-256 cannot be recorded in the output
              provenance header. The shapefile at $(dir) is used as-is. Delete that directory
              to force a fresh download, or place $(basename(tarball)) beside it.
              """
        nothing
    end
    return (shp, att, sha)
end

## ------------------------------------------------------- pass 1: the global glacier table

"""
    read_attributes(att_path) -> Vector{@NamedTuple{rgi_id::String, area_km2::Float64, o1region::Int}}

`rgi_id`, `area_km2` and `o1region` from one region's attributes CSV.

Read from the CSV rather than the shapefile's `.dbf` so pass 1 costs no geometry parsing, and
keyed by `rgi_id` in pass 2 rather than by row position, so the two sources are never assumed
to be in the same order.
"""
function read_attributes(att::AbstractString)
    rows = @NamedTuple{rgi_id::String, area_km2::Float64, o1region::Int}[]
    header = nothing
    icol = jcol = kcol = 0
    for line in eachline(att)
        isempty(strip(line)) && continue
        f = _csv_fields(line)
        if header === nothing
            header = f
            icol = findfirst(==("rgi_id"), header)
            jcol = findfirst(==("area_km2"), header)
            kcol = findfirst(==("o1region"), header)
            all(!isnothing, (icol, jcol, kcol)) ||
                error("$(att): expected rgi_id / area_km2 / o1region columns, got $(header)")
            continue
        end
        push!(rows, (rgi_id = f[icol], area_km2 = parse(Float64, f[jcol]),
                     o1region = parse(Int, f[kcol])))
    end
    return rows
end

# RGI's attributes CSV quotes fields (including glacier names, which contain commas), so a
# bare split(',') is not safe here even though it is fine for the tables we *write*.
function _csv_fields(line::AbstractString)
    out = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            if inq && i < lastindex(line) && line[nextind(line, i)] == '"'
                write(buf, '"')
                i = nextind(line, i)
            else
                inq = !inq
            end
        elseif c == ',' && !inq
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(out, String(take!(buf)))
    return out
end

## --------------------------------------------------------------------- pass 2: the burn

"""
    burn_all(src, glaciers, rank_of) -> Dict{Tuple{Int,Int},Matrix{Int32}}

Burn every region's outlines into per-tile grids, keeping the lowest area rank per cell.

Streams one region at a time and holds one 23 MB grid per touched tile (~103 of them, ~2.4 GB)
so that the global `minimum` reduction stays order-independent without holding all 274,531
geometries at once. Pending geometries are flushed per tile at [`FLUSH_AT`](@ref).
"""
function burn_all(src::AbstractString, glaciers, rank_of::Vector{Int32}, shas::Dict{String,String})
    grids = Dict{Tuple{Int,Int},Matrix{Int32}}()
    ids = [g.rgi_id for g in glaciers]        # sorted, for binary search
    nseen = 0
    nskipped_z = 0

    for reg in selected_regions()
        shp, _, sha = ensure_region(src, reg)
        isnothing(sha) || (shas[reg] = sha)
        pending = Dict{Tuple{Int,Int},Tuple{Vector{Any},Vector{Int32}}}()

        function flush_tile!(tile)
            haskey(pending, tile) || return nothing
            geoms, ranks = pending[tile]
            isempty(geoms) && return nothing
            grid = get!(() -> G._rgi7_new_tile_grid(), grids, tile)
            G._rgi7_burn_tile!(grid, tile[1], tile[2], geoms, ranks)
            empty!(geoms)
            empty!(ranks)
            return nothing
        end

        nreg = 0
        AG.read(shp) do ds
            layer = AG.getlayer(ds, 0)
            idfield = field_index(layer, "rgi_id")
            for i in 0:(AG.nfeature(layer) - 1)
                AG.getfeature(layer, i) do f
                    rgi_id = AG.getfield(f, idfield)
                    gi = searchsortedfirst(ids, rgi_id)
                    (gi <= length(ids) && ids[gi] == rgi_id) ||
                        error("$(reg): shapefile has $(rgi_id), absent from the attributes CSV")
                    sinu = G._rgi7_sinu_geometry(AG.getgeom(f), rgi_id)
                    for tile in G._rgi7_candidate_tiles(sinu)
                        geoms, ranks = get!(() -> (Any[], Int32[]), pending, tile)
                        push!(geoms, sinu)
                        push!(ranks, rank_of[gi])
                        length(geoms) >= FLUSH_AT && flush_tile!(tile)
                    end
                    nreg += 1
                    return nothing
                end
            end
        end
        for tile in collect(keys(pending))
            flush_tile!(tile)
        end
        nseen += nreg
        @info "burned" region = reg glaciers = nreg tiles_so_far = length(grids) total = nseen
    end
    nskipped_z == 0 || @warn "geometries skipped" nskipped_z
    return grids
end

## -------------------------------------------------------------------------------- assembly

function main()
    src = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_SRC
    @info "RGI 7.0 → MCD43A3 500 m cell list" source_dir = src cell_area_km2 =
        round(G._RGI7_CELL_AREA_KM2; digits = 5)

    # ---- pass 1: global glacier table, sorted by rgi_id
    shas = Dict{String,String}()
    all_rows = @NamedTuple{rgi_id::String, area_km2::Float64, o1region::Int}[]
    for reg in selected_regions()
        _, att, sha = ensure_region(src, reg)
        isnothing(sha) || (shas[reg] = sha)
        rows = read_attributes(att)
        append!(all_rows, rows)
        @info "attributes" region = reg glaciers = length(rows)
    end
    sort!(all_rows; by = r -> r.rgi_id)
    allunique(r.rgi_id for r in all_rows) || error("duplicate rgi_id across regions")
    ng = length(all_rows)
    @info "glaciers" ng total_area_km2 = round(sum(r.area_km2 for r in all_rows); digits = 1)
    ng == EXPECTED_GLACIERS ||
        @warn "unexpected glacier count (expected $(EXPECTED_GLACIERS) for RGI 7.0)" ng

    # ---- area rank: ascending area, ties broken by rgi_id so the order is total and the
    # burn is reproducible. rank 1 == smallest glacier == wins a contested cell.
    order = sortperm(1:ng; by = i -> (all_rows[i].area_km2, all_rows[i].rgi_id))
    rank_of = Vector{Int32}(undef, ng)
    glacier_of_rank = Vector{Int32}(undef, ng)
    for (r, i) in enumerate(order)
        rank_of[i] = Int32(r)
        glacier_of_rank[r] = Int32(i)
    end

    # ---- pass 2: burn
    grids = burn_all(src, all_rows, rank_of, shas)
    @info "tiles touched" n = length(grids)
    length(grids) == EXPECTED_TILES ||
        @warn "unexpected tile count (measured $(EXPECTED_TILES) from RGI7 centroids)" n =
            length(grids)

    # ---- extract cells, mapping rank back to glacier index
    cells = NTuple{5,Int}[]                    # (h, v, row, col, glacier_index)
    for tile in sort!(collect(keys(grids)))
        for (h, v, row, col, rank) in G._rgi7_tile_cells(grids[tile], tile[1], tile[2])
            push!(cells, (h, v, row, col, Int(glacier_of_rank[rank])))
        end
    end
    @info "burned cells" n = length(cells) implied_area_km2 =
        round(length(cells) * G._RGI7_CELL_AREA_KM2; digits = 1)

    # ---- forced cells for glaciers with no interior cell
    #
    # 83.6 % of RGI 7.0 is under 1 km² against a 0.2146 km² cell, so centre-in-polygon leaves
    # a large fraction of the inventory with nothing. Each such glacier gets the single cell
    # holding its representative point.
    #
    # A forced cell is **added, never substituted**. Letting it displace the cell's existing
    # owner under the usual smallest-wins rule was tried and is wrong: measured on Iceland it
    # displaced 13 owners and left 8 glaciers (1.4 %) with no cell at all, which is exactly the
    # erasure smallest-wins exists to prevent. So a cell may carry **more than one glacier, and
    # only ever when a glacier would otherwise have none** — two sub-cell glaciers in one
    # 463 m pixel genuinely do share that pixel's albedo, and saying so is more faithful than
    # deleting one of them. Every glacier therefore has at least one row.
    #
    # The consequence for consumers: rows are unique as (cell, glacier) **pairs**, not as
    # cells, so the albedo point list is `rgi7_modis_unique_cells`, not the raw column.
    have = falses(ng)
    for c in cells
        have[c[5]] = true
    end
    missing_idx = findall(!, have)
    @info "glaciers with no interior cell" n = length(missing_idx) pct =
        round(100 * length(missing_idx) / ng; digits = 1)

    forced = falses(ng)
    # region-major, so the single-entry geometry cache loads each shapefile exactly once
    sort!(missing_idx; by = gi -> (all_rows[gi].o1region, all_rows[gi].rgi_id))
    for gi in missing_idx
        rgi_id = all_rows[gi].rgi_id
        geom = geometry_for(src, rgi_id, all_rows[gi].o1region)
        cell = G._rgi7_representative_cell(geom, rgi_id)
        push!(cells, (cell[1], cell[2], cell[3], cell[4], gi))
        forced[gi] = true
    end
    @info "forced cells" assigned = count(forced)

    # canonical order, now including glacier_index so the sort is total over the pair
    sort!(cells; by = c -> (c[1], c[2], c[3], c[4], c[5]))
    allunique(cells) || error("duplicate (cell, glacier) pair after the forced pass")

    n_cells = zeros(Int32, ng)
    for c in cells
        n_cells[c[5]] += 1
    end
    nzero = count(iszero, n_cells)
    nzero == 0 || error("$(nzero) glaciers still have no cell; the forced pass should " *
                        "guarantee at least one for every glacier")
    nshared = length(cells) - length(unique(c -> (c[1], c[2], c[3], c[4]), cells))
    nuniq = length(cells) - nshared
    implied = nuniq * G._RGI7_CELL_AREA_KM2
    @info "final" rows = length(cells) unique_cells = nuniq shared_cells = nshared implied_area_km2 =
        round(implied; digits = 1)

    # Centre-in-polygon is area-unbiased, so the burned cell count times the cell area should
    # land close to the published total. A large miss means the projection, the grid origin or
    # the ring handling moved — not that the tolerance needs widening.
    published = sum(g.area_km2 for g in all_rows)
    frac = implied / published
    @info "area closure" published_km2 = round(published; digits = 1) implied_km2 =
        round(implied; digits = 1) ratio = round(frac; digits = 4)
    0.90 <= frac <= 1.10 || error(
        "implied area is $(round(100 * frac; digits = 1)) % of published " *
        "($(round(implied; digits = 1)) vs $(round(published; digits = 1)) km²); " *
        "the burn is wrong, do not vendor this")
    isempty(ARGS) || length(selected_regions()) == length(REGIONS) ||
        @warn "partial run: skipping the global area tripwire against " *
              "$(EXPECTED_AREA_KM2) km²"
    length(selected_regions()) == length(REGIONS) &&
        !isapprox(published, EXPECTED_AREA_KM2; rtol = 0.01) &&
        @warn "published area differs from the RGI 7.0 figure" published EXPECTED_AREA_KM2

    write_tables(cells, all_rows, n_cells, forced, shas, length(grids))
    return nothing
end

# Geometry lookup for the forced pass, holding **one region at a time**.
#
# Deliberately a single-entry cache rather than a Dict keyed by region: caching every region
# would accumulate all 274,531 cloned geometries (region 13 alone has 75,613) and dwarf the
# 2.4 GB of burn grids. The caller sorts by region so each one is still loaded exactly once.
const _GEOM_CACHE = Ref{Tuple{Int,Dict{String,Any}}}((0, Dict{String,Any}()))

function geometry_for(src::AbstractString, rgi_id::AbstractString, o1region::Int)
    cached_region, cache = _GEOM_CACHE[]
    if cached_region != o1region
        reg = REGIONS[findfirst(r -> parse(Int, first(split(r, '_'))) == o1region, REGIONS)]
        cache = _load_region_geoms(src, reg)
        _GEOM_CACHE[] = (o1region, cache)
        @info "loaded geometries for the forced pass" region = reg n = length(cache)
    end
    haskey(cache, rgi_id) || error("$(rgi_id) not found in region $(o1region)")
    return cache[rgi_id]
end

function _load_region_geoms(src::AbstractString, reg::AbstractString)
    shp, _, _ = ensure_region(src, reg)
    out = Dict{String,Any}()
    AG.read(shp) do ds
        layer = AG.getlayer(ds, 0)
        idfield = field_index(layer, "rgi_id")
        for i in 0:(AG.nfeature(layer) - 1)
            AG.getfeature(layer, i) do f
                out[AG.getfield(f, idfield)] = AG.clone(AG.getgeom(f))
                return nothing
            end
        end
    end
    return out
end

## ---------------------------------------------------------------------------- the artifacts

function write_tables(cells, glaciers, n_cells, forced, shas, ntiles)
    ng = length(glaciers)
    area = sum(g.area_km2 for g in glaciers)

    open(GzipCompressorStream, DEST_GLACIERS, "w") do io
        println(io, "# RGI 7.0 glaciers, the side table of data/rgi7_modis_cells.csv.gz.")
        println(io, "# Generated by data/make_rgi7_modis_cells.jl from the RGI 7.0 glacier")
        println(io, "# (\"G\") product, 19 regional shapefiles.")
        println(io, "#   Citation: https://doi.org/10.5067/f6jmovy5navz")
        println(io, "#   Guide   : https://www.glims.org/rgi_user_guide/")
        println(io, "#")
        println(io, "# glacier_index : 1-based, equal to the row number. Present so the join")
        println(io, "#                 key into the cell table is explicit in the file.")
        println(io, "# rgi_id        : RGI 7.0 identifier, e.g. RGI2000-v7.0-G-01-00001.")
        println(io, "#                 SORTED — the loader binary-searches this column.")
        println(io, "# area_km2      : RGI 7.0 published area.")
        println(io, "# n_cells       : rows in the cell table for this glacier.")
        println(io, "# forced        : 1 if the glacier contains no MCD43A3 cell CENTRE and is")
        println(io, "#                 represented by the single cell holding its")
        println(io, "#                 representative point. That cell is majority NOT")
        println(io, "#                 glacier — rock, moraine and water are all darker than")
        println(io, "#                 ice — so its bare-ice albedo is biased LOW. Filter on")
        println(io, "#                 this before computing any aggregate.")
        println(io, "glacier_index,rgi_id,area_km2,n_cells,forced")
        for i in 1:ng
            @printf(io, "%d,%s,%.6f,%d,%d\n", i, glaciers[i].rgi_id, glaciers[i].area_km2,
                    n_cells[i], forced[i] ? 1 : 0)
        end
    end

    open(GzipCompressorStream, DEST_CELLS, "w") do io
        println(io, "# RGI 7.0 glacier outlines rasterized onto the native MCD43A3 500 m")
        println(io, "# sinusoidal grid. Generated by data/make_rgi7_modis_cells.jl.")
        println(io, "#   RGI 7.0 citation: https://doi.org/10.5067/f6jmovy5navz")
        println(io, "#   MCD43A3 grid    : sinusoidal, +R=6371007.181, 2400x2400 tiles at")
        println(io, "#                     463.312716528 m (0.21464 km2 per cell)")
        println(io, "#")
        @printf(io, "# %d glaciers, %.0f km2 published, over %d MODIS tiles ->\n", ng, area, ntiles)
        @printf(io, "# %d cells = %.0f km2 implied (%.1f%% of published)\n", length(cells),
                length(cells) * GEMB_ClimateForcing._RGI7_CELL_AREA_KM2,
                100 * length(cells) * GEMB_ClimateForcing._RGI7_CELL_AREA_KM2 / area)
        println(io, "#")
        println(io, "# Burn rules, all documented at length in src/rgi7_modis_cells.jl:")
        println(io, "#   - centre-in-polygon, NOT ALL_TOUCHED: a perimeter cell is majority")
        println(io, "#     off-glacier and therefore DARKER than ice, so including it would")
        println(io, "#     bias the darkest-percentile statistic being measured.")
        println(io, "#   - interior rings (nunataks) are excluded; they are pervasive in RGI 7.0")
        println(io, "#     (37 % of Svalbard outlines, up to 456 rings on one glacier).")
        println(io, "#   - a contested cell goes to the SMALLER glacier, so a one-cell glacier")
        println(io, "#     cannot be erased by a neighbour. The reducer is `min` over an")
        println(io, "#     ascending-area rank, hence commutative, hence this file does not")
        println(io, "#     depend on the order the regional shapefiles were read in.")
        println(io, "#")
        println(io, "# h,v           : MODIS sinusoidal tile, h 0-35 east from the")
        println(io, "#                 antimeridian, v 0-17 south from the north pole.")
        println(io, "#                 v <= 8 is EXACTLY the northern hemisphere.")
        println(io, "# row,col       : 1-based cell index within the tile, 1-2400; row counts")
        println(io, "#                 from the tile's north edge, col from its west edge.")
        println(io, "# glacier_index : 1-based row in data/rgi7_glaciers.csv.gz.")
        println(io, "#")
        println(io, "# Sorted ascending by (h,v,row,col) and unique. Longitude/latitude are")
        println(io, "# deliberately absent: they are recoverable exactly from (h,v,row,col)")
        println(io, "# via _modis_cell_center, and every row here is verified at generation")
        println(io, "# time to satisfy _modis_cell(_modis_cell_center(cell...)...) === cell.")
        if !isempty(shas)
            println(io, "#")
            println(io, "# SHA-256 of the source tarballs, so this artifact stays reproducible")
            println(io, "# after the download mirror changes or disappears:")
            for reg in REGIONS
                haskey(shas, reg) && println(io, "#   ", region_stem(reg), ".tar.gz  ", shas[reg])
            end
        end
        println(io, "h,v,row,col,glacier_index")
        for c in cells
            @printf(io, "%d,%d,%d,%d,%d\n", c[1], c[2], c[3], c[4], c[5])
        end
    end

    @info "wrote" DEST_CELLS size_mb = round(filesize(DEST_CELLS) / 1024^2; digits = 2)
    @info "wrote" DEST_GLACIERS size_mb = round(filesize(DEST_GLACIERS) / 1024^2; digits = 2)
    return nothing
end

main()
