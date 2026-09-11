#!/usr/bin/env julia
# Turn the per-(year, hemisphere) intermediates written by data/run_rgi7_ice_albedo_modis.jl
# into the published point-based **GeoParquet** product.
#
# Run manually, in a SCRATCH ENVIRONMENT, because this is the only step that needs
# GeoParquet.jl — which is deliberately kept out of Project.toml for exactly the reason HDF5.jl
# is (see data/make_shaw2025_decoupling.jl): the package itself never reads or writes Parquet,
# so a heavy dependency chain (GeoParquet -> DataFrames + Parquet2 + WellKnownGeometry) has no
# business in the load path of anyone who only wants climate forcing.
#
#   julia --project=/some/scratch/env -e 'using Pkg
#       Pkg.add(["GeoParquet", "DataFrames", "GeoFormatTypes"])
#       Pkg.develop(path="/path/to/GEMB_ClimateForcing.jl")'
#   julia --project=/some/scratch/env data/make_rgi7_ice_albedo_parquet.jl 2019 [2020 ...]
#
# GEMB_ClimateForcing is needed only for the vendored cell list and `_modis_cell_center`.
#
# Fallback if GeoParquet.jl ever breaks on this row count: this repo's GDAL_jll exports
# RegisterOGRParquet and ships libarrow/libparquet, so ArchGDAL could write the file with no
# new dependency at all — slower to write, but already in the package's graph.

using GEMB_ClimateForcing
using CodecZlib
using DataFrames
using GeoFormatTypes
using GeoParquet
using Printf

const G = GEMB_ClimateForcing

# Overridable so a scale test cannot leave a synthetic product sitting in data/ where it
# would be indistinguishable from the real one.
const OUTPUT_DIR = get(ENV, "RGI7_ALBEDO_DIR", joinpath(@__DIR__))

# Melt-season windows, echoed per row so the file is self-describing about which season
# produced each value. `doy_start > doy_end` marks the southern window, which wraps New Year.
const DOY_NORTH = G._MODIS_MELT_SEASON_NORTH
const DOY_SOUTH = G._MODIS_MELT_SEASON_SOUTH

## ------------------------------------------------------------------------------------ input

intermediate_path(year::Integer, hemi::Symbol) =
    joinpath(OUTPUT_DIR, @sprintf("rgi7_ice_albedo_modis_%d_%s.csv.gz", year, hemi))

"""
    read_intermediate!(bsa, wsa, nb, nw, uniq, path) -> Int

Scatter one intermediate's rows into per-unique-cell vectors, returning the row count.

Rows are placed by `searchsortedfirst` on the global unique-cell list rather than by assuming
the file's order matches it. The two are in fact built from the same canonical sort, but a
north file and a south file **interleave** in that order (h is the outer key, so h17v02 sorts
before h17v10), so concatenating them would be wrong and an order assumption would be a silent
mis-attribution rather than an error.
"""
function read_intermediate!(bsa, wsa, nb, nw, uniq::Vector{NTuple{4,Int}},
                            path::AbstractString)
    isfile(path) || error("missing intermediate $(path); run data/run_rgi7_ice_albedo_modis.jl first")
    n = 0
    open(GzipDecompressorStream, path) do io
        for line in eachline(io)
            (isempty(line) || startswith(line, '#') || startswith(line, "h,")) && continue
            f = split(line, ',')
            length(f) == 8 || error("malformed row in $(path): $line")
            cell = (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]), parse(Int, f[4]))
            k = searchsortedfirst(uniq, cell)
            (k <= length(uniq) && uniq[k] == cell) ||
                error("$(path) has cell $(cell), absent from the vendored cell list; the " *
                      "intermediate and data/rgi7_modis_cells.csv.gz are out of sync")
            bsa[k] = parse(Float32, f[5])
            wsa[k] = parse(Float32, f[6])
            nb[k] = parse(Int32, f[7])
            nw[k] = parse(Int32, f[8])
            n += 1
        end
    end
    return n
end

## ----------------------------------------------------------------------------------- output

"""
    wkb_point(lon, lat) -> Vector{UInt8}

A 21-byte little-endian WKB `POINT(lon lat)`: byte order `01`, geometry type `1`, then the two
`Float64` coordinates.

Hand-built rather than routed through a geometry library because that is all a point is, and
it keeps the writer's dependency surface to GeoParquet itself. Axis order is (x = lon,
y = lat), which is what OGC:CRS84 means — *not* EPSG:4326's lat-first order.
"""
function wkb_point(lon::Float64, lat::Float64)
    b = Vector{UInt8}(undef, 21)
    b[1] = 0x01                       # little-endian
    b[2:5] = reinterpret(UInt8, [UInt32(1)])
    b[6:13] = reinterpret(UInt8, [lon])
    b[14:21] = reinterpret(UInt8, [lat])
    return b
end

function build_table(years::Vector{Int})
    t = rgi7_modis_cells()
    uniq = rgi7_modis_unique_cells(t)
    nu = length(uniq)
    nrow = length(t)
    @info "vendored cell list" rows = nrow unique_cells = nu glaciers = n_glaciers(t)

    # One row per (cell, glacier) PAIR, not per cell. For 99.6 % of rows those are the same
    # thing; the exception is a pixel shared by two sub-cell glaciers, where both are kept so
    # that no glacier is missing from the product and `rgi_id` stays meaningful on every row.
    # A shared pixel's albedo is therefore repeated — it is one measurement of one pixel.
    lon = Vector{Float64}(undef, nrow)
    lat = Vector{Float64}(undef, nrow)
    cellpos = Vector{Int}(undef, nrow)
    for i in 1:nrow
        c = (Int(t.h[i]), Int(t.v[i]), Int(t.row[i]), Int(t.col[i]))
        la, lo = G._modis_cell_center(c...)
        lat[i] = la
        lon[i] = lo
        k = searchsortedfirst(uniq, c)
        cellpos[i] = k
    end

    df = DataFrame(
        geometry = [GeoFormatTypes.WellKnownBinary(GeoFormatTypes.Geom(),
                                                   wkb_point(lon[i], lat[i]))
                    for i in 1:nrow],
        rgi_id = [t.rgi_id[t.glacier[i]] for i in 1:nrow],
        area_km2 = [t.area_km2[t.glacier[i]] for i in 1:nrow],
        forced = [t.forced[t.glacier[i]] for i in 1:nrow],
        h = t.h, v = t.v, row = t.row, col = t.col,
        longitude = lon, latitude = lat,
        doy_start = [Int16(t.v[i] <= 8 ? DOY_NORTH[1] : DOY_SOUTH[1]) for i in 1:nrow],
        doy_end = [Int16(t.v[i] <= 8 ? DOY_NORTH[2] : DOY_SOUTH[2]) for i in 1:nrow],
    )

    for y in years
        bsa = fill(NaN32, nu)
        wsa = fill(NaN32, nu)
        nb = zeros(Int32, nu)
        nw = zeros(Int32, nu)
        got = 0
        for hemi in (:north, :south)
            got += read_intermediate!(bsa, wsa, nb, nw, uniq, intermediate_path(y, hemi))
        end
        got == nu || error("year $(y): intermediates covered $(got) of $(nu) unique cells")
        df[!, Symbol("albedo_bsa_$(y)")] = bsa[cellpos]
        df[!, Symbol("albedo_wsa_$(y)")] = wsa[cellpos]
        df[!, Symbol("n_valid_bsa_$(y)")] = nb[cellpos]
        df[!, Symbol("n_valid_wsa_$(y)")] = nw[cellpos]
        res = count(!isnan, bsa)
        @info "year folded in" year = y unique_cells = nu resolved = res pct =
            round(100 * res / nu; digits = 1)
    end
    return df
end

function main()
    years = isempty(ARGS) ? [2019] : parse.(Int, ARGS)
    df = build_table(years)

    dest = joinpath(OUTPUT_DIR, @sprintf("rgi7_glacier_ice_albedo_modis_%d_%d.parquet",
                                         minimum(years), maximum(years)))
    @info "writing GeoParquet" dest rows = nrow(df) cols = ncol(df)
    # crs = nothing is spec-legal and means OGC:CRS84 (lon/lat, WGS84), which is what the
    # cell centres are. zstd is GeoParquet.jl's default and the right choice here: the table is
    # sorted by (h, v, row, col), so h/v/row and rgi_id are highly repetitive.
    GeoParquet.write(dest, df, (:geometry,))
    @info "wrote" dest size_mb = round(filesize(dest) / 1024^2; digits = 1) bytes_per_row =
        round(filesize(dest) / nrow(df); digits = 1)

    # Read it back through the same library, so a malformed file fails here rather than in
    # whatever tool the user reaches for.
    back = GeoParquet.read(dest)
    nrow(back) == nrow(df) || error("round trip lost rows: $(nrow(back)) vs $(nrow(df))")
    @info "round trip OK" rows = nrow(back) cols = names(back)[1:min(end, 6)]
    println("\nValidate independently with, e.g.:")
    println("  python -c \"import geopandas; g=geopandas.read_parquet('$(dest)'); print(g.crs, g.shape); print(g.head())\"")
    return nothing
end

main()
