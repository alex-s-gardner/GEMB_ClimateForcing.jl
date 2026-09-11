"""
Tests for the RGI 7.0 → MCD43A3 grid rasterization (`src/rgi7_modis_cells.jl`) and the
multi-year reduction of its albedo products (`src/rgi7_ice_albedo_climatology.jl`).

Fully offline. The burn-kernel testsets need no vendored data at all — they use synthetic
geometry plus two real RGI 7.0 outlines carried in `test/fixtures/rgi7_iceland_outlines.tsv`
— so they run before `data/make_rgi7_modis_cells.jl` has ever been executed. That ordering is
deliberate: a global burn takes hours, and it must not be the debug loop for a sign error.

The one testset that *does* need the vendored tables skips itself with an informative message
when they are absent, so a fresh checkout still passes.
"""

using Test
import ArchGDAL
import CodecZlib: GzipCompressorStream
import DimensionalData: lookup, metadata
using Statistics: mean, median, quantile
const _G = GEMB_ClimateForcing
const _AG = ArchGDAL

const _FIXTURE_TSV = joinpath(@__DIR__, "fixtures", "rgi7_iceland_outlines.tsv")

# Real RGI 7.0 outlines, keyed by rgi_id. Z is already dropped in the fixture.
function _load_fixture_outlines()
    out = Dict{String,NamedTuple}()
    for line in eachline(_FIXTURE_TSV)
        (isempty(line) || startswith(line, '#') || startswith(line, "rgi_id\t")) && continue
        f = split(line, '\t')
        out[f[1]] = (; area_km2 = parse(Float64, f[2]), n_rings = parse(Int, f[3]),
                       cenlon = parse(Float64, f[4]), cenlat = parse(Float64, f[5]),
                       geom = _AG.fromWKT(String(f[6])))
    end
    return out
end

# The outer ring of a sinusoidal polygon, with every interior ring stripped. Used to isolate
# what the holes actually exclude.
function _fill_holes(sinu)
    outer = _AG.getgeom(sinu, 0)
    return _AG.createpolygon([[(_AG.getx(outer, j - 1), _AG.gety(outer, j - 1))
                               for j in 1:_AG.ngeom(outer)]])
end

_cellset(cells) = Set((c[3], c[4]) for c in cells)

@testset "RGI7 → MODIS cells" begin

    # ---------------------------------------------------------------- the burn grid itself
    #
    # This is the load-bearing testset. Rasters evaluates its inside test at
    # `first(lookup) + (i-1)*step`, so if those values are not bit-identical to
    # `_modis_cell_center`'s arithmetic then the burn grid and the cell lookup are two
    # different grids that merely almost agree — and cells would silently shift by one near
    # tile edges.
    @testset "burn grid is bit-identical to the cell grid" begin
        for (h, v) in ((17, 2), (0, 0), (35, 17), (18, 8), (18, 9))
            xd, yd = _G._rgi7_tile_dims(h, v)
            ulx, uly = _G._modis_tile_origin(h, v)
            px = _G._MODIS_PIXEL_M

            @test first(lookup(xd)) === ulx + 0.5 * px
            @test first(lookup(yd)) === uly - 0.5 * px
            @test step(lookup(xd)) === px
            @test step(lookup(yd)) === -px
            @test length(lookup(xd)) == _G._MODIS_TILE_PIXELS
            @test length(lookup(yd)) == _G._MODIS_TILE_PIXELS

            # every value, not just the first — a Regular span would let a drifting
            # accumulation pass a first/step check while being wrong in the middle
            @test all(c -> lookup(xd)[c] === ulx + (c - 0.5) * px, 1:_G._MODIS_TILE_PIXELS)
            @test all(r -> lookup(yd)[r] === uly - (r - 0.5) * px, 1:_G._MODIS_TILE_PIXELS)
        end
    end

    @testset "tile edges are shared exactly" begin
        # Extends the existing edge test to the tiles the burn actually constructs: a
        # one-ulp gap between neighbours would put a boundary glacier in neither tile.
        for h in 0:(_G._MODIS_H_MAX - 1), v in 0:_G._MODIS_V_MAX
            @test _G._modis_tile_bounds(h, v)[3] === _G._modis_tile_bounds(h + 1, v)[1]
        end
        for v in 0:(_G._MODIS_V_MAX - 1), h in 0:_G._MODIS_H_MAX
            @test _G._modis_tile_bounds(h, v)[2] === _G._modis_tile_bounds(h, v + 1)[4]
        end
    end

    @testset "v ≤ 8 is exactly the northern hemisphere" begin
        # _MODIS_UL_Y == 9 * _MODIS_TILE_SPAN_M exactly, so the v=9 tile edge is the
        # equator to the bit. This is what lets the driver split hemispheres with an
        # integer comparison instead of computing 3.3M cell-centre latitudes.
        @test _G._modis_tile_edge_y(9) === 0.0
        @test 9 * _G._MODIS_TILE_SPAN_M === _G._MODIS_UL_Y
        for v in (0, 4, 8), col in (1, 1200, 2400)
            @test _G._modis_cell_center(18, v, 2400, col)[1] > 0
        end
        for v in (9, 13, 17), col in (1, 1200, 2400)
            @test _G._modis_cell_center(18, v, 1, col)[1] < 0
        end
    end

    # ------------------------------------------------------------------- indexing and burn
    @testset "grid is indexed [col, row]" begin
        # X is the first dimension of _rgi7_tile_dims, so getting this backwards would
        # transpose every glacier — silently, and only detectably as a geographic error.
        g = _G._rgi7_new_tile_grid()
        @test size(g) == (_G._MODIS_TILE_PIXELS, _G._MODIS_TILE_PIXELS)
        @test all(==(_G._RGI7_NO_GLACIER), g)
        g[1500, 1000] = Int32(3)                      # [col, row]
        @test _G._rgi7_tile_cells(g, 17, 2) == [(17, 2, 1000, 1500, 3)]
    end

    @testset "synthetic rectangle burns exactly the enclosed centres" begin
        h, v = 17, 2
        ulx, uly = _G._modis_tile_origin(h, v)
        px = _G._MODIS_PIXEL_M
        r0, c0 = 1000, 1500
        # inset by 1 cm so the rectangle unambiguously contains 5x5 centres and no others
        xlo, xhi = ulx + (c0 - 1) * px + 0.01, ulx + (c0 + 4) * px - 0.01
        yhi, ylo = uly - (r0 - 1) * px - 0.01, uly - (r0 + 4) * px + 0.01
        rect = _AG.createpolygon([[(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi), (xlo, ylo)]])

        grid = _G._rgi7_new_tile_grid()
        _G._rgi7_burn_tile!(grid, h, v, [rect], Int32[7])
        cells = _G._rgi7_tile_cells(grid, h, v)

        @test length(cells) == 25
        @test _cellset(cells) == Set((r, c) for r in r0:(r0 + 4), c in c0:(c0 + 4))
        @test all(c -> c[5] == 7, cells)              # rank carried through
        @test all(c -> c[1] == h && c[2] == v, cells)
        # canonical (row, col) ascending order, which the vendored file relies on
        @test issorted(cells, by = c -> (c[3], c[4]))
    end

    @testset "a sub-cell rectangle burns nothing" begin
        # The forced-cell rule exists because of exactly this: centre-in-polygon gives a
        # glacier smaller than 0.21464 km² no cell at all.
        h, v = 17, 2
        ulx, uly = _G._modis_tile_origin(h, v)
        px = _G._MODIS_PIXEL_M
        # 0.4 cell across, placed to straddle a centre without containing one
        x0 = ulx + 1500 * px + 0.55 * px
        y0 = uly - 1000 * px - 0.55 * px
        tiny = _AG.createpolygon([[(x0, y0 - 0.4px), (x0 + 0.4px, y0 - 0.4px),
                                   (x0 + 0.4px, y0), (x0, y0), (x0, y0 - 0.4px)]])
        grid = _G._rgi7_new_tile_grid()
        _G._rgi7_burn_tile!(grid, h, v, [tiny], Int32[1])
        @test isempty(_G._rgi7_tile_cells(grid, h, v))
    end

    @testset "tie-break: smaller glacier wins, order-independently" begin
        # Two properties in one. Smallest-wins so a one-cell glacier cannot be erased by a
        # neighbour; order-independent so the *vendored file* is not a function of the order
        # the 19 regional shapefiles happen to be read in.
        h, v = 17, 2
        ulx, uly = _G._modis_tile_origin(h, v)
        px = _G._MODIS_PIXEL_M
        r0, c0 = 1000, 1500
        xlo = ulx + (c0 - 1) * px + 0.01
        yhi, ylo = uly - (r0 - 1) * px - 0.01, uly - (r0 + 4) * px + 0.01
        rect(xhi) = _AG.createpolygon([[(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi),
                                        (xlo, ylo)]])
        big   = rect(ulx + (c0 + 4) * px - 0.01)      # 5 columns
        small = rect(ulx + (c0 + 2) * px - 0.01)      # 3 columns, so cols c0..c0+1 overlap

        a = _G._rgi7_new_tile_grid()
        _G._rgi7_burn_tile!(a, h, v, [small, big], Int32[1, 2])   # rank 1 == smaller area
        b = _G._rgi7_new_tile_grid()
        _G._rgi7_burn_tile!(b, h, v, [big, small], Int32[2, 1])
        @test a == b                                              # commutative reducer

        shared = [c for c in _G._rgi7_tile_cells(a, h, v) if c[4] in c0:(c0 + 1)]
        @test !isempty(shared)
        @test all(c -> c[5] == 1, shared)                         # smaller glacier keeps them
    end

    # -------------------------------------------------------- interior rings, on real data
    @testset "interior rings are excluded (real RGI7 outline)" begin
        # Nunataks are pervasive in RGI 7.0, not exotic: 23.5 % of Arctic Canada North
        # outlines have interior rings and one has 456 of them. Filling them silently would
        # inflate cell counts in exactly the big-ice-cap regions that dominate the area
        # total, so this is pinned against a real outline rather than a synthetic square.
        outlines = _load_fixture_outlines()
        f = outlines["RGI2000-v7.0-G-06-00241"]
        @test f.n_rings == 14                                     # 13 nunataks

        sinu = _G._rgi7_sinu_geometry(f.geom, "RGI2000-v7.0-G-06-00241")
        tiles = _G._rgi7_candidate_tiles(sinu)
        @test tiles == [(17, 2)]

        h, v = only(tiles)
        holed = _G._rgi7_new_tile_grid()
        _G._rgi7_burn_tile!(holed, h, v, [sinu], Int32[1])
        filled = _G._rgi7_new_tile_grid()
        _G._rgi7_burn_tile!(filled, h, v, [_fill_holes(sinu)], Int32[1])

        ch = _cellset(_G._rgi7_tile_cells(holed, h, v))
        cf = _cellset(_G._rgi7_tile_cells(filled, h, v))
        @test issubset(ch, cf)                        # holes can only remove cells
        @test length(cf) > length(ch)                 # and they do remove some
        @test length(ch) == 147                       # measured; a change means the burn moved
        @test length(cf) == 154

        # 32.2293 km² / 0.21464 km² = 150.2 cells expected — centre-in-polygon is
        # area-unbiased, so this ties the burn to the published area independently.
        @test isapprox(length(ch) * _G._RGI7_CELL_AREA_KM2, f.area_km2; rtol = 0.05)
    end

    @testset "a sub-cell real glacier needs the forced cell" begin
        outlines = _load_fixture_outlines()
        id = "RGI2000-v7.0-G-06-00147"
        f = outlines[id]
        @test f.area_km2 < _G._RGI7_CELL_AREA_KM2     # 0.0443 km² vs a 0.21464 km² cell

        sinu = _G._rgi7_sinu_geometry(f.geom, id)
        burned = 0
        for (h, v) in _G._rgi7_candidate_tiles(sinu)
            grid = _G._rgi7_new_tile_grid()
            _G._rgi7_burn_tile!(grid, h, v, [sinu], Int32[1])
            burned += length(_G._rgi7_tile_cells(grid, h, v))
        end
        @test burned == 0                             # nothing, hence the rule

        cell = _G._rgi7_representative_cell(f.geom, id)
        @test cell isa NTuple{4,Int}
        @test 0 <= cell[1] <= _G._MODIS_H_MAX && 0 <= cell[2] <= _G._MODIS_V_MAX
        @test 1 <= cell[3] <= _G._MODIS_TILE_PIXELS && 1 <= cell[4] <= _G._MODIS_TILE_PIXELS
        # PointOnSurface is inside the polygon, so the forced cell must be one the
        # glacier's own envelope reaches
        @test (cell[1], cell[2]) in _G._rgi7_candidate_tiles(sinu)
        # and it must survive the same round trip every emitted cell does
        @test _G._modis_cell(_G._modis_cell_center(cell...)...) === cell
    end

    # ------------------------------------------------------------------ geometry handling
    @testset "RGI7 is PolygonZ, and both polygon kinds are accepted" begin
        # The shipped shapefiles are wkbPolygon25D (coorddim 3), not wkbPolygon. A kernel
        # that only matched the 2D type would reject every single RGI7 outline.
        ring = [(0.0, 60.0), (0.1, 60.0), (0.1, 60.1), (0.0, 60.1), (0.0, 60.0)]
        flat = _AG.createpolygon(ring)
        @test _G._rgi7_sinu_geometry(flat, "test") isa _AG.IGeometry

        z = _AG.fromWKT("POLYGON Z ((0 60 10,0.1 60 10,0.1 60.1 10,0 60.1 10,0 60 10))")
        @test _AG.getgeomtype(z) == _AG.wkbPolygon25D
        zs = _G._rgi7_sinu_geometry(z, "test")
        @test zs isa _AG.IGeometry
        # Z is dropped, and the 2D footprint matches the flat version's
        @test isapprox(_AG.geomarea(zs), _AG.geomarea(_G._rgi7_sinu_geometry(flat, "test"));
                       rtol = 1e-12)

        mp = _AG.fromWKT("MULTIPOLYGON (((0 60,0.1 60,0.1 60.1,0 60.1,0 60))," *
                         "((1 60,1.1 60,1.1 60.1,1 60.1,1 60)))")
        @test _G._rgi7_sinu_geometry(mp, "test") isa _AG.IGeometry

        @test_throws ArgumentError _G._rgi7_sinu_geometry(
            _AG.fromWKT("LINESTRING (0 60,1 61)"), "test")
    end

    @testset "antimeridian outlines are rejected, not silently mangled" begin
        # x = R·λ·cos φ with λ ∈ (−180°,180°] sends 179.9 and −179.9 to opposite ends of
        # the grid, so such a ring becomes a planet-spanning polygon burning millions of
        # spurious cells. Fail loud instead.
        bad = _AG.createpolygon([[(179.9, 60.0), (-179.9, 60.0), (-179.9, 61.0),
                                  (179.9, 61.0), (179.9, 60.0)]])
        err = try
            _G._rgi7_sinu_geometry(bad, "RGI2000-v7.0-G-99-99999")
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        @test occursin("RGI2000-v7.0-G-99-99999", err)   # names the glacier
        @test occursin("antimeridian", err)              # says what to do

        # a normal high-latitude outline spanning a few degrees is fine
        ok = _AG.createpolygon([[(-30.0, 80.0), (-25.0, 80.0), (-25.0, 80.5),
                                 (-30.0, 80.5), (-30.0, 80.0)]])
        @test _G._rgi7_sinu_geometry(ok, "fine") isa _AG.IGeometry
    end

    @testset "candidate tiles come from the sinusoidal envelope" begin
        # Taking the envelope in WGS84 and reprojecting its corners would be wrong:
        # sinusoidal does not preserve axis-aligned boxes.
        outlines = _load_fixture_outlines()
        for (id, f) in outlines
            sinu = _G._rgi7_sinu_geometry(f.geom, id)
            tiles = _G._rgi7_candidate_tiles(sinu)
            @test !isempty(tiles)
            @test allunique(tiles)
            @test all(t -> 0 <= t[1] <= _G._MODIS_H_MAX && 0 <= t[2] <= _G._MODIS_V_MAX, tiles)
            # the centroid's own cell must be among the candidates
            @test (_G._modis_cell(f.cenlat, f.cenlon)[1],
                   _G._modis_cell(f.cenlat, f.cenlon)[2]) in tiles
        end
    end

    @testset "closed-form sinusoidal agrees with PROJ" begin
        # The only place PROJ appears in this workflow. Source must be longlat on the MODIS
        # *sphere*: a WGS84 source would make PROJ apply an ellipsoid→sphere datum shift,
        # and the test would measure that shift rather than the projection math.
        src = _AG.importPROJ4("+proj=longlat +a=6371007.181 +b=6371007.181 +no_defs")
        dst = _AG.importPROJ4(_G._MODIS_SINU_PROJ)
        worst = _AG.createcoordtrans(src, dst) do ct
            w = 0.0
            for lat in -85.0:5.0:85.0, lon in -175.0:11.0:175.0
                pt = _AG.createpoint(lon, lat)
                _AG.transform!(pt, ct)
                xc, yc = _G._modis_lonlat_to_sinu(lon, lat)
                w = max(w, hypot(_AG.getx(pt, 0) - xc, _AG.gety(pt, 0) - yc))
            end
            w
        end
        @test worst < 1e-6            # measured 4e-9 m; the docstring claims 0.0 m
    end

    # -------------------------------------------------------------- the vendored artifact
    @testset "vendored cell list" begin
        if !isfile(_G._RGI7_CELLS_CSV) || !isfile(_G._RGI7_GLACIERS_CSV)
            @info "skipping vendored-table tests; generate them with " *
                  "data/make_rgi7_modis_cells.jl" cells = _G._RGI7_CELLS_CSV
            @test_skip false
        else
            t = _G.rgi7_modis_cells()
            @test t === _G.rgi7_modis_cells()          # cached, not re-parsed

            n = length(t)
            ng = _G.n_glaciers(t)
            @test ng == 274_531
            @test all(h -> 0 <= h <= _G._MODIS_H_MAX, t.h)
            @test all(v -> 0 <= v <= _G._MODIS_V_MAX, t.v)
            @test all(r -> 1 <= r <= _G._MODIS_TILE_PIXELS, t.row)
            @test all(c -> 1 <= c <= _G._MODIS_TILE_PIXELS, t.col)

            # Canonical order, and uniqueness of the (cell, glacier) PAIR — not of the cell.
            # A cell repeats only where a sub-cell glacier would otherwise have none, so the
            # albedo point list is rgi7_modis_unique_cells, not the raw columns.
            pair = collect(zip(t.h, t.v, t.row, t.col, t.glacier))
            @test issorted(pair)
            @test allunique(pair)

            uniq = _G.rgi7_modis_unique_cells(t)
            @test issorted(uniq)
            @test allunique(uniq)
            @test length(uniq) <= n
            # sharing is rare; if it ever became common the tie-break would need revisiting
            @test (n - length(uniq)) / n < 0.01
            # every repeated cell must be explained by a forced glacier
            seen = Dict{NTuple{4,Int},Int}()
            for i in 1:n
                c = (Int(t.h[i]), Int(t.v[i]), Int(t.row[i]), Int(t.col[i]))
                if haskey(seen, c)
                    @test t.forced[t.glacier[i]] || t.forced[seen[c]]
                else
                    seen[c] = Int(t.glacier[i])
                end
            end

            # 103 tiles, measured from all 274,531 centroids. A change here means the
            # outlines or the grid math moved; re-derive the number, do not relax the test.
            @test length(unique(zip(t.h, t.v))) == 103

            # every emitted cell round-trips through lon/lat, which is how the albedo API
            # consumes them, and no cell is an off-Earth corner of the sinusoidal grid
            @test all(1:n) do i
                c = (Int(t.h[i]), Int(t.v[i]), Int(t.row[i]), Int(t.col[i]))
                lat, lon = _G._modis_cell_center(c...)
                -90 <= lat <= 90 && -180 <= lon <= 180 && _G._modis_cell(lat, lon) === c
            end

            # glacier bookkeeping
            @test all(g -> 1 <= g <= ng, t.glacier)
            @test sum(t.n_cells) == n
            @test all(>=(1), t.n_cells)                # nobody was silently dropped
            # the forced pass adds one cell and never displaces, so a forced glacier has
            # exactly one and no glacier has zero
            @test all(i -> !t.forced[i] || t.n_cells[i] == 1, 1:ng)
            @test count(t.forced) > 0
            @test issorted(t.rgi_id)                   # binary-searched by rgi7_glacier_cells
            @test all(id -> occursin(r"^RGI2000-v7\.0-G-\d{2}-\d{5}$", id), t.rgi_id)

            # total area: centre-in-polygon is area-unbiased, so cells x cell area should
            # land near the published 706,744 km²
            @test isapprox(length(uniq) * _G._RGI7_CELL_AREA_KM2, 706_744; rtol = 0.10)

            # hemisphere split partitions the cells exactly, by the integer v test
            north, south = _G.rgi7_hemisphere_split(t)
            @test sort(vcat(north, south)) == 1:n
            @test all(i -> t.v[i] <= 8, north)
            @test all(i -> t.v[i] >= 9, south)
            @test !isempty(north) && !isempty(south)

            # and :melt_season must resolve to a real window per half — the whole reason
            # the driver splits, since it returns `nothing` for a straddling list
            latn = [_G._modis_cell_center(Int(t.h[i]), Int(t.v[i]), Int(t.row[i]),
                                          Int(t.col[i]))[1] for i in north[1:min(end, 5000)]]
            lats = [_G._modis_cell_center(Int(t.h[i]), Int(t.v[i]), Int(t.row[i]),
                                          Int(t.col[i]))[1] for i in south[1:min(end, 5000)]]
            @test _G._resolve_doy_range(:melt_season, latn) == _G._MODIS_MELT_SEASON_NORTH
            @test _G._resolve_doy_range(:melt_season, lats) == _G._MODIS_MELT_SEASON_SOUTH
            @test _G._resolve_doy_range(:melt_season, vcat(latn, lats)) === nothing

            # lookup by id returns that glacier's cells and nothing else
            probe = t.rgi_id[fld(ng, 2)]
            idx = _G.rgi7_glacier_cells(t, probe)
            gi = searchsortedfirst(t.rgi_id, probe)
            @test !isempty(idx)
            @test length(idx) == t.n_cells[gi]
            @test all(i -> t.glacier[i] == gi, idx)
        end
    end
end

# ---------------------------------------------------------------- multi-year climatology
#
# Offline and *independent of the vendored tables*: the reduction aligns to whatever
# `RGI7ModisCells` it is handed, so a 40-cell synthetic table exercises every path in
# milliseconds instead of scattering 3.36M rows nine times over. Answers are analytic, so the
# reduction is checked against arithmetic rather than against a previous run of itself.
@testset "RGI7 bare-ice albedo climatology" begin
    # 20 northern (v ≤ 8) and 20 southern (v ≥ 9) cells, in canonical order so the loader's
    # searchsorted placement is exercised exactly as it is on the real table.
    cells_n = [(17, 2, 100 + i, 200) for i in 1:20]
    cells_s = [(17, 12, 100 + i, 200) for i in 1:20]
    allc = sort(vcat(cells_n, cells_s))
    ng = length(allc)
    tbl = _G.RGI7ModisCells(
        UInt8[c[1] for c in allc], UInt8[c[2] for c in allc],
        UInt16[c[3] for c in allc], UInt16[c[4] for c in allc],
        Int32.(1:ng),
        ["RGI2000-v7.0-G-06-" * lpad(i, 5, '0') for i in 1:ng],
        fill(1.0f0, ng), fill(Int32(1), ng), falses(ng))

    uniq = _G.rgi7_modis_unique_cells(tbl)
    @test uniq == allc
    north, south = _G.rgi7_hemisphere_split(tbl)
    @test length(north) == 20 && length(south) == 20
    years = 2001:2005

    # Cell k in the yi-th year gets 0.30 + 0.01*yi, so over five years median == mean == 0.33,
    # minimum 0.31, maximum 0.35. The FIRST cell of each hemisphere is NaN for the first three
    # years, so its n_years is 2 and its median 0.345 — the case proving NaN never reaches
    # `reduction`.
    function write_year(dir, y, yi, cs, hemi)
        open(GzipCompressorStream, _G.rgi7_ice_albedo_path(dir, y, hemi), "w") do io
            println(io, "# synthetic test data, not a retrieval")
            println(io, "h,v,row,col,albedo_bsa,albedo_wsa,n_valid_bsa,n_valid_wsa")
            for (i, c) in enumerate(cs)
                gap = (i == 1 && yi <= 3)
                a = 0.30 + 0.01 * yi
                println(io, join((c[1], c[2], c[3], c[4],
                                  gap ? "NaN" : string(round(a; digits = 4)),
                                  gap ? "NaN" : string(round(a + 0.02; digits = 4)),
                                  gap ? 0 : 40, gap ? 0 : 40), ','))
            end
        end
    end

    mktempdir() do dir
        for (yi, y) in enumerate(years)
            write_year(dir, y, yi, cells_n, :north)
            write_year(dir, y, yi, cells_s, :south)
        end

        clim = _G.rgi7_ice_albedo_climatology(years; dir = dir, cells = tbl, verbose = false)
        @test issetequal(keys(clim), (:albedo_bsa, :albedo_wsa, :n_years_bsa, :n_years_wsa,
                                      :latitude, :longitude))
        @test length(clim[:albedo_bsa]) == ng

        k1n = searchsortedfirst(uniq, cells_n[1])
        k1s = searchsortedfirst(uniq, cells_s[1])
        full = findfirst(i -> i != k1n && i != k1s, 1:ng)

        @test clim[:n_years_bsa][full] == 5
        @test isapprox(clim[:albedo_bsa][full], 0.33f0; atol = 1e-6)
        @test isapprox(clim[:albedo_wsa][full], 0.35f0; atol = 1e-6)

        # NaN years dropped, not propagated
        @test clim[:n_years_bsa][k1n] == 2 && clim[:n_years_bsa][k1s] == 2
        @test isapprox(clim[:albedo_bsa][k1n], 0.345f0; atol = 1e-6)
        @test !isnan(clim[:albedo_bsa][k1n])

        # `reduction` takes a function handle, anonymous ones included
        for (f, expect, name) in ((median, 0.33, "median"), (mean, 0.33, "mean"),
                                  (minimum, 0.31, "minimum"), (maximum, 0.35, "maximum"),
                                  (x -> quantile(x, 0.25), 0.32, "custom"))
            c = _G.rgi7_ice_albedo_climatology(years; dir = dir, cells = tbl,
                                               reduction = f, verbose = false)
            @test isapprox(c[:albedo_bsa][full], Float32(expect); atol = 1e-5)
            @test metadata(c[:albedo_bsa])["reduction"] == name
        end

        # min_years gates the reduction but the count is still reported
        c = _G.rgi7_ice_albedo_climatology(years; dir = dir, cells = tbl,
                                           min_years = 3, verbose = false)
        @test isnan(c[:albedo_bsa][k1n])
        @test c[:n_years_bsa][k1n] == 2
        @test !isnan(c[:albedo_bsa][full])

        # One hemisphere is a supported mode, so it must NOT warn about the other half being
        # absent — a warning fired every year would bury the one that means something.
        # @test_nowarn rather than @test_logs with a LogLevel: that stdlib is not in the test
        # environment either, and "this must not print a warning" is the assertion wanted.
        c = @test_nowarn _G.rgi7_ice_albedo_climatology(
            years; dir = dir, cells = tbl, hemispheres = (:north,), verbose = false)
        @test c[:n_years_bsa][searchsortedfirst(uniq, cells_n[2])] == 5
        @test c[:n_years_bsa][searchsortedfirst(uniq, cells_s[2])] == 0
        @test isnan(c[:albedo_bsa][searchsortedfirst(uniq, cells_s[2])])

        # cell centres carried through, and agreeing with the cell list's own
        latref, lonref = _G.rgi7_modis_cell_points(tbl; index = [full])
        @test clim[:latitude][full] ≈ latref[1]
        @test clim[:longitude][full] ≈ lonref[1]

        # provenance, including the distinction the docstring warns about
        md = metadata(clim)
        @test md["reduction"] == "median"
        @test md["years"] == collect(years)
        @test md["n_cells"] == ng
        @test occursin("not a percentile pooled across years", md["note"])
    end

    # a missing year points at the driver rather than failing obscurely
    mktempdir() do empty_dir
        err = try
            _G.rgi7_ice_albedo_climatology([2099]; dir = empty_dir, cells = tbl, verbose = false)
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing && occursin("run_rgi7_ice_albedo_modis.jl", err)
    end

    # a changed column layout fails loudly rather than transposing albedo into a count
    mktempdir() do bad
        open(GzipCompressorStream, _G.rgi7_ice_albedo_path(bad, 2001, :north), "w") do io
            println(io, "h,v,row,col,albedo_bsa,albedo_wsa,n_valid_bsa")
        end
        err = try
            _G.rgi7_ice_albedo_climatology([2001]; dir = bad, cells = tbl,
                                           hemispheres = (:north,), verbose = false)
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing && occursin("expected", err)
    end

    # validation, all before a file is opened
    @test_throws ArgumentError _G.rgi7_ice_albedo_climatology(Int[]; cells = tbl)
    @test_throws ArgumentError _G.rgi7_ice_albedo_climatology(2001; min_years = 0, cells = tbl)
    @test_throws ArgumentError _G.rgi7_ice_albedo_climatology(2001; hemispheres = (), cells = tbl)
end
