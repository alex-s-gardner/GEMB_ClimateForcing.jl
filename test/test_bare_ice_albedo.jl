using Test
using GEMB_ClimateForcing
using DimensionalData
using Rasters
using Statistics
using CodecZlib
import GeoInterface as GI

const G = GEMB_ClimateForcing

# A synthetic pooled table, so the offline tests do not depend on the vendored run being
# present. Written in the same canonical cell order and with the same header as
# `data/run_rgi7_pooled_albedo.jl`'s `write_pooled`, because the reader validates neither — it
# locates rows by `searchsortedfirst` on the cell list and would silently mis-read a
# differently ordered file.
function _write_test_pooled(dir, hemi, rows)
    path = bare_ice_albedo_path(dir, hemi)
    open(GzipCompressorStream, path, "w") do io
        println(io, "# synthetic test fixture")
        println(io, "h,v,row,col,albedo_bsa,albedo_wsa,n_valid_bsa,n_valid_wsa,k_used_bsa,k_used_wsa")
        for (cell, bsa, wsa, nb, nw, kb, kw) in rows
            println(io, join((cell..., bsa, wsa, nb, nw, kb, kw), ','))
        end
    end
    return path
end

@testset "bare_ice_albedo" begin

    @testset "geometry classification and projection" begin
        pt = GI.Point(-30.0, 70.0)
        @test G._bia_shape(pt) === :point
        @test G._bia_shape(GI.MultiPoint([(-30.0, 70.0), (-30.1, 70.1)])) === :point
        @test G._bia_shape(GI.LineString([(-30.0, 70.0), (-29.9, 70.1)])) === :line
        ring = GI.LinearRing([(-30.0, 70.0), (-29.9, 70.0), (-29.9, 70.1), (-30.0, 70.0)])
        @test G._bia_shape(GI.Polygon([ring])) === :polygon
        @test G._bia_shape(GI.MultiPolygon([GI.Polygon([ring])])) === :polygon

        # Projection must be the same closed form the vendored cell table was burned with, or a
        # burned cell and a looked-up cell can disagree.
        sp = G._bia_project(pt)
        @test GI.geomtrait(sp) isa GI.PointTrait
        @test (GI.x(sp), GI.y(sp)) === G._modis_lonlat_to_sinu(-30.0, 70.0)
        # Structure is preserved, including interior rings.
        hole = GI.LinearRing([(-29.98, 70.02), (-29.95, 70.02), (-29.95, 70.05), (-29.98, 70.02)])
        pp = G._bia_project(GI.Polygon([ring, hole]))
        @test GI.geomtrait(pp) isa GI.PolygonTrait
        @test GI.nring(pp) == 2
    end

    @testset "a bare point is not iterable as points" begin
        # GI.getpoint has no PointTrait method; _bia_points exists solely for this.
        @test length(collect(G._bia_points(GI.Point(1.0, 2.0)))) == 1
        @test length(collect(G._bia_points(GI.LineString([(1.0, 2.0), (3.0, 4.0)])))) == 2
    end

    @testset "burn target is bit-identical to the cell grid" begin
        # The whole design rests on this: boolmask evaluates its inside test at the lookup
        # values, so those must be the exact cell centres `_modis_cell_center` reports.
        #
        # The window slice is pure arithmetic and holds for any tile. The geographic round trip
        # is asserted only for cells inside the projection's valid domain: away from the
        # equator a sinusoidal parallel is far shorter than the 36-tile grid is wide, so most
        # of a high-`v` tile row projects from `|lon| > 180` and has no geographic preimage.
        # Those tiles hold no data and never reach `_bia_burn`.
        for (h, v) in ((17, 2), (0, 0), (35, 17), (16, 1), (21, 9), (12, 15))
            dx, dy = G._rgi7_tile_dims(h, v)
            xs, ys = parent(parent(dx)), parent(parent(dy))
            wx, wy = G._bia_window_dims(h, v, 10:14, 20:25)
            @test parent(parent(wx)) == xs[20:25]
            @test parent(parent(wy)) == ys[10:14]
            for r in 10:14, c in 20:25
                lat, lon = G._modis_cell_center(h, v, r, c)
                abs(lon) <= 180 || continue
                x, y = G._modis_lonlat_to_sinu(lon, lat)
                @test x ≈ xs[c] atol = 1e-6
                @test y ≈ ys[r] atol = 1e-6
                @test G._modis_cell(lat, lon) === (h, v, r, c)
            end
        end
    end

    @testset "round trip holds at every real glacier cell" begin
        # The in-domain guard above must not be hiding a real failure, so pin the identity on
        # cells that actually carry data. Tiles are taken from the vendored table rather than
        # written down: a hand-picked list is easy to get wrong, since a tile can be entirely
        # outside the projection at its own latitude (h01v12 spans 30–40°S, where a sinusoidal
        # parallel is only 16.4 Mm wide against that tile's 17.8 Mm of x).
        cellfile = joinpath(@__DIR__, "..", "data", "rgi7_modis_cells.csv.gz")
        if isfile(cellfile)
            sampled = NTuple{4,Int}[]
            open(GzipDecompressorStream, cellfile) do io
                k = 0
                for line in eachline(io)
                    (isempty(line) || startswith(line, '#') || startswith(line, "h,")) && continue
                    k += 1
                    # Every 40000th row, so the sample spans all 103 tiles cheaply.
                    k % 40_000 == 0 || continue
                    f = split(line, ',')
                    push!(sampled, (parse(Int, f[1]), parse(Int, f[2]),
                                    parse(Int, f[3]), parse(Int, f[4])))
                end
            end
            @test length(sampled) > 50
            @test length(unique((c[1], c[2]) for c in sampled)) > 20
            for cell in sampled
                lat, lon = G._modis_cell_center(cell...)
                @test abs(lon) <= 180
                @test G._modis_cell(lat, lon) === cell
                # ... and the burn puts a point at that centre back in the same cell.
                @test G._bia_burn(G._bia_project(GI.Point(lon, lat)), :point, :center) == [cell]
            end
        else
            @info "skipping real-cell round trip; run data/make_rgi7_modis_cells.jl first"
        end
    end

    @testset "a point selects the cell containing it" begin
        for cell in ((17, 2, 20, 460), (17, 2, 1, 1), (23, 5, 2400, 2400), (1, 9, 1200, 700))
            lat, lon = G._modis_cell_center(cell...)
            burned = G._bia_burn(G._bia_project(GI.Point(lon, lat)), :point, :center)
            @test burned == [cell]
        end
    end

    @testset "burn agrees with Rasters.extract where extract works" begin
        # extract cannot index a 3-D (X, Y, Ti) raster and throws on a multipoint in Rasters
        # 0.15, which is why boolmask is the primitive — but where both work they must agree.
        cell = (17, 2, 20, 460)
        lat, lon = G._modis_cell_center(cell...)
        ring = GI.LinearRing([(lon - 0.12, lat - 0.06), (lon + 0.12, lat - 0.06),
                              (lon + 0.12, lat + 0.06), (lon - 0.12, lat + 0.06),
                              (lon - 0.12, lat - 0.06)])
        cases = ((GI.Point(lon, lat), :point),
                 (GI.LineString([(lon - 0.15, lat - 0.05), (lon + 0.15, lat + 0.05)]), :line),
                 (GI.Polygon([ring]), :polygon))
        for (geom, shape) in cases
            sinu = G._bia_project(geom)
            mine = G._bia_burn(sinu, shape, :center)
            theirs = NTuple{4,Int}[]
            for (h, v, rows, cols) in G._bia_tile_windows(sinu)
                dims = G._bia_window_dims(h, v, rows, cols)
                tmpl = Raster(zeros(Int, length(cols), length(rows)), dims; missingval = missing)
                got = extract(tmpl, sinu; index = true, geometry = false, skipmissing = false)
                got = got isa NamedTuple ? [got] : collect(got)
                for r in got
                    r.index === missing && continue
                    push!(theirs, (h, v, rows[r.index[2]], cols[r.index[1]]))
                end
            end
            @test mine == sort(theirs)
            @test !isempty(mine)
        end
    end

    @testset "multipoint burns, where extract throws" begin
        c1, c2 = (17, 2, 20, 460), (17, 2, 1, 546)
        p1 = G._modis_cell_center(c1...)
        p2 = G._modis_cell_center(c2...)
        mp = GI.MultiPoint([(p1[2], p1[1]), (p2[2], p2[1])])
        @test G._bia_burn(G._bia_project(mp), :point, :center) == sort([c1, c2])
        # Two points in one cell select it once, matching _modis_dedup_points.
        dup = GI.MultiPoint([(p1[2], p1[1]), (p1[2], p1[1])])
        @test G._bia_burn(G._bia_project(dup), :point, :center) == [c1]
    end

    @testset "cross-tile geometry, no duplicates" begin
        # The h17|h18 edge is x = 0, i.e. lon = 0, at any latitude.
        box = GI.Polygon([GI.LinearRing([(-0.05, 69.9), (0.05, 69.9),
                                         (0.05, 69.95), (-0.05, 69.95), (-0.05, 69.9)])])
        cells = G._bia_burn(G._bia_project(box), :polygon, :center)
        @test length(unique(cells)) == length(cells)
        @test Set((c[1], c[2]) for c in cells) == Set([(17, 2), (18, 2)])
        @test issorted(cells)
    end

    @testset "boundary widens the selection monotonically" begin
        cell = (17, 2, 20, 460)
        lat, lon = G._modis_cell_center(cell...)
        ring = GI.LinearRing([(lon - 0.05, lat - 0.03), (lon + 0.05, lat - 0.03),
                              (lon + 0.05, lat + 0.03), (lon - 0.05, lat + 0.03),
                              (lon - 0.05, lat - 0.03)])
        poly = GI.Polygon([ring])
        sinu = G._bia_project(poly)
        inside = G._bia_burn(sinu, :polygon, :inside)
        centre = G._bia_burn(sinu, :polygon, :center)
        touches = G._bia_burn(sinu, :polygon, :touches)
        @test length(inside) <= length(centre) <= length(touches)
        @test issubset(Set(centre), Set(touches))
    end

    @testset "hemisphere split is exact at the equator" begin
        @test G._modis_tile_edge_y(9) === 0.0
        @test all(G._bia_hemisphere(v) === :north for v in 0:8)
        @test all(G._bia_hemisphere(v) === :south for v in 9:17)
    end

    @testset "field parser" begin
        buf = Vector{UInt8}("17,2,1632,2171,NaN,NaN,0,0")
        stop = length(buf)
        h, i = G._bia_next_uint(buf, 1, stop)
        v, i = G._bia_next_uint(buf, i, stop)
        r, i = G._bia_next_uint(buf, i, stop)
        c, i = G._bia_next_uint(buf, i, stop)
        @test (h, v, r, c) === (17, 2, 1632, 2171)
        @test String(buf[i:stop]) == "NaN,NaN,0,0"
        # Multi-digit and single-digit fields, and a value at the very end.
        b2 = Vector{UInt8}("2400")
        @test G._bia_next_uint(b2, 1, length(b2))[1] == 2400
    end

    @testset "input validation" begin
        pt = GI.Point(-30.0, 70.0)
        @test_throws "expected a GeoInterface geometry" bare_ice_albedo("not a geometry")
        @test_throws "boundary must be :center, :touches or :inside" bare_ice_albedo(pt; boundary = :middle)
        @test_throws "min_cells must be >= 1" bare_ice_albedo(pt, median; min_cells = 0)
        @test_throws "cannot represent as one shape" bare_ice_albedo(
            GI.LineString([(-179.0, -70.0), (179.0, -70.0)]))
        # A missing pooled table must say how to make it, and must say that doing so needs no
        # download — that is the whole reason the sample cache is kept.
        mktempdir() do empty
            @test_throws "no pooled bare-ice albedo file" bare_ice_albedo(pt; dir = empty)
            @test_throws "downloads nothing" bare_ice_albedo(pt; dir = empty)
        end
    end

    @testset "read and metadata, against a synthetic pooled table" begin
        mktempdir() do dir
            # Three cells in h17v02; the middle one is present but unresolved (no retrieval
            # survived QC), and cells elsewhere in the query polygon are absent from the
            # product entirely.
            a, b, c = (17, 2, 100, 100), (17, 2, 100, 101), (17, 2, 100, 102)
            _write_test_pooled(dir, :north,
                               [(a, "0.4000", "0.3900", 400, 410, 20, 21),
                                (b, "NaN", "NaN", 0, 0, 0, 0),
                                (c, "0.6000", "0.5900", 500, 500, 25, 25)])

            latb, lonb = G._modis_cell_center(b...)
            pt = bare_ice_albedo(GI.Point(lonb, latb); dir = dir)
            @test dims(pt) == (Dim{:cell}(1:1),)
            @test isnan(pt[:albedo_bsa][1])
            # An unresolved cell is still in the product: that is what distinguishes it from
            # an off-glacier cell.
            @test pt[:in_product][1]
            @test pt[:n_valid_bsa][1] == 0
            @test pt[:k_used_bsa][1] == 0

            lata, lona = G._modis_cell_center(a...)
            at = bare_ice_albedo(GI.Point(lona, lata); dir = dir)
            @test at[:albedo_bsa][1] === 0.4f0
            @test at[:albedo_wsa][1] === 0.39f0
            @test at[:n_valid_bsa][1] == 400
            @test at[:n_valid_wsa][1] == 410
            @test at[:k_used_bsa][1] == 20
            @test at[:cell_id][1] == G._modis_cell_id(a)

            # A cell off the RGI outlines is selected but flagged, not silently NaN-as-cloud.
            lat0, lon0 = G._modis_cell_center(17, 2, 500, 500)
            off = (@test_logs (:warn, r"off-glacier") bare_ice_albedo(
                GI.Point(lon0, lat0); dir = dir))
            @test !off[:in_product][1]
            @test DimensionalData.metadata(off)["n_cells_in_product"] == 0

            # A line across all three cells keeps them in canonical order.
            latc, lonc = G._modis_cell_center(c...)
            ln = GI.LineString([(lona, lata), (lonc, latc)])
            per = bare_ice_albedo(ln; dir = dir)
            @test collect(per[:cell_id]) == [G._modis_cell_id(x) for x in (a, b, c)]
            @test collect(per[:latitude]) ≈ [G._modis_cell_center(x...)[1] for x in (a, b, c)]
            @test collect(per[:in_product]) == [true, true, true]

            md = DimensionalData.metadata(per)
            @test md["n_cells"] == 3
            @test md["shape"] == "line"
            @test md["boundary"] == "center"
            @test md["grid"] == G._MODIS_SINU_PROJ
        end
    end

    @testset "reducer returns scalars across cells" begin
        mktempdir() do dir
            a, b, c = (17, 2, 100, 100), (17, 2, 100, 101), (17, 2, 100, 102)
            _write_test_pooled(dir, :north,
                               [(a, "0.4000", "0.4000", 400, 400, 20, 20),
                                (b, "NaN", "NaN", 0, 0, 0, 0),
                                (c, "0.6000", "0.6000", 500, 500, 25, 25)])
            lata, lona = G._modis_cell_center(a...)
            latc, lonc = G._modis_cell_center(c...)
            ln = GI.LineString([(lona, lata), (lonc, latc)])

            for (f, want) in ((minimum, 0.4f0), (maximum, 0.6f0), (mean, 0.5f0),
                              (median, 0.5f0), (x -> quantile(x, 0.0), 0.4f0))
                r = bare_ice_albedo(ln, f; dir = dir)
                @test r isa NamedTuple
                @test r.albedo_bsa ≈ want
                # The reducer never sees the unresolved cell, which is what makes the count
                # meaningful.
                @test r.n_cells_bsa == 2
                @test r.n_cells == 3
                @test r.n_valid_bsa == 900
                @test r.k_used_bsa == 45
            end
            @test bare_ice_albedo(ln, median; dir = dir).reduction == "median"
            @test bare_ice_albedo(ln, x -> quantile(x, 0.5); dir = dir).reduction == "custom"

            # min_cells reports too-thin a sample rather than reducing over almost nothing.
            @test isnan(bare_ice_albedo(ln, median; dir = dir, min_cells = 3).albedo_bsa)
        end
    end

    @testset "southern geometry reads only the southern table" begin
        mktempdir() do dir
            s = (17, 12, 500, 500)
            _write_test_pooled(dir, :south, [(s, "0.7000", "0.6900", 300, 300, 15, 15)])
            lat, lon = G._modis_cell_center(s...)
            # No northern table exists in `dir`; a northern read would throw.
            r = bare_ice_albedo(GI.Point(lon, lat); dir = dir)
            @test r[:albedo_bsa][1] === 0.7f0
            @test r[:n_valid_bsa][1] == 300
        end
    end

    @testset "pooled product, if present" begin
        dir = joinpath(@__DIR__, "..", "data")
        if isfile(bare_ice_albedo_path(dir, :north))
            cell = (17, 2, 20, 460)
            lat, lon = G._modis_cell_center(cell...)
            r = bare_ice_albedo(GI.Point(lon, lat); dir = dir)
            @test r[:cell_id][1] == "h17v02_r0020_c0460"
            @test r[:in_product][1]
            # One value per cell, no time dimension, and the counts qualify it.
            @test dims(r) == (Dim{:cell}(1:1),)
            n = r[:n_valid_bsa][1]
            @test n >= 0
            if n > 0
                @test !isnan(r[:albedo_bsa][1])
                @test r[:albedo_bsa][1] >= POOLED_ICE_ALBEDO_RANGE[1]
                @test r[:k_used_bsa][1] == G._pooled_k_used(n, 0.05)
            end

            # A reduction agrees with reducing the per-cell stack by hand.
            ring = GI.LinearRing([(lon - 0.12, lat - 0.06), (lon + 0.12, lat - 0.06),
                                  (lon + 0.12, lat + 0.06), (lon - 0.12, lat + 0.06),
                                  (lon - 0.12, lat - 0.06)])
            poly = GI.Polygon([ring])
            per = bare_ice_albedo(poly; dir = dir)
            red = bare_ice_albedo(poly, median; dir = dir)
            keep = filter(!isnan, collect(per[:albedo_bsa]))
            if !isempty(keep)
                @test red.albedo_bsa === Float32(median(keep))
                @test red.n_cells_bsa == length(keep)
            end
        else
            @info "skipping pooled-product testset; run data/run_rgi7_pooled_albedo.jl first"
        end
    end
end
