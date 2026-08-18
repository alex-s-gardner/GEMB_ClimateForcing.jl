# Tests for the MODIS MCD43A3 product layer (src/datasets/modis_albedo.jl) and the
# point-based bare-ice albedo driver (src/glacier_ice_albedo_modis.jl).
#
# Everything here is offline except the block gated on GEMB_TEST_MODIS_ALBEDO=1, which
# needs an Earthdata token and downloads a granule.

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData
import ArchGDAL

using GEMB_ClimateForcing: _modis_tile_origin, _modis_tile_bounds, _modis_lonlat_to_sinu,
    _modis_cell, _modis_cell_center, _modis_granule_date, _modis_granule_tile,
    _modis_subdataset, _modis_tile_bbox, _modis_cmr_temporal, _modis_read_window,
    _mcd43a3_quality_layer, _modis_granule_filename, _modis_cached_granule,
    _MODIS_TILE_SPAN_M, _MODIS_PIXEL_M, _MODIS_TILE_PIXELS, _MODIS_UL_X, _MODIS_UL_Y,
    _MODIS_SPHERE_R, _MCD43A3_SCALE, _modis_dedup_points, _low_percentile_mean,
    _LowPercentileTopK, _accumulate!, _finalize

@testset "MODIS MCD43A3 albedo" begin

    @testset "Layer registry" begin
        # Twenty albedo layers: two sky conditions × ten wavelengths.
        @test length(MCD43A3_LAYERS) == 20
        @test :Albedo_BSA_shortwave in MCD43A3_LAYERS
        @test :Albedo_WSA_shortwave in MCD43A3_LAYERS
        @test all(in(MCD43A3_LAYERS), MCD43A3_ALBEDO_LAYERS)
        @test length(MCD43A3_ALBEDO_LAYERS) == 2

        # Black-sky and white-sky share ONE quality band per wavelength.
        @test _mcd43a3_quality_layer(:Albedo_BSA_shortwave) ==
              _mcd43a3_quality_layer(:Albedo_WSA_shortwave) ==
              :BRDF_Albedo_Band_Mandatory_Quality_shortwave
        @test _mcd43a3_quality_layer(:Albedo_WSA_Band3) ==
              :BRDF_Albedo_Band_Mandatory_Quality_Band3
        @test_throws ArgumentError _mcd43a3_quality_layer(:QFLAG)
    end

    @testset "Tile geometry" begin
        # The grid origin is tile (0, 0)'s upper-left corner by definition.
        @test _modis_tile_origin(0, 0) == (_MODIS_UL_X, _MODIS_UL_Y)

        # Hand-computed against the published constants.
        ulx, uly = _modis_tile_origin(17, 2)
        @test ulx ≈ _MODIS_UL_X + 17 * _MODIS_TILE_SPAN_M
        @test uly ≈ _MODIS_UL_Y - 2 * _MODIS_TILE_SPAN_M

        # Adjacent tiles share an edge exactly — a gap or overlap here would put points
        # near a tile boundary in the wrong granule.
        for h in 0:34
            @test _modis_tile_bounds(h, 5)[3] === _modis_tile_bounds(h + 1, 5)[1]
        end
        for v in 0:16
            @test _modis_tile_bounds(9, v)[2] === _modis_tile_bounds(9, v + 1)[4]
        end

        # A tile is exactly 2400 pixels across.
        xmin, ymin, xmax, ymax = _modis_tile_bounds(15, 2)
        @test (xmax - xmin) ≈ _MODIS_TILE_PIXELS * _MODIS_PIXEL_M
        @test (ymax - ymin) ≈ _MODIS_TILE_PIXELS * _MODIS_PIXEL_M

        @test_throws ArgumentError _modis_tile_origin(36, 0)
        @test_throws ArgumentError _modis_tile_origin(0, 18)
        @test_throws ArgumentError _modis_tile_origin(-1, 0)
    end

    @testset "Sinusoidal projection" begin
        # x = R·λ·cos φ, y = R·φ. Validated against PROJ to 0.0 m; asserted here in closed
        # form so the test needs no Proj dependency.
        @test all(_modis_lonlat_to_sinu(0.0, 0.0) .≈ (0.0, 0.0))
        x, y = _modis_lonlat_to_sinu(-50.0, 67.0)
        @test y ≈ _MODIS_SPHERE_R * deg2rad(67.0)
        @test x ≈ _MODIS_SPHERE_R * deg2rad(-50.0) * cos(deg2rad(67.0))
        # Meridian convergence: the same longitude spans less x at higher latitude.
        @test abs(_modis_lonlat_to_sinu(-50.0, 80.0)[1]) <
              abs(_modis_lonlat_to_sinu(-50.0, 20.0)[1])
        # 0–360 longitude is accepted (via _wrap_longitude), like the rest of the package.
        @test all(_modis_lonlat_to_sinu(310.0, 67.0) .≈ _modis_lonlat_to_sinu(-50.0, 67.0))
    end

    @testset "Cell mapping" begin
        # (0, 0) is the north-west corner of tile (18, 9) — the grid's centre seam.
        @test _modis_cell(0.0, 0.0) == (18, 9, 1, 1)

        # A real Greenland point, cross-checked against the granule's GDAL geotransform.
        @test _modis_cell(67.09, -50.05) == (16, 2, 699, 124)

        # Row/col are 1-based and in range everywhere, including the poles and the
        # antimeridian, where the projection is degenerate.
        for (lat, lon) in ((90.0, 0.0), (-90.0, 0.0), (0.0, 180.0), (0.0, -180.0),
                           (89.99, -179.99), (-89.99, 179.99))
            h, v, row, col = _modis_cell(lat, lon)
            @test 0 <= h <= 35 && 0 <= v <= 17
            @test 1 <= row <= _MODIS_TILE_PIXELS && 1 <= col <= _MODIS_TILE_PIXELS
        end

        @test_throws ArgumentError _modis_cell(91.0, 0.0)
        @test_throws ArgumentError _modis_cell(-90.5, 0.0)

        # A nudge well inside the 463 m pixel stays put; one well past it moves. Measured
        # from the CELL CENTRE, not from an arbitrary point: 67.09°N happens to sit within
        # 100 m of a row edge, where a 100 m nudge legitimately does cross into the next row.
        base = _modis_cell(67.09, -50.05)
        clat, clon = _modis_cell_center(base...)
        dlat_100m = 100 / (_MODIS_SPHERE_R * π / 180)
        @test _modis_cell(clat + dlat_100m, clon) == base
        @test _modis_cell(clat + 6 * dlat_100m, clon) != base
        # ...and the move is exactly one row northwards, not an arbitrary jump.
        @test _modis_cell(clat + 6 * dlat_100m, clon)[3] == base[3] - 1

        # Crossing a tile's west edge lands in the neighbouring tile at the far column.
        ulx, _ = _modis_tile_origin(16, 2)
        lat = 67.0
        # Longitude of the tile edge at this latitude, then a pixel either side of it.
        lon_edge = rad2deg(ulx / (_MODIS_SPHERE_R * cos(deg2rad(lat))))
        dlon = rad2deg(_MODIS_PIXEL_M / (_MODIS_SPHERE_R * cos(deg2rad(lat))))
        inside = _modis_cell(lat, lon_edge + dlon / 2)
        outside = _modis_cell(lat, lon_edge - dlon / 2)
        @test inside[1] == 16 && inside[4] == 1
        @test outside[1] == 15 && outside[4] == _MODIS_TILE_PIXELS
        @test inside[3] == outside[3]   # same row: only the tile column changed

        # Cell centres round-trip back into the same cell, within half a pixel.
        for (lat, lon) in ((67.09, -50.05), (-77.5, 163.0), (0.0, 0.0), (45.0, 90.0))
            cell = _modis_cell(lat, lon)
            clat, clon = _modis_cell_center(cell...)
            @test _modis_cell(clat, clon) == cell
            # `<=`, not `<`: a request landing exactly on a row edge (0.0, 45.0 both do) is
            # half a pixel from the centre by construction, which is the bound, not a miss.
            @test abs(clat - lat) <= 0.5 * _MODIS_PIXEL_M / (_MODIS_SPHERE_R * π / 180)
        end
    end

    @testset "Granule identity" begin
        id = "MCD43A3.A2019182.h17v02.061.2019191033013.hdf"
        # Day 182 of 2019 is 1 July — the NOMINAL date, which `time_start` does not give.
        @test _modis_granule_date(id) == Date(2019, 7, 1)
        @test _modis_granule_tile(id) == (17, 2)
        # Works without the extension, as CMR reports producer_granule_id.
        @test _modis_granule_date(id[1:end-4]) == Date(2019, 7, 1)
        # Day 1 and day 366 (2020 was a leap year) both parse.
        @test _modis_granule_date("MCD43A3.A2020001.h17v02.061.x") == Date(2020, 1, 1)
        @test _modis_granule_date("MCD43A3.A2020366.h17v02.061.x") == Date(2020, 12, 31)

        @test_throws ArgumentError _modis_granule_date("MCD43A3.h17v02.061.x.hdf")
        @test_throws ArgumentError _modis_granule_tile("MCD43A3.A2019182.061.x.hdf")
        @test_throws ArgumentError _modis_granule_tile("MCD43A3.A2019182.h99v02.061.x")

        @test _modis_granule_filename("MCD43A3.A2019182.h17v02.061.x") ==
              "MCD43A3.A2019182.h17v02.061.x.hdf"
        @test _modis_granule_filename("a.hdf") == "a.hdf"
    end

    @testset "Cache lookup by nominal date and tile" begin
        # The production-timestamp field of a granule id is unpredictable, so a cached file
        # must be found by the fields that are known.
        mktempdir() do dir
            @test isnothing(_modis_cached_granule(dir, Date(2019, 7, 1), (17, 2)))
            name = "MCD43A3.A2019182.h17v02.061.2019191033013.hdf"
            touch(joinpath(dir, name))
            @test _modis_cached_granule(dir, Date(2019, 7, 1), (17, 2)) ==
                  joinpath(dir, name)
            # A different date or tile must not match, and neither must a .part file.
            @test isnothing(_modis_cached_granule(dir, Date(2019, 7, 2), (17, 2)))
            @test isnothing(_modis_cached_granule(dir, Date(2019, 7, 1), (16, 2)))
            touch(joinpath(dir, "MCD43A3.A2019183.h17v02.061.x.hdf.part"))
            @test isnothing(_modis_cached_granule(dir, Date(2019, 7, 2), (17, 2)))
        end
    end

    @testset "Subdataset spec never uses /vsicurl/" begin
        spec = _modis_subdataset("/data/MCD43A3.A2019182.h17v02.hdf", :Albedo_BSA_shortwave)
        @test spec ==
              "HDF4_EOS:EOS_GRID:\"/data/MCD43A3.A2019182.h17v02.hdf\":MOD_Grid_BRDF:Albedo_BSA_shortwave"
        @test !occursin("/vsicurl", spec)

        # Regression guard for the documented trap: GDAL's HDF4 driver cannot read through
        # a virtual filesystem, and a bare /vsicurl/ open is a FALSE POSITIVE that only
        # fails on first metadata access. Refusing the path here turns that into a clear
        # error instead of an unrelated GDAL message much later.
        @test_throws ArgumentError _modis_subdataset(
            "/vsicurl/https://example.com/g.hdf", :Albedo_BSA_shortwave)
        @test_throws ArgumentError _modis_subdataset(
            "https://example.com/g.hdf", :Albedo_BSA_shortwave)
    end

    @testset "CMR query construction" begin
        # A one-day temporal window is all that is needed — the extra 15 nominal dates it
        # returns are filtered out by the granule id.
        @test _modis_cmr_temporal(Date(2019, 7, 1)) ==
              "2019-07-01T00:00:00Z,2019-07-01T23:59:59Z"

        # The bounding box must CONTAIN the tile: over-selecting is filtered away, but
        # under-selecting silently drops a tile.
        bbox = parse.(Float64, split(_modis_tile_bbox(15, 2), ","))
        w, s, e, n = bbox
        @test w < e && s < n
        @test -180 <= w && e <= 180
        for (lat, lon) in ((67.09, -50.05), (60.5, -45.0))
            _modis_cell(lat, lon)[1:2] == (15, 2) || continue
            @test w <= lon <= e && s <= lat <= n
        end
        # Every tile that has land must produce a valid box.
        for h in 0:35, v in 3:14
            b = parse.(Float64, split(_modis_tile_bbox(h, v), ","))
            @test b[1] <= b[3] && b[2] <= b[4]
        end
    end

    @testset "Point deduplication" begin
        # THE core new behaviour: points sharing a 500 m cell are derived once.
        # Three points within one cell (sub-pixel offsets), then two in another.
        d = 20 / (_MODIS_SPHERE_R * π / 180)     # ~20 m in degrees latitude
        lat = [67.09, 67.09 + d, 67.09 - d, 66.50, 66.50 + d]
        lon = [-50.05, -50.05, -50.05, -49.80, -49.80]

        cells, cell_of_point = _modis_dedup_points(lat, lon)
        @test length(cells) == 2
        @test cell_of_point == [1, 1, 1, 2, 2]
        @test cells[1] == _modis_cell(lat[1], lon[1])
        @test cells[2] == _modis_cell(lat[4], lon[4])

        # Unique cells appear in first-encounter order, so the mapping is deterministic.
        cells_r, map_r = _modis_dedup_points(reverse(lat), reverse(lon))
        @test map_r == [1, 1, 2, 2, 2]
        @test cells_r[1] == cells[2]

        # A single point, and an all-distinct list, are both handled.
        @test _modis_dedup_points([67.0], [-50.0]) == ([_modis_cell(67.0, -50.0)], [1])
        lat3 = [67.0, 68.0, 69.0]
        lon3 = [-50.0, -48.0, -46.0]
        cells3, map3 = _modis_dedup_points(lat3, lon3)
        @test length(cells3) == 3 && map3 == [1, 2, 3]
    end

    @testset "Reduction over a cell vector" begin
        # The point path reuses the C3S streaming top-k accumulator verbatim with an
        # (n_unique, 1) "grid", so it must stay bit-identical to the batch reference.
        n = 7
        n_obs = 36
        obs = [Float32.(0.3 .+ 0.5 .* rand(n, 1)) for _ in 1:n_obs]
        for percentile in (0.05, 0.1, 0.25)
            acc = _LowPercentileTopK((n, 1), n_obs, percentile, 1)
            for o in obs
                _accumulate!(acc, o)
            end
            streamed, counts = _finalize(acc)
            batch = _low_percentile_mean(obs, percentile, 1)
            @test streamed == batch[1]
            @test counts == batch[2]
        end
    end

    @testset "Driver input validation (pre-network)" begin
        lat = [67.09, 66.5]
        lon = [-50.05, -49.8]
        # Nothing below may reach the network: each must throw on argument checking alone.
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, [-50.0], 2019)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(Float64[], Float64[], 2019)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis([91.0], [-50.0], 2019)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis([67.0], [-400.0], 2019)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, Int[])
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 1999)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2100)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   percentile=0.0)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   percentile=1.5)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   min_samples=0)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   albedo_range=(1.0, 0.3))
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   doy_range=(0, 100))
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   doy_range=(100, 400))
        # A window covering the whole year the long way round is a mistake, not a wrap.
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   doy_range=(101, 100))
        # The only symbol accepted is :melt_season; a typo must not silently mean "all".
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   doy_range=:summer)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   stride=0)
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   layers=(:QFLAG,))
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2019;
                                                                   layers=())
        # A doy_range that selects no date in the (partial) first year of the record.
        @test_throws ArgumentError compute_glacier_ice_albedo_modis(lat, lon, 2000;
                                                                   doy_range=(1, 30))
        # The tuple convenience method validates identically.
        @test_throws ArgumentError compute_glacier_ice_albedo_modis([(91.0, -50.0)], 2019)
    end

    @testset "Date list construction" begin
        using GEMB_ClimateForcing: _modis_dates
        # A full year is daily.
        @test length(_modis_dates(2019, nothing, 1)) == 365
        @test length(_modis_dates(2020, nothing, 1)) == 366
        # `doy_range` is the recommended volume knob: inclusive, in day-of-year.
        d = _modis_dates(2019, (120, 290), 1)
        @test length(d) == 290 - 120 + 1
        @test dayofyear(first(d)) == 120 && dayofyear(last(d)) == 290
        # `stride` thins whatever survives, keeping the first date.
        s = _modis_dates(2019, (120, 290), 5)
        @test first(s) == first(d)
        @test length(s) == length(120:5:290)
        # 2000 is truncated to the start of the record, not padded with empty dates.
        y2000 = _modis_dates(2000, nothing, 1)
        @test first(y2000) == Date(2000, 2, 16)
        @test isempty(_modis_dates(2000, (1, 30), 1))
    end

    @testset "Wrapping doy_range (southern hemisphere)" begin
        using GEMB_ClimateForcing: _modis_dates, _modis_doys, _melt_season_doy_range,
            _resolve_doy_range, _MODIS_MELT_SEASON_NORTH, _MODIS_MELT_SEASON_SOUTH

        # `lo > hi` wraps New Year — the southern melt season. Both arcs are present and the
        # count is the sum of their lengths, not a negative difference.
        d = _modis_dates(2019, (300, 110), 1)
        @test length(d) == 110 + (365 - 300 + 1)
        @test dayofyear(first(d)) == 1        # dates come back ascending, not seasonal
        @test dayofyear(last(d)) == 365
        # The gap is the austral winter: nothing between 111 and 299.
        @test !any(x -> 110 < dayofyear(x) <= 299, d)
        @test allunique(d)                    # the two arcs must not overlap
        # Leap years extend the December arc, not the January one.
        @test length(_modis_dates(2020, (300, 110), 1)) == 110 + (366 - 300 + 1)

        # A wrapping window is contiguous *in season*: 31 Dec and 1 Jan are both selected.
        @test Date(2019, 12, 31) in d && Date(2019, 1, 1) in d

        # `stride` thins the selected days, so it applies across the wrap rather than
        # restarting at day 1.
        s = _modis_doys((300, 110), 365)[1:5:end]
        @test _modis_dates(2019, (300, 110), 5) == [Date(2019) + Day(x - 1) for x in s]

        # Non-wrapping behaviour is unchanged.
        @test _modis_doys((120, 290), 365) == collect(120:290)
        @test _modis_doys(nothing, 365) == collect(1:365)
        # `hi` past the end of a short year is clipped, in both forms.
        @test _modis_doys((120, 400), 365) == collect(120:365)
        @test last(_modis_doys((300, 400), 365)) == 365

        # Hemisphere resolution is per point, from the sign of the latitude.
        @test _melt_season_doy_range(67.0) == _MODIS_MELT_SEASON_NORTH
        @test _melt_season_doy_range(-49.3) == _MODIS_MELT_SEASON_SOUTH
        @test _melt_season_doy_range(0.0) == _MODIS_MELT_SEASON_NORTH   # tie → northern
        @test _resolve_doy_range(:melt_season, [67.0, 66.9]) == _MODIS_MELT_SEASON_NORTH
        @test _resolve_doy_range(:melt_season, [-49.3, -50.5]) == _MODIS_MELT_SEASON_SOUTH
        # A list straddling the equator has no single melt season, so it must fall back to
        # the whole year rather than sample the wrong season for half the points.
        @test _resolve_doy_range(:melt_season, [67.0, -49.3]) === nothing
        # Explicit values pass through untouched.
        @test _resolve_doy_range((180, 220), [-49.3]) == (180, 220)
        @test _resolve_doy_range(nothing, [-49.3]) === nothing

        # The southern window really is southern: it must exclude austral midwinter and
        # include midsummer. Getting this backwards is the bug the auto window prevents.
        south = _modis_dates(2019, _MODIS_MELT_SEASON_SOUTH, 1)
        @test Date(2019, 1, 15) in south      # austral high summer
        @test !(Date(2019, 7, 15) in south)   # austral midwinter
        north = _modis_dates(2019, _MODIS_MELT_SEASON_NORTH, 1)
        @test Date(2019, 7, 15) in north
        @test !(Date(2019, 1, 15) in north)
    end
end

# ---------------------------------------------------------------------------- live tests

if get(ENV, "GEMB_TEST_MODIS_ALBEDO", "") == "1"
    @testset "MODIS MCD43A3 (live)" begin
        using GEMB_ClimateForcing: _cmr_granules, mcd43a3_granules, MCD43A3_SHORT_NAME,
            MCD43A3_VERSION, _modis_assert_hdf4_driver,
            _granule_id, _granule_data_url, _granule_bytes
        import EarthData

        # Fail loudly and early if the running GDAL cannot read HDF4 at all.
        @test isnothing(_modis_assert_hdf4_driver())

        date = Date(2019, 7, 1)
        tile = (16, 2)

        # CMR returns the 16-day window's worth of granules; exactly one survives the
        # nominal-date + tile filter.
        granules = _cmr_granules(; short_name=MCD43A3_SHORT_NAME, version=MCD43A3_VERSION,
                                 temporal=_modis_cmr_temporal(date),
                                 bounding_box=_modis_tile_bbox(tile...), verbose=false)
        @test !isempty(granules)
        matching = filter(g -> _modis_granule_date(g.id) == date &&
                               _modis_granule_tile(g.id) == tile, granules)
        @test length(matching) == 1
        # The trap, live: more nominal dates come back than were asked for.
        @test length(unique(_modis_granule_date.(getproperty.(granules, :id)))) > 1

        # `_granule_id` reads `DataGranule.Identifiers`' ProducerGranuleId; for MCD43A3 that
        # equals `GranuleUR`, but the equality is a per-provider convention rather than a
        # schema guarantee, so it is asserted here against live records instead of assumed.
        # If LP DAAC ever diverges, the typed identifier is still the right field and the
        # `A<YYYYDDD>` / `h##v##` filters keep working — this test just stops being silent.
        raw = EarthData.request(EarthData.granule_url,
            Dict("short_name" => MCD43A3_SHORT_NAME, "version" => MCD43A3_VERSION,
                 "temporal" => _modis_cmr_temporal(date),
                 "bounding_box" => _modis_tile_bbox(tile...)),
            EarthData.Granules.UMM_G; page_size=10)
        @test !isempty(raw)
        for g in raw
            @test _granule_id(g) == String(g.GranuleUR)
            # And the id is the archive filename, which is what the filters parse.
            @test occursin(r"^MCD43A3\.A\d{7}\.h\d\dv\d\d\.061\.\d+$", _granule_id(g))
            # A `GET DATA` URL exists, and the size has been converted out of `SizeUnit`
            # rather than taken literally as bytes. The bound is deliberately loose at the
            # bottom: measured over one tile-day, granules run 1.04 MB (h15v03, almost all
            # ocean) to 87 MB (h16v01), so a "~70 MB granule" floor would be wrong. What is
            # being asserted is only that MB was not read as B — an unconverted 53.9037
            # rounds to 54.
            @test !isnothing(_granule_data_url(g))
            @test 100_000 < _granule_bytes(g) < 500_000_000
        end

        mktempdir() do cache
            paths = mcd43a3_granules(date, [tile]; cache_path=cache, verbose=false)
            @test haskey(paths, tile)
            path = paths[tile]
            @test isfile(path) && filesize(path) > 10_000_000
            # No .part file survives a successful download.
            @test !isfile(path * ".part")

            # THE highest-value live assertion: the closed-form tile arithmetic the whole
            # point-sampling design rests on, checked against the granule's own
            # geotransform.
            ArchGDAL.read(_modis_subdataset(path, :Albedo_BSA_shortwave)) do ds
                gt = ArchGDAL.getgeotransform(ds)
                ulx, uly = _modis_tile_origin(tile...)
                @test abs(gt[1] - ulx) < 1.0
                @test abs(gt[4] - uly) < 1.0
                @test abs(gt[2] - _MODIS_PIXEL_M) < 1e-6
                @test abs(-gt[6] - _MODIS_PIXEL_M) < 1e-6
                @test ArchGDAL.width(ds) == ArchGDAL.height(ds) == _MODIS_TILE_PIXELS
            end

            # Every documented layer is present, and values are in the raw Int16 domain —
            # GDAL reports the 0.001 scale but does NOT apply it on read.
            for layer in MCD43A3_LAYERS
                w = _modis_read_window(path, layer, 700:701, 120:121)
                @test eltype(w) === Int16
                @test all(v -> 0 <= v <= 32767, w)
            end
            qa = _modis_read_window(path, _mcd43a3_quality_layer(:Albedo_BSA_shortwave),
                                    700:701, 120:121)
            @test eltype(qa) === UInt8

            # A bare /vsicurl/ open of HDF4 is a FALSE POSITIVE: it "succeeds" and then
            # fails on first metadata access. Documented here so nobody removes the
            # download step.
            remote = "/vsicurl/" * granules[1].url
            @test_throws Exception ArchGDAL.read(
                "HDF4_EOS:EOS_GRID:\"$(remote)\":MOD_Grid_BRDF:Albedo_BSA_shortwave") do ds
                ArchGDAL.getgeotransform(ds)
            end

            # End to end, on points that deliberately share a cell.
            d = 20 / (_MODIS_SPHERE_R * π / 180)
            lat = [67.09, 67.09 + d, 66.95]
            lon = [-50.05, -50.05, -49.90]
            ice = compute_glacier_ice_albedo_modis(lat, lon, 2019;
                doy_range=(180, 190), min_samples=1, cache_path=cache,
                keep_granules=true, verbose=false, progress=false)

            @test issetequal(keys(ice), (:albedo_bsa, :albedo_wsa,
                                         :n_valid_observations_bsa,
                                         :n_valid_observations_wsa,
                                         :latitude, :longitude, :cell_id))
            bsa = ice[:albedo_bsa]
            @test size(bsa) == (3, 1)
            @test all(v -> isnan(v) || 0.3 <= v <= 1.0, bsa)
            @test all(<=(11), ice[:n_valid_observations_bsa])

            # Points in the same cell share their cell id AND their value, bit-for-bit.
            @test ice[:cell_id][1] == ice[:cell_id][2]
            @test ice[:cell_id][1] != ice[:cell_id][3]
            @test isequal(bsa[1, 1], bsa[2, 1])
            @test isequal(ice[:albedo_wsa][1, 1], ice[:albedo_wsa][2, 1])
            @test ice[:n_valid_observations_bsa][1, 1] ==
                  ice[:n_valid_observations_bsa][2, 1]

            # Black-sky is generally darker than white-sky over ice; a sanity check on
            # layer identity, not a physical law.
            wsa = ice[:albedo_wsa]
            both = [(bsa[i], wsa[i]) for i in eachindex(bsa)
                    if !isnan(bsa[i]) && !isnan(wsa[i])]
            isempty(both) || @test count(p -> p[1] <= p[2], both) >= length(both) / 2

            # A second identical call hits the per-(date, tile) cache and is bit-identical.
            ice2 = compute_glacier_ice_albedo_modis(lat, lon, 2019;
                doy_range=(180, 190), min_samples=1, cache_path=cache,
                keep_granules=true, verbose=false, progress=false)
            @test isequal(parent(ice2[:albedo_bsa]), parent(bsa))
            @test ice2[:n_valid_observations_bsa] == ice[:n_valid_observations_bsa]
        end

        # `keep_granules=false` is the DEFAULT and the only path that removes files, so it
        # gets its own cache directory: no granule may survive, the sample cache must, and
        # the values must match the keep_granules=true run above.
        mktempdir() do cache
            d = 20 / (_MODIS_SPHERE_R * π / 180)
            ice = compute_glacier_ice_albedo_modis([67.09, 67.09 + d], [-50.05, -50.05],
                2019; doy_range=(180, 182), min_samples=1, cache_path=cache,
                verbose=false, progress=false)
            @test isempty(filter(endswith(".hdf"), readdir(cache)))
            @test !isempty(readdir(joinpath(cache, "samples")))
            @test isequal(ice[:albedo_bsa][1, 1], ice[:albedo_bsa][2, 1])
            # Having discarded the granules, a re-run must still complete from the sample
            # cache alone — that is what makes the destructive default affordable.
            ice2 = compute_glacier_ice_albedo_modis([67.09, 67.09 + d], [-50.05, -50.05],
                2019; doy_range=(180, 182), min_samples=1, cache_path=cache,
                verbose=false, progress=false)
            @test isequal(parent(ice2[:albedo_bsa]), parent(ice[:albedo_bsa]))
        end
    end
else
    @info "Skipping live MODIS albedo tests. Set GEMB_TEST_MODIS_ALBEDO=1 and EARTHDATA_TOKEN to enable."
end
