"""
Tests for the C3S satellite surface albedo reader (`satellite_albedo`).

The offline set covers everything that does not need CDS: the 10-daily calendar, the
cross-product cost model and request chunking, request construction, cache keying,
filename parsing, and input validation. These always run.

The integration set orders real data and is opt-in — it submits CDS jobs that queue
server-side for minutes:

    GEMB_TEST_SATELLITE_ALBEDO=1 CDS_API_KEY=... julia --project=. -e 'using Pkg; Pkg.test()'
"""

using Test
using GEMB_ClimateForcing
using Dates
using Rasters
using DimensionalData
import JSON3
import NCDatasets

using GEMB_ClimateForcing: _albedo_nominal_days, _albedo_timesteps, _albedo_validate_years,
    _albedo_normalize_variables, _albedo_cost, _albedo_chunk_timesteps, _albedo_request,
    _albedo_request_hash, _albedo_area, _albedo_area_tag, _albedo_date_from_path,
    _albedo_variable_from_path, _albedo_is_zip, _albedo_era_tag, _default_albedo_cache,
    _ALBEDO_COST_LIMIT, _ALBEDO_YEARS, cds_estimate_cost, _ALBEDO_DEFAULT_LAYERS,
    _ALBEDO_NON_DATA_LAYERS, _albedo_resolve_layer

@testset "Satellite albedo (C3S)" begin

    @testset "Nominal days (10 / 20 / end-of-month)" begin
        @test _albedo_nominal_days(2019, 2) == [10, 20, 28]
        @test _albedo_nominal_days(2020, 2) == [10, 20, 29]   # leap year
        @test _albedo_nominal_days(2019, 6) == [10, 20, 30]
        @test _albedo_nominal_days(2019, 7) == [10, 20, 31]

        # Every day we can ever generate must be one the server accepts. The CDS
        # constraints endpoint reports nominal_day = ["10","20","28","29","30","31"].
        server_days = Set([10, 20, 28, 29, 30, 31])
        generated = Set(d for y in 2018:2024, m in 1:12 for d in _albedo_nominal_days(y, m))
        @test generated ⊆ server_days
    end

    @testset "Timestep enumeration" begin
        ts = _albedo_timesteps((DateTime(2019, 6, 1), DateTime(2019, 6, 30)))
        @test ts == [Date(2019, 6, 10), Date(2019, 6, 20), Date(2019, 6, 30)]

        # Endpoints are inclusive and partial months are trimmed.
        ts = _albedo_timesteps((DateTime(2019, 6, 15), DateTime(2019, 7, 20)))
        @test ts == [Date(2019, 6, 20), Date(2019, 6, 30),
                     Date(2019, 7, 10), Date(2019, 7, 20)]

        # A full year has 36 timesteps (3 per month), and they are sorted and unique.
        year_ts = _albedo_timesteps((DateTime(2019, 1, 1), DateTime(2019, 12, 31)))
        @test length(year_ts) == 36
        @test issorted(year_ts)
        @test allunique(year_ts)

        # nominal_day override replaces the standard cycle.
        ts = _albedo_timesteps((DateTime(2019, 6, 1), DateTime(2019, 6, 30));
                               nominal_day=[15])
        @test ts == [Date(2019, 6, 15)]

        # Reversed range is rejected.
        @test_throws ArgumentError _albedo_timesteps((DateTime(2019, 6, 1), DateTime(2019, 5, 1)))
    end

    @testset "Coverage validation" begin
        @test _albedo_validate_years([Date(2019, 6, 10)]) === nothing
        # Before the Sentinel-3 era (AVHRR/VGT territory, not supported).
        @test_throws ArgumentError _albedo_validate_years([Date(2005, 6, 10)])
        # A range containing no 10-daily timestep at all.
        @test_throws ArgumentError _albedo_timesteps((DateTime(2019, 6, 1), DateTime(2019, 6, 5))) |>
                                  _albedo_validate_years
    end

    @testset "Variable validation" begin
        @test _albedo_normalize_variables(:albb_dh) == [:albb_dh]
        @test _albedo_normalize_variables([:albb_dh, :albb_bh]) == [:albb_dh, :albb_bh]
        @test _albedo_normalize_variables([:albb_dh, :albb_dh]) == [:albb_dh]   # deduped
        @test_throws ArgumentError _albedo_normalize_variables(:not_a_variable)
        @test_throws ArgumentError _albedo_normalize_variables(Symbol[])
        @test Set(SATELLITE_ALBEDO_VARIABLES) ==
              Set((:albb_bh, :albb_dh, :alsp_bh, :alsp_dh))
    end

    @testset "Cost model is the cross product, not the date count" begin
        # Within one month, cost == number of timesteps.
        @test _albedo_cost(1, [Date(2019, 6, 10)]) == 1
        @test _albedo_cost(1, [Date(2019, 6, 10), Date(2019, 6, 20)]) == 2
        @test _albedo_cost(4, [Date(2019, 6, 10), Date(2019, 6, 20)]) == 8

        # Straddling a month boundary costs the cross product (2 months × 2 days = 4),
        # NOT the 2 dates actually wanted. Under-counting here is what gets a request
        # rejected server-side.
        @test _albedo_cost(1, [Date(2019, 5, 31), Date(2019, 6, 10)]) == 4

        # Verified against the live costing endpoint: a full non-leap year spans 12
        # months × 5 distinct nominal days (10, 20, 28, 30, 31 — a non-leap February
        # contributes 28 and never 29), so 1 variable costs 60 and 4 cost 240.
        year_ts = _albedo_timesteps((DateTime(2019, 1, 1), DateTime(2019, 12, 31)))
        @test _albedo_cost(1, year_ts) == 60
        @test _albedo_cost(4, year_ts) == 240

        # A leap year swaps day 28 for 29 rather than adding it, so the distinct-day
        # count — and hence the cost — is unchanged.
        leap_ts = _albedo_timesteps((DateTime(2020, 1, 1), DateTime(2020, 12, 31)))
        @test sort(unique(day.(leap_ts))) == [10, 20, 29, 30, 31]
        @test _albedo_cost(1, leap_ts) == 60

        @test _albedo_cost(1, Date[]) == 0
    end

    @testset "Chunking respects the size limit" begin
        year_ts = _albedo_timesteps((DateTime(2019, 1, 1), DateTime(2019, 12, 31)))

        for n_vars in (1, 2, 4)
            chunks = _albedo_chunk_timesteps(year_ts, n_vars)
            # Every chunk is submittable.
            @test all(_albedo_cost(n_vars, c) <= _ALBEDO_COST_LIMIT for c in chunks)
            # Exact partition: nothing lost, nothing duplicated.
            @test sort(vcat(chunks...)) == year_ts
            @test all(!isempty(c) for c in chunks)
        end

        # Multi-year requests split too.
        five_yr = _albedo_timesteps((DateTime(2019, 1, 1), DateTime(2023, 12, 31)))
        chunks = _albedo_chunk_timesteps(five_yr, 1)
        @test sort(vcat(chunks...)) == five_yr
        @test all(_albedo_cost(1, c) <= _ALBEDO_COST_LIMIT for c in chunks)

        # More variables than the limit can never fit, even for one timestep.
        @test_throws ArgumentError _albedo_chunk_timesteps(year_ts, _ALBEDO_COST_LIMIT + 1)
        @test_throws ArgumentError _albedo_chunk_timesteps(year_ts, 0)
    end

    @testset "Request construction" begin
        dates = [Date(2019, 6, 10), Date(2019, 6, 20)]
        req = _albedo_request([:albb_dh], dates)

        # `sensor` is a bare String while everything else is an array — a CDS quirk.
        @test req["sensor"] isa AbstractString
        @test req["variable"] == ["albb_dh"]
        @test req["satellite"] == ["sentinel_3"]
        @test req["product_version"] == ["v3_1"]
        @test req["horizontal_resolution"] == ["300m"]

        # Zero-padded string selectors.
        @test req["year"] == ["2019"]
        @test req["month"] == ["06"]
        @test req["nominal_day"] == ["10", "20"]

        # No `area` key for a global request.
        @test !haskey(req, "area")

        # `area` is [north, west, south, east].
        area = _albedo_area((X=(-48.0, -47.5), Y=(66.5, 67.0)))
        @test area == [67.0, -48.0, 66.5, -47.5]
        req_area = _albedo_request([:albb_dh], dates; area=area)
        @test req_area["area"] == [67.0, -48.0, 66.5, -47.5]

        # Extents.Extent works as well as a NamedTuple.
        @test _albedo_area(Extent(X=(-48.0, -47.5), Y=(66.5, 67.0))) == area
        @test _albedo_area(nothing) === nothing

        # Out-of-range coordinates are caught locally.
        @test_throws ArgumentError _albedo_area((X=(-200.0, -47.5), Y=(66.5, 67.0)))
        @test_throws ArgumentError _albedo_area((X=(-48.0, -47.5), Y=(66.5, 95.0)))

        # Multiple variables are all carried through.
        multi = _albedo_request([:albb_dh, :albb_bh], dates)
        @test Set(multi["variable"]) == Set(["albb_dh", "albb_bh"])
    end

    @testset "Cache keying is stable" begin
        dates = [Date(2019, 6, 10), Date(2019, 6, 20)]
        req = _albedo_request([:albb_dh], dates)

        # Deterministic across calls (must also be stable across sessions — hence SHA,
        # not Base.hash).
        @test _albedo_request_hash(req) == _albedo_request_hash(_albedo_request([:albb_dh], dates))

        # Invariant to key insertion order and array ordering.
        shuffled = Dict{String,Any}(reverse(collect(req)))
        shuffled["nominal_day"] = reverse(req["nominal_day"])
        @test _albedo_request_hash(shuffled) == _albedo_request_hash(req)

        # Sensitive to what actually changes the data.
        area = _albedo_area((X=(-48.0, -47.5), Y=(66.5, 67.0)))
        @test _albedo_request_hash(_albedo_request([:albb_dh], dates; area=area)) !=
              _albedo_request_hash(req)
        @test _albedo_request_hash(_albedo_request([:albb_bh], dates)) !=
              _albedo_request_hash(req)
        @test _albedo_request_hash(_albedo_request([:albb_dh], [Date(2019, 6, 10)])) !=
              _albedo_request_hash(req)

        # Area tags separate spatial subsets that would otherwise share filenames.
        @test _albedo_area_tag(nothing) == "global"
        @test _albedo_area_tag(area) != _albedo_area_tag(nothing)
        @test _albedo_area_tag(area) == _albedo_area_tag(copy(area))
        @test _albedo_area_tag(area) != _albedo_area_tag([68.0, -48.0, 66.5, -47.5])

        @test startswith(_default_albedo_cache(), tempdir())
        @test _albedo_era_tag() == "sentinel_3_v3_1_300m"
    end

    @testset "Filename parsing" begin
        name = "c3s_ALBB-DH_20190610000000_GLOBE_OLCI_V3.1.nc"
        @test _albedo_date_from_path(name) == Date(2019, 6, 10)
        @test _albedo_variable_from_path(name) == :albb_dh
        @test _albedo_date_from_path("/tmp/some/dir/$(name)") == Date(2019, 6, 10)

        @test _albedo_variable_from_path("c3s_ALBB-BH_20190620000000_GLOBE_OLCI_V3.1.nc") == :albb_bh
        @test _albedo_variable_from_path("c3s_ALSP-DH_20190620000000_GLOBE_OLCI_V3.1.nc") == :alsp_dh
        @test _albedo_variable_from_path("c3s_ALSP-BH_20190620000000_GLOBE_OLCI_V3.1.nc") == :alsp_bh

        # Unparseable names return nothing rather than throwing, so callers can fall back.
        @test _albedo_date_from_path("no_date_here.nc") === nothing
        @test _albedo_variable_from_path("c3s_20190610_GLOBE.nc") === nothing

        # Shuffled paths sort into correct Ti order by parsed date.
        names = ["c3s_ALBB-DH_20190630000000_GLOBE_OLCI_V3.1.nc",
                 "c3s_ALBB-DH_20190610000000_GLOBE_OLCI_V3.1.nc",
                 "c3s_ALBB-DH_20190620000000_GLOBE_OLCI_V3.1.nc"]
        @test sort(_albedo_date_from_path.(names)) ==
              [Date(2019, 6, 10), Date(2019, 6, 20), Date(2019, 6, 30)]
    end

    @testset "Archive magic-byte sniffing" begin
        dir = mktempdir()
        try
            zip_path = joinpath(dir, "a.zip")
            write(zip_path, UInt8[0x50, 0x4b, 0x03, 0x04, 0x00, 0x00])   # "PK\x03\x04"
            @test _albedo_is_zip(zip_path)

            nc3_path = joinpath(dir, "a.nc")
            write(nc3_path, UInt8[0x43, 0x44, 0x46, 0x01, 0x00, 0x00])   # "CDF\x01"
            @test !_albedo_is_zip(nc3_path)

            hdf5_path = joinpath(dir, "b.nc")
            write(hdf5_path, UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a])  # "\x89HDF"
            @test !_albedo_is_zip(hdf5_path)

            # A file too short to hold a signature must not error.
            short_path = joinpath(dir, "short")
            write(short_path, UInt8[0x50, 0x4b])
            @test !_albedo_is_zip(short_path)
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "Layer selection in multi-layer product files" begin
        # One ordered variable arrives as a file of several layers plus a scalar `crs`.
        # Reading it unnamed picks `crs` (no X/Y), which is the bug this guards.
        @test _ALBEDO_DEFAULT_LAYERS[:albb_dh] == :AL_DH_BB
        @test _ALBEDO_DEFAULT_LAYERS[:albb_bh] == :AL_BH_BB
        @test Set(keys(_ALBEDO_DEFAULT_LAYERS)) == Set(SATELLITE_ALBEDO_VARIABLES)
        @test :crs in _ALBEDO_NON_DATA_LAYERS

        dir = mktempdir()
        try
            path = joinpath(dir, "c3s_ALBB-DH_20190610000000_GLOBE_SENTINEL3_V3.1.0.nc")
            NCDatasets.NCDataset(path, "c") do ds
                NCDatasets.defDim(ds, "lon", 3)
                NCDatasets.defDim(ds, "lat", 2)
                NCDatasets.defVar(ds, "lon", [-48.0, -47.9, -47.8], ("lon",))
                NCDatasets.defVar(ds, "lat", [67.0, 66.9], ("lat",))
                for name in ("AL_DH_BB", "AL_DH_BB_ERR", "AL_DH_VI", "QFLAG")
                    NCDatasets.defVar(ds, name, fill(0.5, 3, 2), ("lon", "lat"))
                end
                NCDatasets.defVar(ds, "crs", 'a', ())   # the scalar decoy
            end

            @test Set(satellite_albedo_layers(path)) ==
                  Set([:AL_DH_BB, :AL_DH_BB_ERR, :AL_DH_VI, :QFLAG])
            @test !(:crs in satellite_albedo_layers(path))

            # Default resolves to the headline broadband layer; explicit choices are
            # honoured; unknown ones fail loudly rather than silently reading something else.
            @test _albedo_resolve_layer(path, :albb_dh, nothing) == :AL_DH_BB
            @test _albedo_resolve_layer(path, :albb_dh, :AL_DH_VI) == :AL_DH_VI
            @test_throws ArgumentError _albedo_resolve_layer(path, :albb_dh, :AL_DH_NI)
            @test_throws ArgumentError _albedo_resolve_layer(path, :albb_dh, :crs)

            # If the expected default is absent the fallback warns and picks a real layer.
            @test (@test_logs (:warn,) match_mode=:any _albedo_resolve_layer(
                      path, :albb_bh, nothing)) in satellite_albedo_layers(path)
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "Input validation surfaces before any network call" begin
        # Year outside the Sentinel-3 era: rejected without needing a token.
        @test_throws ArgumentError satellite_albedo(
            time_range=(DateTime(2005, 6, 1), DateTime(2005, 6, 30)),
            extent=Extent(X=(-48.0, -47.5), Y=(66.5, 67.0)), verbose=false)

        # Unknown variable.
        @test_throws ArgumentError satellite_albedo(
            time_range=(DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
            variable=:not_a_variable, verbose=false)

        # Reversed time range.
        @test_throws ArgumentError satellite_albedo(
            time_range=(DateTime(2019, 6, 30), DateTime(2019, 6, 1)), verbose=false)

        # Out-of-range extent.
        @test_throws ArgumentError satellite_albedo(
            time_range=(DateTime(2019, 6, 1), DateTime(2019, 6, 30)),
            extent=Extent(X=(-48.0, -47.5), Y=(66.5, 200.0)), verbose=false)
    end

    # ------------------------------------------------------------------------------
    # Integration: orders real data from CDS. Opt-in — jobs queue for minutes.
    # ------------------------------------------------------------------------------
    if get(ENV, "GEMB_TEST_SATELLITE_ALBEDO", "0") in ("1", "true", "TRUE") &&
       haskey(ENV, "CDS_API_KEY")

        @testset "Integration: order from CDS" begin
            cache = mktempdir()
            token = ENV["CDS_API_KEY"]
            # A small box over the Greenland ice sheet interior — bright, so a missing
            # scale_factor shows up as an out-of-range value.
            extent = Extent(X=(-48.0, -47.5), Y=(66.5, 67.0))
            try
                @testset "Costing endpoint (no job submitted)" begin
                    req = _albedo_request([:albb_dh], [Date(2019, 6, 10)];
                                          area=_albedo_area(extent))
                    cost, limit = cds_estimate_cost("satellite-albedo", req; token=token)
                    @test cost == 1.0
                    @test limit == float(_ALBEDO_COST_LIMIT)
                end

                @testset "Single timestep → lazy RasterSeries" begin
                    series = satellite_albedo(
                        time_range=(DateTime(2019, 6, 10), DateTime(2019, 6, 10)),
                        extent=extent, variable=:albb_dh, cache_path=cache)

                    @test series isa Rasters.AbstractRasterSeries
                    @test length(series) == 1
                    @test hasdim(series, Ti)
                    @test lookup(series, Ti)[1] == Date(2019, 6, 10)

                    layer = series[1]
                    @test layer isa Rasters.AbstractRaster
                    @test hasdim(layer, X)
                    @test hasdim(layer, Y)
                    # Still disk-backed: nothing read yet.
                    @test !(parent(layer) isa Array)
                    @test metadata(layer)["product_version"] == "v3_1"
                    @test metadata(layer)["variable"] == :albb_dh
                    # A real albedo layer was read, not the scalar `crs`.
                    @test metadata(layer)["layer"] == :AL_DH_BB
                    @test Rasters.name(layer) == :AL_DH_BB
                    @test eltype(layer) <: Union{Missing,Real}
                    # The singleton time axis is dropped — Ti lives on the series.
                    @test !hasdim(layer, Ti)
                    @test ndims(layer) == 2

                    store0 = joinpath(cache, "files", _albedo_era_tag(),
                                      _albedo_area_tag(_albedo_area(extent)))
                    nc0 = first(joinpath.(store0,
                        filter(f -> endswith(f, ".nc"), readdir(store0))))
                    @test :AL_DH_BB in satellite_albedo_layers(nc0)
                    @test !(:crs in satellite_albedo_layers(nc0))

                    # A NetCDF landed in the per-timestep store, and the transient job
                    # archive was cleaned up.
                    store = joinpath(cache, "files", _albedo_era_tag(),
                                     _albedo_area_tag(_albedo_area(extent)))
                    ncs = filter(f -> endswith(f, ".nc"), readdir(store))
                    @test length(ncs) >= 1
                    @test !any(f -> endswith(f, ".zip") || f == "download",
                               readdir(joinpath(cache, "jobs"); join=false))

                    # Physically plausible albedo: dimensionless 0–1, and bright over ice.
                    vals = collect(skipmissing(read(layer)))
                    @test !isempty(vals)
                    @test minimum(vals) >= 0.0
                    @test maximum(vals) <= 1.0
                    @test maximum(vals) > 0.3   # catches an unapplied scale_factor
                end

                @testset "Cached: no job resubmitted" begin
                    store = joinpath(cache, "files", _albedo_era_tag(),
                                     _albedo_area_tag(_albedo_area(extent)))
                    files = joinpath.(store, filter(f -> endswith(f, ".nc"), readdir(store)))
                    mtimes_before = mtime.(files)

                    series = satellite_albedo(
                        time_range=(DateTime(2019, 6, 10), DateTime(2019, 6, 10)),
                        extent=extent, variable=:albb_dh, cache_path=cache)

                    @test length(series) == 1
                    @test mtime.(files) == mtimes_before
                end
            finally
                rm(cache; recursive=true, force=true)
            end
        end
    else
        @info "Skipping satellite albedo integration tests. Set GEMB_TEST_SATELLITE_ALBEDO=1 and CDS_API_KEY to enable."
    end
end
