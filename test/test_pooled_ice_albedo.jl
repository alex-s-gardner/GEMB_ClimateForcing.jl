using Test
using GEMB_ClimateForcing
using DimensionalData
using Dates
using Statistics

const G = GEMB_ClimateForcing

using GEMB_ClimateForcing: _modis_cells_key, _modis_sample_cache_path, _pooled_k_used,
    _pooled_cached_dates, _pooled_layer_names, _MCD43A3_SCALE

const _CACHE_COLS = [:Albedo_BSA_shortwave, :Albedo_WSA_shortwave,
                     :BRDF_Albedo_Band_Mandatory_Quality_shortwave]

"""
    _write_cache(root, cells, rows_by_date) -> String

Build a synthetic sample cache under `root` for `cells`.

`rows_by_date` maps a `Date` to one `(bsa_dn, wsa_dn, qa)` tuple per cell, in cell order — the
raw form the cache actually holds. The cache is keyed by a SHA of the whole cell list, so the
key is derived here the same way the reader derives it.
"""
function _write_cache(root, cells, rows_by_date)
    key = _modis_cells_key(cells)
    dir = joinpath(root, "samples", key)
    mkpath(dir)
    for (date, rows) in rows_by_date
        open(joinpath(dir, string(date, ".tsv")), "w") do io
            println(io, join(_CACHE_COLS, "\t"))
            for r in rows
                println(io, join(r, "\t"))
            end
        end
    end
    return key
end

# The statistic, computed independently: mean of the darkest ceil(pct*n) of all pooled values.
function _brute(vals, pct)
    isempty(vals) && return (NaN32, 0, 0)
    s = sort(Float32.(vals))
    n = length(s)
    k = max(1, ceil(Int, pct * n))
    return (Float32(mean(view(s, 1:k))), n, k)
end

@testset "pooled bare-ice albedo" begin

    @testset "k_used is exactly what _finalize_into! averages" begin
        # Derived rather than returned from the accumulator, so the identity has to hold for
        # every count — a wrong k_used would misreport how much evidence is behind a value.
        for pct in (0.05, 0.1, 0.25, 1.0), n in 0:400
            @test _pooled_k_used(n, pct) == (n < 1 ? 0 : max(1, ceil(Int, pct * n)))
        end
        @test _pooled_k_used(0, 0.05) == 0
        @test _pooled_k_used(1, 0.05) == 1        # never 0 for a cell with any retrieval
        @test _pooled_k_used(4446, 0.05) == 223   # the northern record, i.e. kmax
    end

    @testset "layer naming matches the per-year stack" begin
        @test _pooled_layer_names(:Albedo_BSA_shortwave) ==
            (:albedo_bsa, :n_valid_observations_bsa, :k_used_bsa)
        @test _pooled_layer_names(:Albedo_WSA_shortwave) ==
            (:albedo_wsa, :n_valid_observations_wsa, :k_used_wsa)
    end

    @testset "pooled == brute force over all dates" begin
        mktempdir() do root
            cells = sort([(17, 2, 100, 100), (17, 2, 100, 101), (17, 2, 100, 102)])
            dates = Date(2001, 1, 1) .+ Day.(0:39)
            truth1 = Float64[]
            truth3 = Float64[]
            rows_by_date = Pair{Date,Vector{Tuple{Int,Int,Int}}}[]
            for (i, d) in enumerate(dates)
                a1 = 300 + 10 * (i - 1)          # cell 1: 0.30 … 0.69, QA 0 (kept)
                a3 = 200 + 7 * (i - 1)           # cell 3: straddles the 0.25 floor
                q3 = iseven(i) ? 0 : 1           # ... and alternates QA 0 / QA 1 (rejected)
                push!(truth1, a1 * 0.001)
                if q3 == 0 && 0.25 <= a3 * 0.001 <= 1.0
                    push!(truth3, a3 * 0.001)
                end
                # cell 2 is fill on every date, so it must come back NaN with n_valid 0.
                push!(rows_by_date, d => [(a1, a1, 0), (32767, 32767, 255), (a3, a3, q3)])
            end
            _write_cache(root, cells, rows_by_date)

            pct = 0.05
            ice = pool_ice_albedo_from_cache(cells; cache_path = root, percentile = pct,
                                            albedo_range = (0.25, 1.0), progress = false)
            b1 = _brute(truth1, pct)
            b3 = _brute(truth3, pct)

            @test ice[:albedo_bsa][1] === b1[1]
            @test ice[:n_valid_observations_bsa][1] == b1[2]
            @test ice[:k_used_bsa][1] == b1[3]

            @test isnan(ice[:albedo_bsa][2])
            @test ice[:n_valid_observations_bsa][2] == 0
            @test ice[:k_used_bsa][2] == 0

            @test ice[:albedo_bsa][3] === b3[1]
            @test ice[:n_valid_observations_bsa][3] == b3[2]
            @test ice[:k_used_bsa][3] == b3[3]

            # No sample-count threshold: a cell with a single retrieval still reports a value,
            # qualified by n_valid rather than suppressed.
            @test all(>=(0), collect(ice[:n_valid_observations_bsa]))

            @test collect(ice[:cell_id]) == [G._modis_cell_id(c) for c in cells]
            @test collect(ice[:latitude]) ≈ [G._modis_cell_center(c...)[1] for c in cells]

            md = DimensionalData.metadata(ice)
            @test md["n_dates"] == length(dates)
            @test md["percentile"] == pct
            @test md["albedo_range"] == [0.25, 1.0]
            @test md["n_cells"] == 3

            # Pooling is order independent — the accumulator retains a set, and `_finalize`
            # sorts it before averaging.
            shuffled = pool_ice_albedo_from_cache(cells; cache_path = root, percentile = pct,
                                                 albedo_range = (0.25, 1.0),
                                                 dates = reverse(dates), progress = false)
            @test isequal(collect(shuffled[:albedo_bsa]), collect(ice[:albedo_bsa]))

            # The floor is load-bearing, not cosmetic: raising it to the per-year run's 0.30
            # drops real dark retrievals and brightens the answer.
            hi = pool_ice_albedo_from_cache(cells; cache_path = root, percentile = pct,
                                           albedo_range = (0.30, 1.0), progress = false)
            @test hi[:n_valid_observations_bsa][3] < ice[:n_valid_observations_bsa][3]
            @test hi[:albedo_bsa][3] > ice[:albedo_bsa][3]

            # Restricting the dates restricts the statistic, and kmax follows the date count.
            half = pool_ice_albedo_from_cache(cells; cache_path = root, percentile = pct,
                                             albedo_range = (0.25, 1.0),
                                             dates = dates[1:20], progress = false)
            @test half[:n_valid_observations_bsa][1] == 20
            @test half[:albedo_bsa][1] === _brute(truth1[1:20], pct)[1]
        end
    end

    @testset "QA whitelist rejects, and is a whitelist not a reject-list" begin
        mktempdir() do root
            cells = [(17, 2, 100, 100)]
            dates = Date(2001, 1, 1) .+ Day.(0:7)
            # QA classes 0,2,4,6 are full inversions (kept); 1,3,5,7 are magnitude (rejected).
            rows = [dates[i] => [(400, 400, i - 1)] for i in 1:8]
            _write_cache(root, cells, rows)
            ice = pool_ice_albedo_from_cache(cells; cache_path = root, progress = false)
            @test ice[:n_valid_observations_bsa][1] == 4
            # An unknown future class must be rejected by a whitelist, which a reject-list
            # could not do.
            ice9 = pool_ice_albedo_from_cache(cells; cache_path = root, qa_keep = [0],
                                             progress = false)
            @test ice9[:n_valid_observations_bsa][1] == 1
        end
    end

    @testset "the 0.001 scale is ours, and fill is rejected by arithmetic" begin
        mktempdir() do root
            cells = [(17, 2, 100, 100)]
            # 32767 * 0.001 = 32.767, rejected by the range check rather than by luck; 1200
            # scales to 1.2, legitimately above 1.0 but outside the QC range.
            rows = [Date(2001, 1, 1) => [(32767, 32767, 0)],
                    Date(2001, 1, 2) => [(1200, 1200, 0)],
                    Date(2001, 1, 3) => [(400, 400, 0)]]
            _write_cache(root, cells, rows)
            ice = pool_ice_albedo_from_cache(cells; cache_path = root, progress = false)
            @test ice[:n_valid_observations_bsa][1] == 1
            @test ice[:albedo_bsa][1] ≈ 400 * _MCD43A3_SCALE
        end
    end

    @testset "input validation and cache integrity" begin
        mktempdir() do root
            cells = [(17, 2, 100, 100)]
            _write_cache(root, cells, [Date(2001, 1, 1) => [(400, 400, 0)]])

            @test_throws "no cells supplied" pool_ice_albedo_from_cache(
                NTuple{4,Int}[]; cache_path = root, progress = false)
            @test_throws "percentile must be in (0, 1]" pool_ice_albedo_from_cache(
                cells; cache_path = root, percentile = 0.0, progress = false)
            @test_throws "percentile must be in (0, 1]" pool_ice_albedo_from_cache(
                cells; cache_path = root, percentile = 1.5, progress = false)
            @test_throws "albedo_range must be" pool_ice_albedo_from_cache(
                cells; cache_path = root, albedo_range = (1.0, 0.25), progress = false)
            @test_throws "is not an albedo layer" pool_ice_albedo_from_cache(
                cells; cache_path = root, layer = :QFLAG, progress = false)
            # A requested date that is not cached is an error, never a silent gap: a pooled
            # statistic over an unknown subset of the record is not interpretable.
            @test_throws "not cached" pool_ice_albedo_from_cache(
                cells; cache_path = root, dates = [Date(1999, 1, 1)], progress = false)
            @test_throws "no dates to fold" pool_ice_albedo_from_cache(
                cells; cache_path = root, dates = Date[], progress = false)
        end
        # A cell list that does not match the cache resolves to a different key, so the
        # directory is absent rather than half-matching.
        mktempdir() do root
            _write_cache(root, [(17, 2, 100, 100)], [Date(2001, 1, 1) => [(400, 400, 0)]])
            @test_throws "no sample cache at" pool_ice_albedo_from_cache(
                [(17, 2, 100, 999)]; cache_path = root, progress = false)
        end
        @test_throws "no sample cache at" _pooled_cached_dates(mktempdir(), "deadbeefdeadbeef")
    end

    @testset "cached dates are discovered ascending" begin
        mktempdir() do root
            cells = [(17, 2, 100, 100)]
            dates = [Date(2003, 5, 1), Date(2001, 1, 2), Date(2010, 12, 31)]
            key = _write_cache(root, cells, [d => [(400, 400, 0)] for d in dates])
            @test _pooled_cached_dates(root, key) == sort(dates)
        end
    end
end
