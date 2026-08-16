"""
Tests for `compute_glacier_ice_albedo` — glacier bare-ice albedo as the mean of the
darkest few percent of a year's C3S albedo observations.

Everything here is offline. The reduction kernel (`_low_percentile_mean`), the
per-observation QC mask (`_valid_albedo`), the QFLAG legend parsing (`_qflag_rejects`)
and the driver's input validation need no CDS access; the CDS-ordering path they sit on
top of is covered by `test_copernicus_albedo.jl`, which gates its live orders behind
`GEMB_TEST_SATELLITE_ALBEDO=1`.
"""

using Test
using GEMB_ClimateForcing
using Statistics
using Dates
using Rasters
using NCDatasets

using GEMB_ClimateForcing: _low_percentile_mean, _LowPercentileTopK, _accumulate!,
    _finalize, _valid_albedo, _qflag_rejects, _qflag_table, _albedo_qflag_legend,
    _albedo_isbad, _ICE_ALBEDO_RANGE, GLACIER_ICE_ALBEDO_QFLAG_REJECT

@testset "Glacier bare-ice albedo" begin

    @testset "Low-percentile mean" begin
        # 20 observations 0.41 … 0.60. At 5 % that is ceil(0.05*20) = 1 observation:
        # the darkest, 0.41.
        obs = [fill(Float32(0.40 + 0.01i), 2, 3) for i in 1:20]
        stat, counts = _low_percentile_mean(obs, 0.05, 10)
        @test all(≈(0.41f0), stat)
        @test all(==(Int32(20)), counts)

        # NaN observations are skipped, and the count reflects that. With 19 valid,
        # ceil(0.05*19) = 1 → still the darkest survivor, now 0.42.
        obs[1][1, 1] = NaN32
        stat, counts = _low_percentile_mean(obs, 0.05, 10)
        @test stat[1, 1] ≈ 0.42f0
        @test counts[1, 1] == 19
        @test stat[2, 2] ≈ 0.41f0
        @test counts[2, 2] == 20

        # percentile = 1 averages everything, i.e. the plain annual mean.
        stat, _ = _low_percentile_mean(obs, 1.0, 1)
        @test stat[2, 2] ≈ Float32(mean(0.41:0.01:0.60))

        # At least one observation is always averaged, even when percentile*n rounds to 0.
        stat, _ = _low_percentile_mean(obs, 0.001, 1)
        @test stat[2, 2] ≈ 0.41f0

        # 25 % of 20 = 5 darkest: 0.41 … 0.45.
        stat, _ = _low_percentile_mean(obs, 0.25, 10)
        @test stat[2, 2] ≈ Float32(mean(0.41:0.01:0.45))

        # Under-sampled pixel-years are NaN, but their count is still reported so a
        # caller can see *why* the pixel is missing.
        stat, counts = _low_percentile_mean(obs, 0.05, 25)
        @test all(isnan, stat)
        @test all(>(0), counts)

        # An all-NaN pixel yields NaN with a zero count, never an error.
        blank = [fill(NaN32, 2, 2) for _ in 1:12]
        stat, counts = _low_percentile_mean(blank, 0.05, 1)
        @test all(isnan, stat)
        @test all(==(Int32(0)), counts)

        # Argument validation.
        @test_throws ArgumentError _low_percentile_mean(Matrix{Float32}[], 0.05, 10)
        @test_throws ArgumentError _low_percentile_mean(obs, 0.0, 10)
        @test_throws ArgumentError _low_percentile_mean(obs, 1.5, 10)
        @test_throws ArgumentError _low_percentile_mean(
            [fill(0.5f0, 2, 2), fill(0.5f0, 3, 3)], 0.05, 1)
    end

    @testset "Streaming reduction matches the batch form" begin
        # The driver streams timesteps through the accumulator instead of holding the
        # year's grids; it must agree bit-for-bit with sorting every value, since it keeps
        # the same k smallest values and sorts them the same way before averaging.
        obs = [Float32[0.4+0.013i 0.9-0.02i; 0.55 0.3+0.01i] for i in 1:36]
        obs[3][1, 1] = NaN32
        obs[7][2, 2] = NaN32
        for p in (0.05, 0.1, 0.25, 1.0)
            acc = _LowPercentileTopK(size(first(obs)), length(obs), p, 5)
            for m in obs
                _accumulate!(acc, m)
            end
            stat, counts = _finalize(acc)
            ref, ref_counts = _low_percentile_mean(obs, p, 5)
            @test all(((a, b),) -> (isnan(a) && isnan(b)) || a === b, zip(stat, ref))
            @test counts == ref_counts
        end

        # Feeding more grids than declared is a programming error, not silent truncation.
        acc = _LowPercentileTopK((2, 2), 1, 0.5, 1)
        _accumulate!(acc, fill(0.5f0, 2, 2))
        @test_throws ArgumentError _accumulate!(acc, fill(0.5f0, 2, 2))
        @test_throws ArgumentError _accumulate!(
            _LowPercentileTopK((2, 2), 2, 0.5, 1), fill(0.5f0, 3, 3))
        @test_throws ArgumentError _LowPercentileTopK((2, 2), 0, 0.5, 1)
    end

    @testset "Per-observation QC" begin
        lo, hi = _ICE_ALBEDO_RANGE

        # Missing / non-finite retrievals are rejected however NCDatasets surfaces them.
        @test _albedo_isbad(missing)
        @test _albedo_isbad(NaN)
        @test !_albedo_isbad(0.5)
        m = _valid_albedo([0.5 missing; NaN 0.7], nothing, nothing)
        @test m[1, 1] ≈ 0.5f0
        @test isnan(m[1, 2]) && isnan(m[2, 1])
        @test m[2, 2] ≈ 0.7f0
        @test eltype(m) === Float32

        # Out-of-range albedo is rejected at both ends of `albedo_range`.
        m = _valid_albedo([lo - 0.05 lo; hi hi + 0.05], nothing, nothing)
        @test isnan(m[1, 1]) && isnan(m[2, 2])
        @test m[1, 2] ≈ Float32(lo) && m[2, 1] ≈ Float32(hi)

        # A lower floor admits darker ice — the 0.3 default is glaciological, not physical.
        m = _valid_albedo([0.2 0.5], nothing, nothing; albedo_range=(0.15, 1.0))
        @test m[1] ≈ 0.2f0

        # QFLAG: bitwise rejection, with unflagged and unrelated bits kept.
        m = _valid_albedo([0.5 0.5 0.5], [0 2 5], nothing; bitmask=2)
        @test m[1] ≈ 0.5f0        # no bits set
        @test isnan(m[2])         # bit 1 set → rejected
        @test m[3] ≈ 0.5f0        # bits 0 and 2 set, neither in the mask

        # QFLAG: exact-value rejection (flag_values products), and a missing flag.
        m = _valid_albedo([0.5 0.5 0.5], [1 7 missing], nothing; values=[7])
        @test m[1] ≈ 0.5f0
        @test isnan(m[2])
        @test isnan(m[3])

        # Uncertainty: oversized rejected, missing tolerated, threshold off entirely.
        m = _valid_albedo([0.5 0.5 0.5], nothing, [0.1 0.9 missing]; max_error=0.2)
        @test m[1] ≈ 0.5f0
        @test isnan(m[2])
        @test m[3] ≈ 0.5f0
        m = _valid_albedo([0.5 0.5], nothing, [0.1 0.9]; max_error=nothing)
        @test all(≈(0.5f0), m)
    end

    @testset "QFLAG legend → reject mask" begin
        # The real v3.1 legend shape: snow_presence, then per-band no-obs and BRDF bits.
        table = [(meaning="snow_presence", mask=1, is_value=false),
                 (meaning="no_obs_in_last_decade_for_Oa03", mask=2, is_value=false),
                 (meaning="no_obs_in_last_decade_for_Oa04", mask=4, is_value=false),
                 (meaning="brdf_warning_for_Oa03", mask=1024, is_value=false)]

        bitmask, values, matched = _qflag_rejects(table)
        @test bitmask == 2 | 4
        @test isempty(values)
        @test length(matched) == 2
        # snow_presence is the signal, not an error — it must never be rejected.
        @test !any(contains("snow"), matched)
        @test bitmask & 1 == 0
        # brdf_warning is off by default, opt-in only.
        @test bitmask & 1024 == 0
        bitmask_brdf, _, matched_brdf = _qflag_rejects(table; brdf_warning=true)
        @test bitmask_brdf == 2 | 4 | 1024
        @test length(matched_brdf) == 3

        # flag_values products come back as exact values rather than a mask.
        vtable = [(meaning="cloud", mask=3, is_value=true),
                  (meaning="clear", mask=0, is_value=true)]
        bitmask, values, matched = _qflag_rejects(vtable)
        @test bitmask == 0
        @test values == [3]
        @test matched == ["cloud"]

        # An unreadable/absent legend degrades to no rejection, not an error.
        @test _qflag_rejects(NamedTuple[]) == (0, Int[], String[])

        # The default pattern list must not accidentally match snow_presence.
        @test !any(p -> occursin(p, "snow_presence"), GLACIER_ICE_ALBEDO_QFLAG_REJECT)
    end

    @testset "QFLAG legend read from a file" begin
        # The legend is parsed from the file's own CF attributes, so exercise that against
        # a real NetCDF laid out like a delivered product file rather than a hand-built
        # table. The driver reaches the file via the `"file"` key `satellite_albedo` puts
        # in each layer's metadata, so this covers the path it actually takes.
        dir = mktempdir()
        try
            withlegend = joinpath(dir, "withlegend.nc")
            NCDataset(withlegend, "c") do ds
                defDim(ds, "lon", 3)
                defDim(ds, "lat", 2)
                q = defVar(ds, "QFLAG", Int32, ("lon", "lat"))
                q[:, :] = Int32[0 1; 2 4; 8 0]
                q.attrib["flag_masks"] = Int32[1, 2, 4]
                q.attrib["flag_meanings"] =
                    "snow_presence no_obs_in_last_decade_for_Oa03 brdf_warning_for_Oa03"
                a = defVar(ds, "AL_DH_BB", Float32, ("lon", "lat"))
                a[:, :] = fill(0.5f0, 3, 2)
            end

            table = _qflag_table(withlegend)
            @test length(table) == 3
            @test table[1].meaning == "snow_presence"
            @test all(t -> !t.is_value, table)      # flag_masks, not flag_values
            bitmask, values, matched = _qflag_rejects(table)
            @test bitmask == 2
            @test isempty(values)
            @test matched == ["no_obs_in_last_decade_for_Oa03"]

            # A file with no QFLAG variable degrades to an empty legend, not an error.
            nolegend = joinpath(dir, "nolegend.nc")
            NCDataset(nolegend, "c") do ds
                defDim(ds, "lon", 2)
                defDim(ds, "lat", 2)
                defVar(ds, "AL_DH_BB", Float32, ("lon", "lat"))[:, :] = fill(0.5f0, 2, 2)
            end
            @test isempty(_qflag_table(nolegend))

            # End to end from a file path, as the driver calls it — including the degraded
            # paths (no legend in the file, no file at all), which warn instead of erroring.
            @test _albedo_qflag_legend(withlegend) == (2, Int[],
                                                       ["no_obs_in_last_decade_for_Oa03"])
            @test _albedo_qflag_legend(withlegend; brdf_warning=true)[1] == 2 | 4
            # A caller-supplied pattern list overrides the default (it used to be inert).
            @test _albedo_qflag_legend(withlegend; patterns=["snow"])[1] == 1
            @test (@test_logs (:warn,) _albedo_qflag_legend(nolegend)) ==
                  (0, Int[], String[])
            @test (@test_logs (:warn,) _albedo_qflag_legend(nothing)) ==
                  (0, Int[], String[])
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "Input validation" begin
        # These all throw before any network access, so they are safe without a token.
        @test_throws ArgumentError compute_glacier_ice_albedo(Int[])
        @test_throws ArgumentError compute_glacier_ice_albedo(2019; percentile=0.0)
        @test_throws ArgumentError compute_glacier_ice_albedo(2019; percentile=1.5)
        @test_throws ArgumentError compute_glacier_ice_albedo(2019; min_samples=0)
        @test_throws ArgumentError compute_glacier_ice_albedo(2019;
                                                             albedo_range=(1.0, 0.3))

        # A single year is accepted as a scalar as well as a range.
        @test hasmethod(compute_glacier_ice_albedo, Tuple{Int})
        @test GLACIER_ICE_ALBEDO_YEARS == 2018:2024
    end
end
