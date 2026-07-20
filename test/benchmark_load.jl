#!/usr/bin/env julia
# Benchmark & correctness harness for GEMB_ClimateForcing.jl I/O + allocation work.
#
# Standalone script — NOT included by runtests.jl. Drives the load hot path and the
# CPU-side conversion helpers, and provides a golden-output correctness gate used to
# verify that optimizations do not change numerical results.
#
# Usage:
#   export CDS_API_KEY=...           # or have ~/.cdsapirc with `key:`
#   julia --project=. -t auto test/benchmark_load.jl              # full run + record baseline golden
#   julia --project=. -t auto test/benchmark_load.jl --check      # run + assert vs saved golden
#   julia --project=. -t auto test/benchmark_load.jl --cpu-only   # skip network, CPU micro-bench only
#
# The golden snapshot is written to test/.golden_forcing.jls (gitignored).

using GEMB_ClimateForcing
using GEMB_ClimateForcing: dewpoint_to_vapor_pressure, find_nearest_index
using BenchmarkTools
using Serialization
using DimensionalData
using Dates
using Printf

const ARGS_SET = Set(ARGS)
const DO_CHECK = "--check" in ARGS_SET
const CPU_ONLY = "--cpu-only" in ARGS_SET
const GOLDEN_PATH = joinpath(@__DIR__, ".golden_forcing.jls")

# ---------------------------------------------------------------------------
# Credential resolution: ENV first, then ~/.cdsapirc (mirrors get_cds_api_key)
# ---------------------------------------------------------------------------
function resolve_token()
    haskey(ENV, "CDS_API_KEY") && !isempty(strip(ENV["CDS_API_KEY"])) && return String(strip(ENV["CDS_API_KEY"]))
    rc = joinpath(homedir(), ".cdsapirc")
    if isfile(rc)
        for line in readlines(rc)
            startswith(line, "key:") && return String(strip(split(line, ":", limit=2)[2]))
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
const LAT = 72.58   # Summit, Greenland
const LON = -38.46
const RANGE_7D  = (DateTime(2020, 1, 1), DateTime(2020, 1, 8))
const RANGE_1YR = (DateTime(2020, 1, 1), DateTime(2020, 12, 31))

const VARS = (:temperature_air, :pressure_air, :vapor_pressure, :wind_speed,
              :precipitation, :shortwave_downward, :longwave_downward)
const META_KEYS = ("latitude", "longitude", "temperature_air_mean",
                   "wind_speed_mean", "precipitation_mean")

# Reduce a DimStack to a comparable, serializable summary (values + metadata means).
function summarize(stack)
    d = Dict{Symbol,Any}()
    for v in VARS
        d[v] = collect(parent(stack[v]))
    end
    md = DimensionalData.metadata(stack)
    for k in META_KEYS
        d[Symbol(k)] = md[k]
    end
    return d
end

function compare_golden(summary, golden; rtol=1e-10)
    problems = String[]
    for v in VARS
        a, b = summary[v], golden[v]
        if length(a) != length(b)
            push!(problems, "$v: length $(length(a)) vs golden $(length(b))")
        elseif !isapprox(a, b; rtol=rtol, nans=true)
            push!(problems, "$v: values differ beyond rtol=$rtol (max abs Δ = $(maximum(abs.(a .- b))))")
        end
    end
    for k in META_KEYS
        a, b = summary[Symbol(k)], golden[Symbol(k)]
        if !(a isa Number && isapprox(a, b; rtol=rtol)) && a != b
            push!(problems, "$k: $a vs golden $b")
        end
    end
    return problems
end

# Original (pre-Step-4) conversion chain, kept for A/B comparison.
function convert_chain_old(prec_raw, u10, v10, sw_raw)
    precipitation = Float64.(prec_raw) .* 1000.0
    precipitation[precipitation .< 0] .= 0.0
    wind_speed = hypot.(Float64.(u10), Float64.(v10))
    shortwave = Float64.(sw_raw) ./ 3600.0
    shortwave[shortwave .< 0] .= 0.0
    return precipitation, wind_speed, shortwave
end

# Step 4: fused, native-dtype input, no BitArray masks (mirrors the loader).
function convert_chain(prec_raw, u10, v10, sw_raw)
    precipitation = @. max(Float64(prec_raw) * 1000.0, 0.0)
    wind_speed = @. hypot(Float64(u10), Float64(v10))
    shortwave = @. max(Float64(sw_raw) / 3600.0, 0.0)
    return precipitation, wind_speed, shortwave
end

# ---------------------------------------------------------------------------
# CPU micro-benchmarks (network-independent, deterministic allocation counts)
# ---------------------------------------------------------------------------
function cpu_benchmarks()
    println("\n=== CPU micro-benchmarks (synthetic, network-independent) ===")
    n = 8760  # one year of hourly steps
    # Synthetic native-dtype inputs resembling ERA5-Land raw slices.
    dewpoint = Float32.(range(230.0f0, 273.0f0; length=n))
    prec_raw = Float32.(range(-1.0f-5, 5.0f-4; length=n))
    u10      = Float32.(range(-10.0f0, 10.0f0; length=n))
    v10      = Float32.(range(-8.0f0, 8.0f0; length=n))
    sw_raw   = Float32.(range(-100.0f0, 3.0f6; length=n))
    lats     = collect(range(90.0, -90.0; length=1801))

    println("\n-- dewpoint_to_vapor_pressure (n=$n)")
    @btime dewpoint_to_vapor_pressure($dewpoint)

    println("\n-- find_nearest_index (n=$(length(lats)))")
    @btime find_nearest_index($lats, $(-38.46))

    println("\n-- conversion chain OLD (Float64 copies + BitArray masks, n=$n)")
    @btime convert_chain_old($prec_raw, $u10, $v10, $sw_raw)
    println("\n-- conversion chain NEW (fused, native input, n=$n)")
    @btime convert_chain($prec_raw, $u10, $v10, $sw_raw)
    println()
end

# ---------------------------------------------------------------------------
# Network benchmarks
# ---------------------------------------------------------------------------
function network_benchmarks(token)
    cache = mktempdir()
    println("\n=== Network benchmarks (cache_path=$cache) ===")

    # --- 7-day load: cold (compile + network) then warm cache ---
    println("\n-- 7-day load, COLD (compile + network) --")
    @time s7_cold = climate_forcing(:era5land, LAT, LON;
        time_range=RANGE_7D, token=token, chunk_strategy=:geo, cache_path=cache)

    println("\n-- 7-day load, WARM cache (isolates CPU/alloc from network) --")
    GC.gc()
    a_warm = @allocated (global s7_warm = climate_forcing(:era5land, LAT, LON;
        time_range=RANGE_7D, token=token, chunk_strategy=:geo, cache_path=cache))
    @time climate_forcing(:era5land, LAT, LON;
        time_range=RANGE_7D, token=token, chunk_strategy=:geo, cache_path=cache)
    @printf("   @allocated (warm 7-day full load): %.2f MiB\n", a_warm / 2^20)

    # --- 1-year load: warm-cache wall time + allocations ---
    println("\n-- 1-year load, COLD --")
    @time climate_forcing(:era5land, LAT, LON;
        time_range=RANGE_1YR, token=token, chunk_strategy=:geo, cache_path=cache)
    println("\n-- 1-year load, WARM cache --")
    GC.gc()
    a1 = @allocated climate_forcing(:era5land, LAT, LON;
        time_range=RANGE_1YR, token=token, chunk_strategy=:geo, cache_path=cache)
    t1 = @elapsed (global s1yr = climate_forcing(:era5land, LAT, LON;
        time_range=RANGE_1YR, token=token, chunk_strategy=:geo, cache_path=cache))
    @printf("   1-year warm: %.3f s, %.2f MiB allocated\n", t1, a1 / 2^20)

    return s7_warm, s1yr
end

# ---------------------------------------------------------------------------
# Golden gate
# ---------------------------------------------------------------------------
function golden_gate(stack)
    summary = summarize(stack)
    if DO_CHECK
        if !isfile(GOLDEN_PATH)
            error("--check requested but no golden snapshot at $GOLDEN_PATH; run without --check first.")
        end
        golden = deserialize(GOLDEN_PATH)
        problems = compare_golden(summary, golden)
        if isempty(problems)
            println("\n✅ CORRECTNESS GATE PASSED — output matches golden snapshot.")
        else
            println("\n❌ CORRECTNESS GATE FAILED:")
            for p in problems
                println("   - ", p)
            end
            error("Output diverged from golden snapshot.")
        end
    else
        serialize(GOLDEN_PATH, summary)
        println("\n💾 Golden snapshot written to $GOLDEN_PATH (baseline for future --check runs).")
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    println("GEMB_ClimateForcing.jl benchmark harness")
    println("threads = $(Threads.nthreads()), check=$DO_CHECK, cpu_only=$CPU_ONLY")

    cpu_benchmarks()

    if CPU_ONLY
        println("\n--cpu-only: skipping network benchmarks and golden gate.")
        return
    end

    token = resolve_token()
    if isnothing(token)
        println("\n⚠️  No CDS credentials (ENV[CDS_API_KEY] or ~/.cdsapirc). Skipping network path.")
        return
    end

    _, s1yr = network_benchmarks(token)
    golden_gate(s1yr)
end

main()
