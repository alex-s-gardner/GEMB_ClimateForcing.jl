# Test caching performance
#
# This script demonstrates the performance improvement from chunk caching.
# Requires CDS_API_KEY to be set.

using GEMB_ClimateForcing
using Dates

# Get CDS API key
token = get(ENV, "CDS_API_KEY", nothing)
if isnothing(token) || isempty(strip(token))
    println("⚠️  CDS_API_KEY not set. Skipping caching test.")
    println("Set export CDS_API_KEY='your-token-here' to run this test.")
    exit(0)
end

# Test location and time range
lat, lon = 72.58, -38.46  # Summit, Greenland
time_range = (DateTime(2020, 1, 1), DateTime(2020, 1, 2))  # 1 day

println("="^70)
println("Testing Chunk Caching Performance")
println("="^70)

# First load (cold cache)
println("\n1. First load (cold cache)...")
t1_start = time()
data1 = climate_forcing(:era5land, lat, lon; time_range=time_range, token=token, cache_size=128)
t1_end = time()
t1_elapsed = t1_end - t1_start
println("   Time: $(round(t1_elapsed, digits=2))s")

# Second load (warm cache - should be much faster)
println("\n2. Second load (warm cache)...")
t2_start = time()
data2 = climate_forcing(:era5land, lat, lon; time_range=time_range, token=token, cache_size=128)
t2_end = time()
t2_elapsed = t2_end - t2_start
println("   Time: $(round(t2_elapsed, digits=2))s")

# Verify data is identical
println("\n3. Verifying data consistency...")
@assert data1[:temperature_air] == data2[:temperature_air]
@assert data1[:pressure_air] == data2[:pressure_air]
@assert data1[:wind_speed] == data2[:wind_speed]
println("   ✓ Data identical between loads")

# Calculate speedup
speedup = t1_elapsed / t2_elapsed
println("\n" * "="^70)
println("Results:")
println("  First load:  $(round(t1_elapsed, digits=2))s")
println("  Second load: $(round(t2_elapsed, digits=2))s")
println("  Speedup:     $(round(speedup, digits=1))x")
println("="^70)

if speedup > 2.0
    println("\n✓ Caching provides significant speedup!")
else
    println("\n⚠️  Speedup lower than expected ($(round(speedup, digits=1))x)")
    println("    This may be due to HTTP connection caching masking chunk cache benefits.")
end
