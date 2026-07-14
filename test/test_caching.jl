# Test disk caching performance using Zarr.CachingStore
#
# Demonstrates persistent caching: first load fetches from network,
# second load reads from local disk cache.
#
# Requires CDS_API_KEY to be set.

using GEMB_ClimateForcing
using Dates

# Get CDS API key
token = get(ENV, "CDS_API_KEY", nothing)
if isnothing(token) || isempty(strip(token))
    println("CDS_API_KEY not set. Skipping caching test.")
    exit(0)
end

# Test location and time range
lat, lon = 72.58, -38.46  # Summit, Greenland
time_range = (DateTime(2020, 1, 1), DateTime(2020, 1, 2))  # 1 day

# Use a temporary directory for the cache
cache_path = mktempdir()

println("="^70)
println("Testing Disk Caching (Zarr.CachingStore)")
println("="^70)
println("  Cache path: $cache_path")

# First load (cold cache - fetches from network)
println("\n1. First load (cold cache - network fetch)...")
t1_start = time()
data1 = climate_forcing(:era5land, lat, lon; time_range=time_range, token=token, cache_path=cache_path)
t1_elapsed = time() - t1_start
println("   Time: $(round(t1_elapsed, digits=2))s")

# Second load (warm cache - reads from disk)
println("\n2. Second load (warm cache - disk read)...")
t2_start = time()
data2 = climate_forcing(:era5land, lat, lon; time_range=time_range, token=token, cache_path=cache_path)
t2_elapsed = time() - t2_start
println("   Time: $(round(t2_elapsed, digits=2))s")

# Verify data is identical
println("\n3. Verifying data consistency...")
@assert data1[:temperature_air] == data2[:temperature_air]
@assert data1[:pressure_air] == data2[:pressure_air]
@assert data1[:wind_speed] == data2[:wind_speed]
println("   Data identical between loads")

# Calculate speedup
speedup = t1_elapsed / t2_elapsed
println("\n" * "="^70)
println("Results:")
println("  First load (network): $(round(t1_elapsed, digits=2))s")
println("  Second load (cache):  $(round(t2_elapsed, digits=2))s")
println("  Speedup:              $(round(speedup, digits=1))x")
println("="^70)
