# Test authenticated HTTP access (without needing CDS token)
#
# This demonstrates the AuthenticatedHTTPStore working with a test endpoint

using GEMB_ClimateForcing
using HTTP

println("Testing AuthenticatedHTTPStore")
println("="^60)
println()

# Create store with test token
println("1. Creating authenticated store...")
store = GEMB_ClimateForcing.AuthenticatedHTTPStore(
    "https://httpbin.org";
    token="test-token-12345"
)

println("   ✓ Store created")
println("   URL: ", store.url)
println("   Headers: ", store.headers)
println()

# Test authentication by accessing a protected endpoint
println("2. Testing HTTP GET with Bearer token...")
try
    # httpbin.org/bearer expects Authorization: Bearer <token>
    # It returns the token if authenticated correctly
    response_data = store["bearer"]

    if !isnothing(response_data)
        println("   ✓ Authenticated request successful!")
        println("   Response: ", String(response_data)[1:min(100, length(response_data))], "...")
    else
        println("   ⚠ Received nothing (404 or similar)")
    end
catch e
    println("   Error: ", e)
end
println()

# Show that this is exactly how Zarr.jl's GCStore works
println("3. Pattern comparison:")
println("   This implementation follows Zarr.jl's GCStore pattern:")
println()
println("   GCStore (Google Cloud):")
println("     headers[\"Authorization\"] = \"Bearer \$(access_token)\"")
println("     HTTP.request(\"GET\", url, headers, ...)")
println()
println("   AuthenticatedHTTPStore (ERA5-Land):")
println("     headers[\"Authorization\"] = \"Bearer \$(token)\"")
println("     HTTP.request(\"GET\", url, headers, ...)")
println()
println("   ✓ Same pattern, different cloud provider!")
println()

println("="^60)
println("AuthenticatedHTTPStore is ready for ERA5-Land!")
println("="^60)
println()
println("Next steps:")
println("  1. Get CDS API key: https://cds.climate.copernicus.eu/")
println("  2. export CDS_API_KEY=\"your-key-here\"")
println("  3. Run: julia --project=. examples/era5_land_example.jl")
