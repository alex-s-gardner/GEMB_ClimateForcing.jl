# CI/CD Setup Guide

## Current Status

The GitHub Actions CI is configured but **integration tests are skipped** because the `CDS_API_KEY` secret is not configured in the repository.

## What's Working

✅ **Basic tests** (27 tests):
- Package loading
- Unit validation (physical ranges, units)
- Input validation (lat/lon ranges, time ranges, arguments)
- Error handling (missing/invalid tokens)

✅ **CI runs on**:
- Julia 1.10 (LTS) and Julia 1.x (latest)
- Ubuntu, macOS, and Windows
- Both `main` and `feature/**` branches

## What's Skipped

⏭️ **Integration tests** (2 tests):
- Actual ERA5-Land data loading from ECMWF servers
- Parallel loading consistency checks

These tests require a valid CDS API key to access the ECMWF ARCO Zarr stores.

## Enabling Integration Tests in CI

To enable full integration testing in GitHub Actions:

### 1. Get a CDS API Key

1. Create a free account at https://cds.climate.copernicus.eu/
2. Navigate to your user profile
3. Copy your API key (format: `UID:TOKEN`)

### 2. Add Secret to GitHub Repository

1. Go to your GitHub repository
2. Navigate to **Settings** → **Secrets and variables** → **Actions**
3. Click **New repository secret**
4. Name: `CDS_API_KEY`
5. Value: Paste your API key from step 1
6. Click **Add secret**

### 3. Verify

Once the secret is added:
- Push a new commit or re-run the workflow
- Integration tests will automatically run in CI
- Check the logs for "Running ERA5-Land Data Retrieval Tests"

## Test Behavior

The test suite uses conditional testing:

```julia
# Only run if CDS_API_KEY is available and non-empty
token = get(ENV, "CDS_API_KEY", nothing)

if !isnothing(token) && !isempty(strip(token))
    # Run integration tests
else
    # Skip with informative message
end
```

This ensures:
- ✅ CI doesn't fail when secret is not configured
- ✅ Basic tests always run
- ✅ Integration tests run when credentials are available
- ✅ No accidental API calls with empty/invalid tokens

## Local Testing

To run integration tests locally:

```bash
# Set your API key
export CDS_API_KEY="your-key-here"

# Run all tests (including integration)
julia --project=. -e 'using Pkg; Pkg.test()'

# Or run the example directly
julia --project=. examples/era5_land_example.jl
```

Without `CDS_API_KEY` set, basic tests still run but integration tests are skipped.

## Rate Limiting

The ECMWF ARCO stores have usage limits. To avoid excessive API calls:
- Integration tests use minimal time ranges (1-2 days)
- Tests run in parallel to minimize total requests
- Consider running integration tests only on the `main` branch in CI

To restrict integration tests to specific branches, modify `.github/workflows/CI.yml`:

```yaml
- name: Run tests
  uses: julia-actions/julia-runtest@v1
  env:
    # Only enable integration tests on main branch
    CDS_API_KEY: ${{ github.ref == 'refs/heads/main' && secrets.CDS_API_KEY || '' }}
```

## Troubleshooting

### Tests Skip with "CDS_API_KEY not set"
- Check that the secret is added in GitHub repository settings
- Verify the secret name is exactly `CDS_API_KEY`
- Check workflow logs to see if the environment variable is available

### HTTP 401 "Invalid authorization header format"
- Check that your API key format is correct: `UID:TOKEN`
- Verify there are no extra spaces or newlines in the secret
- Regenerate the key from CDS if needed

### Tests Fail Locally But Pass in CI
- Make sure you've run `julia --project=. -e 'using Pkg; Pkg.instantiate()'` first
- Check that you're using Julia 1.10+ (`julia --version`)
- Verify OpenSSL is available on your system

## Current CI Results

With the recent fixes:
- ✅ Resolved Statistics version conflict (stdlib version tied to Julia version)
- ✅ Integration tests properly skip when credentials unavailable
- ✅ All basic tests pass on Julia 1.10+ across all platforms
