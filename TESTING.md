# Testing GEMB_ClimateForcing.jl

## Quick Start

### Run All Tests (without data retrieval)
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This runs validation tests that don't require CDS API credentials.

### Run Data Retrieval Tests

**Prerequisites:**
1. Get a free CDS API key from https://cds.climate.copernicus.eu/
2. Set the environment variable:
   ```bash
   export CDS_API_KEY="your-token-here"
   ```

**Then run:**
```bash
# Full test suite (includes data retrieval)
julia --project=. -e 'using Pkg; Pkg.test()'

# Or standalone data retrieval test
julia --project=. test/test_data_retrieval.jl
```

## Test Files

### `test/runtests.jl` (Main Test Suite)
Comprehensive test suite covering:
- ✅ Package loading
- ✅ Input validation (lat/lon/time_range/chunk_strategy)
- ✅ Error handling
- ⚠️ ERA5-Land data retrieval (requires CDS_API_KEY)
- ⚠️ GEMB integration (requires CDS_API_KEY)

### `test/test_api_validation.jl`
API validation tests that run WITHOUT credentials:
- Required arguments
- Latitude/longitude ranges
- Time range validation
- Chunk strategy validation
- Dataset support

```bash
julia --project=. test/test_api_validation.jl
```

### `test/test_data_retrieval.jl`
Dedicated data retrieval tests (requires CDS_API_KEY):
1. Minimal 2-day request (Summit, Greenland)
2. Different location (Antarctica)
3. Time-chunked strategy
4. Full GEMB integration test

```bash
export CDS_API_KEY="your-token-here"
julia --project=. test/test_data_retrieval.jl
```

## Test Coverage

### Without CDS API Key ✓
- Package loading
- Input validation
- Error messages
- API structure

### With CDS API Key ✓
- Actual data retrieval from ERA5-Land ARCO
- Physical range validation
- GEMB integration
- Different locations and time ranges
- Both chunk strategies (:geo and :time)

## Expected Test Duration

| Test Suite | Duration | Network Required |
|------------|----------|------------------|
| API Validation | ~5 seconds | No |
| Data Retrieval (2 days) | ~30-60 seconds | Yes |
| Full Integration | ~60-90 seconds | Yes |

## Continuous Integration

For CI environments without CDS credentials:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This will run validation tests and skip data retrieval tests gracefully.

For CI with CDS credentials:
```yaml
env:
  CDS_API_KEY: ${{ secrets.CDS_API_KEY }}
run: |
  julia --project=. -e 'using Pkg; Pkg.test()'
```

## Troubleshooting

### "CDS_API_KEY not set"
- Get a key from https://cds.climate.copernicus.eu/
- Set: `export CDS_API_KEY="your-token"`

### HTTP 401/403 errors
- Check that your CDS API key is valid
- Verify you can log in to the CDS website

### "HTTP timeout" or slow downloads
- First download from a new store takes 30-60 seconds (normal)
- Check your network connection
- ECMWF servers may be under load

### Tests fail with GEMB errors
- Ensure GEMB.jl is installed: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Check GEMB.jl compatibility version in Project.toml
