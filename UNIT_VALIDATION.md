# Unit Validation System

## Overview

The package now includes automatic unit validation for climate forcing variables before constructing the `ClimateForcing` struct. This catches common unit conversion errors early and provides clear error messages.

## Implementation

**Location:** `src/datasets/era5_land.jl`

The `validate_climate_forcing_units()` function checks that all climate variables have physically reasonable values before passing them to `GEMB.initialize_forcing()`.

## Expected Units and Ranges

| Variable | Units | Valid Range | Common Errors Caught |
|----------|-------|-------------|---------------------|
| `temperature_air` | K | [180, 330] | Temperature in °C instead of K |
| `pressure_air` | Pa | [30000, 110000] | Pressure in hPa instead of Pa |
| `precipitation` | kg/m²/hr | [0, 100] | Precipitation in m instead of kg/m² |
| `wind_speed` | m/s | [0, 100] | Negative values, unrealistic speeds |
| `shortwave_downward` | W/m² | [0, 1500] | Radiation in J/m² instead of W/m² |
| `longwave_downward` | W/m² | [50, 500] | Radiation in J/m² instead of W/m² |
| `vapor_pressure` | Pa | [0, 10000] | Negative values, exceeds saturation |

## Validation Process

The validation runs automatically during data loading:

```julia
cf = climate_forcing(
    :era5land, lat, lon;
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"]
)
# Validation happens internally before returning
```

## Error Messages

When validation fails, you'll see clear messages indicating:
1. Which variable(s) failed validation
2. Expected range
3. Actual range found
4. Possible cause (e.g., "possible unit error: should be W/m², not J/m²")

Example error:
```
ArgumentError: Climate forcing unit validation failed:
  shortwave_downward: expected ≤ 1500 W/m², got 720000.0 W/m² (possible unit error: should be W/m², not J/m²)
  longwave_downward: expected ≥ 50 W/m², got 0.08 W/m² (possible unit error: should be W/m², not J/m²)
```

## Testing

### Run all unit validation tests:
```bash
julia --project=. test/test_unit_validation.jl
```

### Test scenarios covered:
- ✓ Valid data passes without errors
- ✓ Temperature in °C caught (180-330 K range)
- ✓ Pressure in hPa caught (30000-110000 Pa range)
- ✓ Negative precipitation caught
- ✓ Precipitation in meters caught (too large)
- ✓ Negative wind speed caught
- ✓ Unrealistic wind speeds caught (>100 m/s)
- ✓ Radiation in J/m² caught (should be W/m²)
- ✓ Negative shortwave radiation caught
- ✓ Vapor pressure errors caught
- ✓ Boundary conditions handled correctly

## Integration

The validation is integrated into the data loading pipeline:

```
Load ERA5-Land data from Zarr
    ↓
Extract and slice time-series
    ↓
Convert units (m → kg/m², J/m² → W/m², etc.)
    ↓
*** VALIDATE UNITS *** ← New step
    ↓
Create ClimateForcing struct
    ↓
Return to user
```

## Benefits

1. **Early error detection**: Catches unit errors before they propagate into GEMB simulations
2. **Clear diagnostics**: Error messages indicate what went wrong and why
3. **Prevents silent failures**: Invalid data won't produce misleading results
4. **Self-documenting**: Expected ranges serve as documentation for required units
5. **Maintainable**: Adding new datasets requires implementing the same validation

## Adding Validation for New Datasets

When implementing a new dataset loader (e.g., `load_merra2()`), call `validate_climate_forcing_units()` before creating the `ClimateForcing` struct:

```julia
# After unit conversions
validate_climate_forcing_units(
    temperature_air=temperature_air,
    pressure_air=pressure_air,
    precipitation=precipitation,
    wind_speed=wind_speed,
    shortwave_downward=shortwave_downward,
    longwave_downward=longwave_downward,
    vapor_pressure=vapor_pressure
)

# Then create ClimateForcing
cf = GEMB.initialize_forcing(...)
```

This ensures consistency across all dataset loaders.
