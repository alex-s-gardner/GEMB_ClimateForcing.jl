# Contributing to GEMB_ClimateForcing.jl

Thank you for considering contributing to GEMB_ClimateForcing.jl! This document provides guidelines for contributing to the project.

## Getting Started

1. **Fork the repository** and clone your fork
2. **Set up the development environment**:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
3. **Get a CDS API key** (for integration tests):
   - Sign up at https://cds.climate.copernicus.eu/
   - Accept the terms and conditions
   - Get your API key from your profile
   - Set environment variable: `export CDS_API_KEY="your-key"`

## Development Workflow

### Making Changes

1. Create a new branch for your feature/fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes, following the code style guidelines below

3. Add tests for new functionality

4. Run the test suite:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```

### Code Style

- Follow standard Julia style conventions
- Use descriptive variable names (e.g., `temperature_air`, not `ta`)
- Add docstrings for public functions
- Include type annotations for clarity
- Keep functions focused and composable

### Testing

- **Unit tests** should not require external API calls
- **Integration tests** can use the CDS API (will be skipped if no key)
- All tests should be deterministic where possible
- Use small time ranges (1-2 days) for integration tests

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb (e.g., "Add", "Fix", "Update", "Remove")
- Reference issues when applicable: "Fix #123"

## Adding New Datasets

To add support for a new reanalysis dataset:

1. Create `src/datasets/your_dataset.jl`
2. Implement `load_your_dataset(lat, lon; time_range, token, kwargs...)`
3. Return a `GEMB.ClimateForcing` struct with all required variables
4. Call `validate_climate_forcing_units()` before creating the struct
5. Add dispatch in `src/interface.jl`
6. Add tests in `test/runtests.jl`
7. Update README.md and documentation

See `src/datasets/era5_land.jl` as a reference implementation.

## Required Variables

All dataset loaders must provide these variables:
- `temperature_air` (K)
- `pressure_air` (Pa)
- `precipitation` (kg/m² per hour)
- `wind_speed` (m/s)
- `shortwave_downward` (W/m²)
- `longwave_downward` (W/m²)
- `vapor_pressure` (Pa)

**Important**: Call `validate_climate_forcing_units()` to ensure correct units!

## Pull Request Process

1. Ensure all tests pass
2. Update documentation if needed
3. Add an entry to CHANGELOG.md (if present)
4. Submit PR with clear description of changes
5. Address review feedback

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and grow
- Follow the Julia Community standards

## Questions?

- Open an issue for bugs or feature requests
- Tag issues appropriately (bug, enhancement, documentation, etc.)
- Provide minimal reproducible examples for bugs

## License

By contributing, you agree that your contributions will be licensed under the same license as the project.
