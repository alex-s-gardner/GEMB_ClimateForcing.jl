"""
GEMB_ClimateForcing

Load climate forcing data from various reanalysis datasets and return DimStack
with climate variables. Supports ERA5-Land, ERA5, MERRA-2, and other datasets.

# Example
```julia
using GEMB_ClimateForcing
using GEMB  # Extension provides DimStack → ClimateForcing conversion

# Load ERA5-Land data for Summit, Greenland
forcing_data = climate_forcing(
    :era5land, 72.58, -38.46;
    time_range=(DateTime(2020,1,1), DateTime(2020,12,31)),
    token=ENV["CDS_API_KEY"]
)

# Convert to GEMB.ClimateForcing (requires GEMB.jl)
cf = GEMB.ClimateForcing(forcing_data)

# Use with GEMB
output = gemb(profile, cf, mp)
```
"""
module GEMB_ClimateForcing

using Dates
using DimensionalData
using PrecompileTools
using Random
using Statistics
using Zarr
using HTTP
using OpenSSL
using CodecZlib

# Export main interface
export climate_forcing, climate_chunk_map, get_cds_api_key, climate_adjust_for_elevation,
    temperature_adjust, precipitation_adjust,
    empirical_lapse_rate, climate_model_invariant, ERA5_LAND_INVARIANT_PARAMETERS,
    GREENLAND_LAPSE_RATE, ARCTIC_LAPSE_RATE, ANTARCTICA_LAPSE_RATE

# Export synthetic forcing and climate fitting
export simulate_climate_forcing, simulation_parameter_sets,
    datetime2decyear, decyear2datetime,
    dewpoint_to_vapor_pressure, vapor_pressure_to_relative_humidity,
    relative_humidity_to_vapor_pressure,
    fit_air_temperature, fit_precipitation, fit_longwave_irradiance_delta,
    fit_seasonal_daily_noise, varname2longname, simulate_coeffs_disp,
    geopotential2height, geoid_undulation, GEOID_MODELS, surface_elevation

# Export satellite observations (ordered from the CDS Retrieve API, not a lazy store)
export satellite_albedo, satellite_albedo_layers, SATELLITE_ALBEDO_VARIABLES

# Export glacier bare-ice albedo derived from the satellite albedo record
export compute_glacier_ice_albedo, GLACIER_ICE_ALBEDO_YEARS,
    GLACIER_ICE_ALBEDO_QFLAG_REJECT

# Export on-glacier temperature decoupling (Shaw et al. 2025 lookup table)
export glacier_decoupling, glacier_decoupling_table, GlacierDecoupling,
    climate_adjust_for_glacier

# Include submodules
include("utils.jl")
include("authenticated_http_store.jl")
# CDS Retrieve (job-based) client — needs get_cds_api_key from utils.jl.
include("cds_retrieve.jl")
include("interface.jl")
include("datasets/era5_land.jl")
include("datasets/copernicus_dem.jl")
include("datasets/copernicus_albedo.jl")
# Glacier bare-ice albedo (darkest-percentile annual reduction of the albedo record).
# Must follow copernicus_albedo.jl — reuses satellite_albedo and its cache-layout helpers.
include("glacier_ice_albedo.jl")
include("geoid.jl")
include("invariant.jl")
include("elevation_adjustment.jl")
# Direct temperature / precipitation perturbations — reuses the saturation-vapor
# and emissivity helpers from elevation_adjustment.jl.
include("climate_adjustment.jl")
# Shaw et al. (2025) per-glacier decoupling factors, read from the vendored
# data/shaw2025_glacier_decoupling.csv.gz lookup table.
include("glacier_decoupling.jl")
# Ambient → on-glacier correction, applying the decoupling factor to a forcing stack.
# Shares _rebuild_forcing / _adjust_longwave (utils.jl) with the two siblings above.
include("glacier_adjustment.jl")

# Synthetic forcing generation
include("simulate/simulate_climate_forcing.jl")

# Climate fitting functions
include("fit_climate/fit_air_temperature.jl")
include("fit_climate/fit_precipitation.jl")
include("fit_climate/fit_longwave_irradiance_delta.jl")
include("fit_climate/fit_seasonal_daily_noise.jl")
include("fit_climate/varname2longname.jl")
include("fit_climate/simulate_coeffs_disp.jl")

# ----------------------------------------------------------------------------
# Precompilation workload
#
# The offline workflows are generic over `DimStack` type parameters that only
# materialize at the first call, so a bare `using` left several seconds of
# inference to the user's first adjustment or fit. Running a small forcing stack
# through the adjustment and fitting entry points here caches those
# specializations into the package image.
#
# Deliberately network-free: nothing here touches a Zarr store, the CDS Retrieve
# API, or a `/vsicurl/` read, so precompilation stays offline and CI-safe. The
# load path (`climate_forcing`) is therefore not covered.
#
# The stack is built on the same layer names, element type and `Ti` lookup as a
# real forcing stack, since those are what the specializations key on. Keep it
# short (a few days of hourly steps) — the workload's own runtime is added to
# every precompile of this package.
# ----------------------------------------------------------------------------
@setup_workload begin
    _times = DateTime(2020, 1, 1):Hour(1):DateTime(2020, 1, 4)
    _n = length(_times)
    _ti = Ti(collect(_times))
    _layer(v) = DimArray(fill(v, _n), (_ti,); metadata=Dict("units" => ""))
    _stack = DimStack((
            temperature_air    = _layer(255.0),
            pressure_air       = _layer(68000.0),
            vapor_pressure     = _layer(120.0),
            wind_speed         = _layer(6.0),
            precipitation      = _layer(0.05),
            shortwave_downward = _layer(90.0),
            longwave_downward  = _layer(190.0),
        ); metadata = Dict{String,Any}(
            "latitude" => 72.58, "longitude" => -38.46, "elevation" => 3200.0,
            "temperature_air_mean" => 255.0, "wind_speed_mean" => 6.0,
            "precipitation_mean" => 438.0,
            "temperature_observation_height" => 2.0,
            "wind_observation_height" => 10.0,
        ))

    @compile_workload begin
        # Adjustment workflows (elevation, direct perturbation, glacier decoupling).
        climate_adjust_for_elevation(_stack, 250.0)
        climate_adjust_for_elevation(_stack, -100.0;
            lapse_rate=GREENLAND_LAPSE_RATE, precip_scaling_method=:clausius_clapeyron)
        temperature_adjust(_stack, 2.0)
        precipitation_adjust(_stack, 1.1)
        climate_adjust_for_glacier(_stack, 0.8)

        # Time-axis conversions and the fitting entry points.
        _dec_year = datetime2decyear(collect(_times))
        decyear2datetime(_dec_year)
        _series = collect(parent(_stack[:wind_speed])) .+ range(0.0, 1.0, _n)
        fit_seasonal_daily_noise(_dec_year, _series)
        fit_air_temperature(_dec_year,
            collect(parent(_stack[:temperature_air])) .+ range(0.0, 2.0, _n), 72.58, 3200.0)
        fit_precipitation(_dec_year, collect(parent(_stack[:precipitation])))
        fit_longwave_irradiance_delta(collect(range(-20.0, 20.0, _n)))

        # Synthetic-forcing kernels. Driven directly on the short series rather than
        # through `simulate_climate_forcing`, which always generates its parameter
        # set's full 32-year hourly record and would add ~70 ms to every precompile.
        _lp, _co = simulation_parameter_sets("test_1")
        _rng = Random.MersenneTwister(_lp.rand_seed)
        simulate_shortwave_irradiance(_dec_year, _lp.latitude)
        simulate_air_temperature(_rng, _dec_year, _lp.latitude, _lp.elevation;
            mean_offset=_co.temperature_air.mean_offset,
            lat_scale=_co.temperature_air.lat_scale,
            daily_amp_scale=_co.temperature_air.daily_amp_scale,
            weather_sigma_scale=_co.temperature_air.weather_sigma_scale,
            weather_corr=_co.temperature_air.weather_corr)
        simulate_air_pressure(_rng, _dec_year, fill(255.0, _n), _lp.latitude, _lp.elevation)
        simulate_seasonal_daily_noise(_rng, _dec_year, _co.wind_speed)
        simulate_precipitation(_rng, _dec_year, _co.precipitation)
        simulate_longwave_irradiance_delta(_rng, _dec_year, _co.longwave_downward)
        simulate_longwave_irradiance(fill(255.0, _n), fill(120.0, _n))

        # Shared unit conversions.
        vapor_pressure_to_relative_humidity(parent(_stack[:vapor_pressure]),
                                            parent(_stack[:temperature_air]))
        relative_humidity_to_vapor_pressure(parent(_stack[:temperature_air]), fill(80.0, _n))
        dewpoint_to_vapor_pressure(fill(250.0, _n))
        empirical_lapse_rate([255.0, 253.0, 251.0], [1000.0, 1300.0, 1600.0])
    end
end

end # module
