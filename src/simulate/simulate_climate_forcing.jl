# Synthetic climate forcing generation
# Translated from MATLAB simulate_* functions for GEMB.jl

using Random  # for MersenneTwister, randn, rand

#=============================================================================
# Parameter Sets
=============================================================================#

"""
    simulation_parameter_sets(set_id::String)

Retrieve predefined parameter sets for climate forcing simulations.
Returns `(location_parameters, coeffs)` as NamedTuples.

Matches MATLAB's `simulation_parameter_sets.m`.
"""
function simulation_parameter_sets(set_id::String)
    valid_sets = ["test_1"]

    if set_id == "test_1"
        location_parameters = (
            description = "parameters estimated using simulation_parameters_estimate_from_data.m as fit to original TEST_INPUT_1.mat data",
            latitude = -73.3307,                         # [deg]
            longitude = 290.6250,                        # [deg]
            elevation = 700.0,                           # [m]
            start_year = 1994,                           # [year]
            end_year = 2025,                             # [year]
            wind_observation_height = 10.0,              # [m]
            temperature_observation_height = 2.0,        # [m]
            temperature_air_mean = 259.4,                # [K]
            precipitation_mean = 1177.3,                 # [kg m-2 yr-1]
            time_step_hours = 1.0,                       # [hours]
            rand_seed = 42,                              # seed for RNG
        )

        coeffs = (
            temperature_air = (
                mean_offset = 9.8847,
                lat_scale = 0.1635,
                daily_amp_scale = 0.0000,
                weather_sigma_scale = 0.8418,
                weather_corr = 0.7315,
            ),
            relative_humidity = (
                beta = [79.9128, -0.2619, 1.3055, 5.0727, 0.5081],
                noise_std = 11.9734,
                noise_lag1 = 0.9896,
                min_max = [0.0, 100.0],
            ),
            longwave_downward = (
                mu = [113.2225, 23.6664],
                sigma = [30.9535, 38.1762],
                P = [0.3590, 0.6410],
                limits = [-142.9820, 222.0189],
                prctile_bounds = [-67.5562, 147.9489],
                min_max = [0.0, Inf],
            ),
            wind_speed = (
                beta = [5.2495, 0.1296, -0.1237, -0.3030, -0.2289],
                noise_std = 2.9324,
                noise_lag1 = 0.9800,
                min_max = [0.0, Inf],
            ),
            precipitation = (
                P01_harmonics = [0.0213, -0.0002, -0.0019],
                P11_harmonics = [0.9443, -0.0013, -0.0036],
                Alpha_harmonics = [1.1357, 0.0610, 0.1906],
                Beta_harmonics = [0.4060, -0.0376, -0.1312],
                wet_threshold = 0.1000,
            ),
        )
    else
        msg = join(valid_sets, ", ")
        error("\"$set_id\" is not a valid simulation parameter set. Valid sets include: [$msg]")
    end

    return location_parameters, coeffs
end

#=============================================================================
# Helper: DateTime <-> Decimal year conversions
=============================================================================#

"""
    decyear2datetime(decyear::AbstractVector{<:Real})

Convert decimal year values to Julia `DateTime` objects.
Matches MATLAB's `decyear2datenum.m` logic (accounts for leap years).
"""
function decyear2datetime(decyear::AbstractVector{<:Real})
    datetimes = Vector{DateTime}(undef, length(decyear))
    for i in eachindex(decyear)
        year_part = floor(Int, decyear[i])
        fractional_year = decyear[i] - year_part

        start_of_year = DateTime(year_part, 1, 1)
        start_of_next_year = DateTime(year_part + 1, 1, 1)
        days_in_year = Dates.value(start_of_next_year - start_of_year) / (1000 * 86400)  # milliseconds to days

        day_offset = fractional_year * days_in_year
        datetimes[i] = start_of_year + Dates.Millisecond(round(Int64, day_offset * 86400 * 1000))
    end
    return datetimes
end

"""
    datetime2decyear(datetimes::AbstractVector{DateTime})

Convert Julia `DateTime` objects to decimal year values.
"""
function datetime2decyear(datetimes::AbstractVector{DateTime})
    dec_year = Vector{Float64}(undef, length(datetimes))
    for i in eachindex(datetimes)
        yr = Dates.year(datetimes[i])
        start_of_year = DateTime(yr, 1, 1)
        start_of_next_year = DateTime(yr + 1, 1, 1)
        days_in_year = Dates.value(start_of_next_year - start_of_year) / (1000.0 * 86400.0)
        day_offset = Dates.value(datetimes[i] - start_of_year) / (1000.0 * 86400.0)
        dec_year[i] = yr + day_offset / days_in_year
    end
    return dec_year
end

#=============================================================================
# simulate_air_temperature
=============================================================================#

"""
    simulate_air_temperature(rng, dec_year, latitude, elevation; kwargs...)

Simulate near-surface air temperature [K] as a function of time, latitude, and elevation.
Uses seasonal cycle, diurnal cycle, and AR(1) weather noise.

Matches MATLAB's `simulate_air_temperature.m`.
"""
function simulate_air_temperature(
    rng::AbstractRNG,
    dec_year::Vector{Float64},
    latitude::Float64,
    elevation::Float64;
    lat_scale::Float64=1.0,
    daily_amp_scale::Float64=1.0,
    weather_sigma_scale::Float64=1.0,
    weather_corr::Float64=0.7,
    mean_offset::Float64=0.0
)
    # 2. Climatology (Mean Temp)
    phi = deg2rad(latitude)
    T_sea_level = 300.0 - 50.0 * sin(phi)^2
    LapseRate = 0.0065
    temperature_air_mean = T_sea_level - (elevation * LapseRate) + mean_offset

    # 3. Seasonal Cycle (Annual Wave)
    T_amp_annual = 3.0 + (22.0 * abs(sin(phi)) * lat_scale)

    # Hemisphere logic
    phase_shift_annual = latitude > 0 ? 0.5 : 0.0

    year_frac = dec_year .- floor.(dec_year)
    seasonal_signal = cos.(2π .* (year_frac .- phase_shift_annual))

    # 4. Diurnal Cycle (Daily Wave)
    day_fraction = mod.(dec_year .* 365.25, 1.0)
    DTR = 10.0 + (elevation / 1000.0)
    T_amp_daily = (DTR / 2.0) * daily_amp_scale
    diurnal_signal = cos.(2π .* (day_fraction .- 0.625))

    # 5. Synoptic Weather (Correlated Noise)
    start_day = floor(Int, minimum(dec_year) * 365.25)
    end_day = ceil(Int, maximum(dec_year) * 365.25)
    num_days = end_day - start_day + 1

    base_sigma = 8.0
    sigma_weather = base_sigma * weather_sigma_scale

    # Generate white noise scaled to maintain variance over time
    white_noise = randn(rng, num_days) .* (sigma_weather * sqrt(1.0 - weather_corr^2))

    # AR(1) process
    daily_noise = zeros(num_days)
    daily_noise[1] = white_noise[1]
    for t in 2:num_days
        daily_noise[t] = (weather_corr * daily_noise[t-1]) + white_noise[t]
    end

    # Interpolate daily noise to user time steps
    user_day_indices = (dec_year .* 365.25) .- start_day .+ 1.0

    # Linear interpolation (matches MATLAB interp1 'linear','extrap')
    weather_signal = _interp1_linear(1:num_days, daily_noise, user_day_indices)

    # 6. Final Combination
    temperature_air = temperature_air_mean .+
        (T_amp_annual .* seasonal_signal) .+
        (T_amp_daily .* diurnal_signal) .+
        weather_signal

    return temperature_air
end

#=============================================================================
# simulate_air_pressure
=============================================================================#

"""
    simulate_air_pressure(rng, dec_year, temperature_air, latitude, elevation)

Simulate screen-level atmospheric pressure [Pa] using the hypsometric equation.
Includes stochastic MSLP fluctuations via AR(1) process.

Matches MATLAB's `simulate_air_pressure.m`.
"""
function simulate_air_pressure(
    rng::AbstractRNG,
    dec_year::Vector{Float64},
    temperature_air::Vector{Float64},
    latitude::Float64,
    elevation::Float64
)
    # 2. Simulate Mean Sea Level Pressure Weather Patterns
    phi = deg2rad(latitude)
    sigma_pressure = 300.0 + (900.0 * sin(phi)^2)

    alpha = 0.85

    start_day = floor(Int, minimum(dec_year) * 365.25)
    end_day = ceil(Int, maximum(dec_year) * 365.25)
    num_days = end_day - start_day + 1

    white_noise = randn(rng, num_days) .* (sigma_pressure * sqrt(1.0 - alpha^2))

    daily_noise = zeros(num_days)
    daily_noise[1] = white_noise[1]
    for t in 2:num_days
        daily_noise[t] = (alpha * daily_noise[t-1]) + white_noise[t]
    end

    # Map daily noise to user time steps
    user_day_indices = (dec_year .* 365.25) .- start_day .+ 1.0
    weather_anomaly = _interp1_linear(1:num_days, daily_noise, user_day_indices)

    # Final Simulated MSLP
    P_msl_std = 101325.0
    P_msl = P_msl_std .+ weather_anomaly

    # 3. Calculate Local Station Pressure (Hypsometric Reduction)
    LapseRate = 0.0065
    T_column_avg = temperature_air .+ (0.5 * LapseRate * elevation)

    g = 9.81
    R = 287.05

    exponent = (-g * elevation) ./ (R .* T_column_avg)
    pressure_air = P_msl .* exp.(exponent)

    return pressure_air
end

#=============================================================================
# simulate_precipitation
=============================================================================#

"""
    simulate_precipitation(rng, dec_year, coeffs)

Generate synthetic precipitation [kg m-2] using a Markov chain occurrence model
with gamma-distributed amounts.

Matches MATLAB's `simulate_precipitation.m`.
"""
function simulate_precipitation(
    rng::AbstractRNG,
    dec_year::Vector{Float64},
    coeffs::NamedTuple
)
    n = length(dec_year)
    precipitation = zeros(n)

    # 1. Reconstruct Time-Varying Parameters
    t_season = mod.(dec_year, 1.0)

    # Basis functions: [1, sin(2*pi*t), cos(2*pi*t)]
    X = hcat(ones(n), sin.(2π .* t_season), cos.(2π .* t_season))

    # Calculate params and clamp to valid ranges
    P01_t = X * coeffs.P01_harmonics
    P11_t = X * coeffs.P11_harmonics
    Alpha_t = X * coeffs.Alpha_harmonics
    Beta_t = X * coeffs.Beta_harmonics

    # Clamping
    P01_t = clamp.(P01_t, 0.0, 1.0)
    P11_t = clamp.(P11_t, 0.0, 1.0)
    Alpha_t = max.(Alpha_t, 0.1)
    Beta_t = max.(Beta_t, 0.1)

    # 2. Simulation Loop
    # Initial state based on first P01
    is_raining = rand(rng) < P01_t[1]

    for i in 1:n
        # Determine transition probability based on previous state
        prob_wet = is_raining ? P11_t[i] : P01_t[i]

        # Update State
        if rand(rng) < prob_wet
            is_raining = true
            # Generate Amount (Gamma Distribution)
            precipitation[i] = _gamrnd(rng, Alpha_t[i], Beta_t[i])
            # Enforce threshold floor
            if precipitation[i] < coeffs.wet_threshold
                precipitation[i] = coeffs.wet_threshold
            end
        else
            is_raining = false
            precipitation[i] = 0.0
        end
    end

    return precipitation
end

#=============================================================================
# simulate_shortwave_irradiance
=============================================================================#

"""
    simulate_shortwave_irradiance(dec_year, latitude)

Simulate clear-sky downwelling shortwave irradiance [W m-2] using the Haurwitz model.

Matches MATLAB's `simulate_shortwave_irradiance.m`.
"""
function simulate_shortwave_irradiance(dec_year::Vector{Float64}, latitude::Float64)
    n = length(dec_year)

    # 1. Time Conversion (Decimal Year -> Day & Hour)
    year_val = floor.(Int, dec_year)
    frac_year = dec_year .- year_val

    # Check for leap year
    is_leap = ((mod.(year_val, 4) .== 0) .& (mod.(year_val, 100) .!= 0)) .|
              (mod.(year_val, 400) .== 0)
    days_in_year = 365.0 .+ Float64.(is_leap)

    # Continuous Day of Year
    doy_continuous = frac_year .* days_in_year .+ 1.0

    # Integer Day Number for Declination
    n_day = floor.(Int, doy_continuous)

    # Solar Hour (0 to 24)
    solar_hour = (doy_continuous .- n_day) .* 24.0

    # 2. Solar Geometry
    phi = deg2rad(latitude)

    # Solar Declination (Cooper's Equation)
    delta = deg2rad.(23.45 .* sind.((360.0 ./ 365.25) .* (284.0 .+ n_day)))

    # Hour Angle
    omega = deg2rad.(15.0 .* (solar_hour .- 12.0))

    # Cosine of Solar Zenith Angle
    cos_theta_z = (sin(phi) .* sin.(delta)) .+ (cos(phi) .* cos.(delta) .* cos.(omega))

    # 3. Calculate Irradiance (Haurwitz Model)
    shortwave_downward = zeros(n)
    daylight_mask = cos_theta_z .> 0.0

    for i in 1:n
        if daylight_mask[i]
            ctz = cos_theta_z[i]
            shortwave_downward[i] = 1098.0 * ctz * exp(-0.057 / ctz)
        end
    end

    return shortwave_downward
end

#=============================================================================
# simulate_longwave_irradiance
=============================================================================#

"""
    simulate_longwave_irradiance(temperature_air, vapor_pressure)

Estimate downward longwave radiation [W m-2] using Brutsaert's (1975) parameterization.

Matches MATLAB's `simulate_longwave_irradiance.m`.
"""
function simulate_longwave_irradiance(temperature_air::Vector{Float64}, vapor_pressure::Vector{Float64})
    sigma = 5.670374419e-8  # Stefan-Boltzmann constant [W m-2 K-4]

    # Convert Pa -> hPa for Brutsaert's coefficient
    e_hPa = vapor_pressure ./ 100.0

    # Effective emissivity (Brutsaert 1975, clear sky)
    epsilon_clear = 1.24 .* (e_hPa ./ temperature_air) .^ (1.0 / 7.0)

    # Stefan-Boltzmann Law
    longwave_downward = epsilon_clear .* sigma .* (temperature_air .^ 4)

    return longwave_downward
end

#=============================================================================
# simulate_longwave_irradiance_delta
=============================================================================#

"""
    simulate_longwave_irradiance_delta(rng, dec_year, coeff)

Generate longwave radiation anomalies from a Gaussian mixture model with truncation.

Matches MATLAB's `simulate_longwave_irradiance_delta.m`.
"""
function simulate_longwave_irradiance_delta(
    rng::AbstractRNG,
    dec_year::Vector{Float64},
    coeff::NamedTuple
)
    n = length(dec_year)
    longwave_downward_delta = zeros(n)

    # Cumulative probability for component selection
    cumP = cumsum(coeff.P)

    # Pre-generate component selectors
    comp_selector = rand(rng, n)

    # Sample and Truncate
    for i in 1:n
        # Determine which component to sample
        idx = findfirst(x -> comp_selector[i] <= x, cumP)

        # Draw until value is within observed 'soft' bounds
        val = coeff.mu[idx] + coeff.sigma[idx] * randn(rng)
        while val < coeff.prctile_bounds[1] || val > coeff.prctile_bounds[2]
            val = coeff.mu[idx] + coeff.sigma[idx] * randn(rng)
        end
        longwave_downward_delta[i] = val
    end

    return longwave_downward_delta
end

#=============================================================================
# simulate_seasonal_daily_noise
=============================================================================#

"""
    simulate_seasonal_daily_noise(rng, dec_year, coeffs)

Generate data from fitted seasonal + daily deterministic signal with AR(1) correlated noise.

Matches MATLAB's `simulate_seasonal_daily_noise.m`.
"""
function simulate_seasonal_daily_noise(
    rng::AbstractRNG,
    dec_year::Vector{Float64},
    coeffs::NamedTuple
)
    t = dec_year
    n = length(t)

    # 1. Reconstruct Deterministic Component
    omega_yr = 2π
    omega_day = 2π * 365.25

    # Design Matrix: [1, cos(omega_day*t), sin(omega_day*t), cos(omega_yr*t), sin(omega_yr*t)]
    X = hcat(
        ones(n),
        cos.(omega_day .* t),
        sin.(omega_day .* t),
        cos.(omega_yr .* t),
        sin.(omega_yr .* t)
    )

    # Calculate pure signal
    y_deterministic = X * coeffs.beta

    # 2. Generate Synthetic Correlated Noise (AR(1))
    sigma = coeffs.noise_std
    phi = coeffs.noise_lag1

    # Scaling for driving white noise to match target std
    white_noise_scale = sigma * sqrt(1.0 - phi^2)

    # Generate white noise
    u = randn(rng, n) .* white_noise_scale

    # AR(1) loop
    noise_loop = zeros(n)
    noise_loop[1] = randn(rng) * sigma
    for i in 2:n
        noise_loop[i] = phi * noise_loop[i-1] + u[i]
    end

    # 3. Combine
    y_sim = y_deterministic .+ noise_loop

    return y_sim
end

#=============================================================================
# Helper: Linear interpolation (matches MATLAB interp1 linear with extrap)
=============================================================================#

"""
    _interp1_linear(x, y, xi)

Linear interpolation with extrapolation, matching MATLAB's `interp1(x, y, xi, 'linear', 'extrap')`.
"""
function _interp1_linear(x::AbstractRange, y::Vector{Float64}, xi::Vector{Float64})
    n = length(x)
    result = similar(xi)
    x1 = Float64(first(x))
    dx = Float64(step(x))

    for i in eachindex(xi)
        # Fractional index (0-based)
        frac_idx = (xi[i] - x1) / dx
        idx_lo = floor(Int, frac_idx) + 1  # Convert to 1-based

        if idx_lo < 1
            # Extrapolate below
            idx_lo = 1
            idx_hi = 2
        elseif idx_lo >= n
            # Extrapolate above
            idx_lo = n - 1
            idx_hi = n
        else
            idx_hi = idx_lo + 1
        end

        # Linear interpolation weight
        t = (xi[i] - (x1 + (idx_lo - 1) * dx)) / dx
        result[i] = y[idx_lo] + t * (y[idx_hi] - y[idx_lo])
    end

    return result
end

#=============================================================================
# Helper: Gamma random variate (shape, scale parameterization)
=============================================================================#

"""
    _gamrnd(rng, shape, scale)

Generate a gamma-distributed random variate with given shape (alpha) and scale (beta).
Uses Marsaglia and Tsang's method for shape >= 1, and transformation for shape < 1.
Matches MATLAB's `gamrnd(shape, scale)`.
"""
function _gamrnd(rng::AbstractRNG, shape::Float64, scale::Float64)
    if shape < 1.0
        # For shape < 1, use: X = Gamma(shape+1) * U^(1/shape)
        return _gamrnd(rng, shape + 1.0, scale) * rand(rng)^(1.0 / shape)
    end

    # Marsaglia and Tsang's method for shape >= 1
    d = shape - 1.0 / 3.0
    c = 1.0 / sqrt(9.0 * d)

    while true
        x = randn(rng)
        v = (1.0 + c * x)^3
        if v > 0.0
            u = rand(rng)
            if u < 1.0 - 0.0331 * x^4 || log(u) < 0.5 * x^2 + d * (1.0 - v + log(v))
                return d * v * scale
            end
        end
    end
end

#=============================================================================
# Main: simulate_climate_forcing
=============================================================================#

"""
    simulate_climate_forcing(set_id::String, time_step_hours::Int=0)

Reproducibly generate synthetic climate forcing data for GEMB simulations
based on predefined parameter sets.

# Arguments
- `set_id::String`: Identifier for the simulation parameter set (e.g., "test_1").
- `time_step_hours::Int`: Temporal resolution in hours. If 0 (default), uses the
  default time step defined by the parameter set.

# Returns
A `DimStack` containing synthetic time series and metadata, suitable for conversion
to `GEMB.ClimateForcing` via `GEMB.ClimateForcing(ds)` when GEMB.jl is loaded.

Matches MATLAB's `simulate_climate_forcing.m`.
"""
function simulate_climate_forcing(set_id::String, time_step_hours::Int=0)
    # Load climate simulation parameter set
    location_parameters, coeffs = simulation_parameter_sets(set_id)

    # Determine time step in hours
    if time_step_hours == 0
        dt_hours = location_parameters.time_step_hours
    else
        dt_hours = Float64(time_step_hours)
    end

    # Create DateTime time vector directly. The interval is half-open:
    # [start_year-01-01, (end_year+1)-01-01), i.e. the trailing boundary instant
    # (end_year+1)-01-01T00:00 is excluded. Including it would append a lone
    # orphan step in a new day/year that pollutes daily output binning; dropping
    # it yields whole days only (e.g. 32 years × 365.25 d × 8 steps).
    dt_start = DateTime(location_parameters.start_year, 1, 1)
    dt_end = DateTime(location_parameters.end_year + 1, 1, 1)

    # Generate time vector using Dates arithmetic
    dt_step = Dates.Millisecond(round(Int64, dt_hours * 3600 * 1000))
    time_vec = collect(dt_start:dt_step:(dt_end - dt_step))

    # Convert DateTime to decimal year for internal simulation functions
    dec_year = datetime2decyear(time_vec)

    # Initialize RNG with seed (Mersenne Twister, same as MATLAB default)
    rng = Random.MersenneTwister(location_parameters.rand_seed)

    # Simulate downward shortwave radiation (deterministic, no RNG needed)
    shortwave_downward = simulate_shortwave_irradiance(dec_year, location_parameters.latitude)

    # Simulate air temperature [K]
    temperature_air = simulate_air_temperature(
        rng, dec_year, location_parameters.latitude, location_parameters.elevation;
        mean_offset=coeffs.temperature_air.mean_offset,
        lat_scale=coeffs.temperature_air.lat_scale,
        daily_amp_scale=coeffs.temperature_air.daily_amp_scale,
        weather_sigma_scale=coeffs.temperature_air.weather_sigma_scale,
        weather_corr=coeffs.temperature_air.weather_corr
    )

    # Simulate screen-level air pressure [Pa]
    pressure_air = simulate_air_pressure(
        rng, dec_year, temperature_air, location_parameters.latitude, location_parameters.elevation
    )

    # Simulate screen-level relative humidity [%]
    relative_humidity = simulate_seasonal_daily_noise(rng, dec_year, coeffs.relative_humidity)
    clamp!(relative_humidity, coeffs.relative_humidity.min_max[1], coeffs.relative_humidity.min_max[2])

    # Screen-level vapor pressure [Pa]
    vapor_pressure = relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)

    # Downward longwave radiation [W m-2]
    longwave_downward = simulate_longwave_irradiance(temperature_air, vapor_pressure)
    longwave_downward .+= simulate_longwave_irradiance_delta(rng, dec_year, coeffs.longwave_downward)
    clamp!(longwave_downward, coeffs.longwave_downward.min_max[1], coeffs.longwave_downward.min_max[2])

    # Screen-level wind speed [m s-1]
    wind_speed = simulate_seasonal_daily_noise(rng, dec_year, coeffs.wind_speed)
    clamp!(wind_speed, coeffs.wind_speed.min_max[1], coeffs.wind_speed.min_max[2])

    # Precipitation [kg m-2 per timestep]
    precipitation = simulate_precipitation(rng, dec_year, coeffs.precipitation)

    # Build DimStack with all climate variables (same shape as ERA5-Land output)
    time_dim = Ti(time_vec)

    stack = DimStack((
        temperature_air = DimArray(temperature_air, (time_dim,);
                                  metadata=Dict("units" => "K", "long_name" => "2m air temperature")),
        pressure_air = DimArray(pressure_air, (time_dim,);
                               metadata=Dict("units" => "Pa", "long_name" => "surface pressure")),
        vapor_pressure = DimArray(vapor_pressure, (time_dim,);
                                 metadata=Dict("units" => "Pa", "long_name" => "vapor pressure")),
        wind_speed = DimArray(wind_speed, (time_dim,);
                             metadata=Dict("units" => "m/s", "long_name" => "wind speed")),
        precipitation = DimArray(precipitation, (time_dim,);
                                metadata=Dict("units" => "kg/m²", "long_name" => "total precipitation per timestep")),
        shortwave_downward = DimArray(shortwave_downward, (time_dim,);
                                     metadata=Dict("units" => "W/m²", "long_name" => "surface solar radiation downward")),
        longwave_downward = DimArray(longwave_downward, (time_dim,);
                                    metadata=Dict("units" => "W/m²", "long_name" => "surface thermal radiation downward"))
    ); metadata = Dict(
        "latitude" => location_parameters.latitude,
        "longitude" => location_parameters.longitude,
        "dataset" => "synthetic",
        "temperature_air_mean" => location_parameters.temperature_air_mean,
        "wind_speed_mean" => Statistics.mean(wind_speed),
        "precipitation_mean" => location_parameters.precipitation_mean,  # annual [kg m-2 yr-1]
        "temperature_observation_height" => location_parameters.temperature_observation_height,
        "wind_observation_height" => location_parameters.wind_observation_height
    ))

    return stack
end
