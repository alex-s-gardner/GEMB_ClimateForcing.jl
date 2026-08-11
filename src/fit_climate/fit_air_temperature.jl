"""
    fit_air_temperature(dec_year, temperature_air, latitude, elevation)

Estimate simulation coefficients for air temperature from observed data.

# Arguments
- `dec_year::Vector{Float64}`: Decimal years (e.g., 2024.0, 2024.002)
- `temperature_air::Vector{Float64}`: Observed air temperatures [K]
- `latitude::Float64`: Latitude [degrees]
- `elevation::Float64`: Elevation [m]

# Returns
- `NamedTuple` with fields:
  - `mean_offset`: Temperature offset from theoretical mean [K]
  - `lat_scale`: Latitude scaling factor for seasonal amplitude
  - `daily_amp_scale`: Daily amplitude scaling factor
  - `weather_sigma_scale`: Weather noise standard deviation scale
  - `weather_corr`: Lag-1 autocorrelation of weather noise

# Details
Fits a harmonic model to observed temperature data to extract parameters
for `simulate_air_temperature`. The model accounts for:
1. Mean temperature (latitude + elevation dependent)
2. Annual cycle (latitude dependent amplitude)
3. Diurnal cycle (elevation dependent amplitude)
4. Weather noise (AR1 process with daily averaging)

Matches MATLAB's `fit_air_temperature.m`.

# References
Original MATLAB implementation in GEMB/src/fit_simulated_climate_to_data/
"""
function fit_air_temperature(dec_year::Vector{Float64},
                             temperature_air::Vector{Float64},
                             latitude::Float64,
                             elevation::Float64)

    # Constants from the simulation model
    BASE_SIGMA = 8.0
    LAPSE_RATE = 0.0065

    ## 1. Calculate Mean Offset
    # The model calculates a theoretical mean based on latitude/elevation.
    # We compare the observed mean to this theoretical mean.

    phi = deg2rad(latitude)
    T_sea_level_theoretical = 300.0 - 50.0 * sin(phi)^2
    T_mean_theoretical = T_sea_level_theoretical - (elevation * LAPSE_RATE)

    T_obs_mean = mean(temperature_air)

    # RESULT 1: mean_offset
    mean_offset = T_obs_mean - T_mean_theoretical

    ## 2. Harmonic Analysis (Least Squares Fit)
    # Remove the mean and fit the Annual and Diurnal cosines

    T_anomaly = temperature_air .- T_obs_mean

    # A. Construct Annual Basis Function
    # Match simulation logic: cos(2*pi * (year_frac - phase))
    # Phase is 0.5 for North, 0 for South
    phase_annual = latitude > 0 ? 0.5 : 0.0

    year_frac = dec_year .- floor.(dec_year)
    basis_annual = cos.(2π .* (year_frac .- phase_annual))

    # B. Construct Diurnal Basis Function
    # Match simulation logic: cos(2*pi * (day_fraction - 0.625))
    day_fraction = mod.(dec_year .* 365.25, 1.0)
    basis_diurnal = cos.(2π .* (day_fraction .- 0.625))

    # C. Solve Linear Regression: Anomaly = c1*Annual + c2*Diurnal
    X = hcat(basis_annual, basis_diurnal)
    betas = X \ T_anomaly

    fitted_amp_annual = betas[1]
    fitted_amp_daily = betas[2]

    ## 3. Derive Scale Coefficients

    # --- lat_scale ---
    # Model: Amp = 3 + (22 * |sin(latitude)| * lat_scale)
    # Inverse: lat_scale = (fitted_amp - 3) / (22 * |sin(latitude)|)
    denom_annual = 22.0 * abs(sin(phi))

    if denom_annual < 1e-4
        # Handle Equator case (sin(0)=0) to avoid Inf
        lat_scale = 1.0
    else
        lat_scale = (fitted_amp_annual - 3.0) / denom_annual
    end
    # Enforce bounds (scale cannot be negative)
    lat_scale = max(0.0, lat_scale)

    # --- daily_amp_scale ---
    # Model: Amp = (DTR_base / 2) * daily_amp_scale
    dtr_base = 10.0 + (elevation / 1000.0)
    base_amp_daily = dtr_base / 2.0

    daily_amp_scale = fitted_amp_daily / base_amp_daily
    daily_amp_scale = max(0.0, daily_amp_scale)

    ## 4. Estimate Weather Parameters (Residual Analysis)
    # Remove the fitted deterministic cycles to get pure weather noise

    T_deterministic = (fitted_amp_annual .* basis_annual) .+ (fitted_amp_daily .* basis_diurnal)
    residuals = T_anomaly .- T_deterministic

    # The simulation assumes weather is a DAILY process (AR1).
    # We must average residuals by day to estimate these parameters accurately.

    day_indices = floor.(Int, dec_year .* 365.25)

    # Compute daily mean of residuals
    unique_days = sort(unique(day_indices))
    daily_res = Float64[]
    for day in unique_days
        mask = day_indices .== day
        push!(daily_res, mean(residuals[mask]))
    end

    # --- weather_sigma_scale ---
    # Model: sigma = 8.0 * scale
    std_daily = std(daily_res)
    weather_sigma_scale = std_daily / BASE_SIGMA

    # --- weather_corr ---
    # Calculate lag-1 autocorrelation
    if length(daily_res) > 2
        # Calculate correlation between x[t] and x[t+1]
        n = length(daily_res) - 1
        x1 = @view daily_res[1:n]
        x2 = @view daily_res[2:end]

        # Pearson correlation
        mean_x1 = mean(x1)
        mean_x2 = mean(x2)

        cov_val = sum((x1 .- mean_x1) .* (x2 .- mean_x2)) / n
        std_x1 = std(x1, corrected=false)
        std_x2 = std(x2, corrected=false)

        raw_corr = cov_val / (std_x1 * std_x2)

        # Clamp between 0 and 0.99
        weather_corr = clamp(raw_corr, 0.0, 0.99)
    else
        # Fallback for insufficient data
        weather_corr = 0.7
    end

    return (
        mean_offset = mean_offset,
        lat_scale = lat_scale,
        daily_amp_scale = daily_amp_scale,
        weather_sigma_scale = weather_sigma_scale,
        weather_corr = weather_corr
    )
end
