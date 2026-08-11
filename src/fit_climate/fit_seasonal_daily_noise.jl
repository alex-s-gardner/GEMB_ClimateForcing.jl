"""
    fit_seasonal_daily_noise(dec_year, y_data)

Fit yearly/daily sinusoids and autoregressive noise statistics to data.

# Arguments
- `dec_year::Vector{Float64}`: Time points in decimal years
- `y_data::Vector{Float64}`: Observed data values

# Returns
- `NamedTuple` with fields:
  - `beta`: 5×1 vector of linear regression weights
  - `noise_std`: Standard deviation of residuals
  - `noise_lag1`: Lag-1 autocorrelation of residuals

# Details
Fits a model representing data as the linear sum of daily and annual sinusoids,
plus noise and a constant offset:

    y = c1 + c2*cos(Yr) + c3*sin(Yr) + c4*cos(Day) + c5*sin(Day) + noise

where:
- Yr = 2π * dec_year (annual frequency)
- Day = 2π * 365.25 * dec_year (daily frequency)

The residuals are analyzed to estimate noise standard deviation and
lag-1 autocorrelation (for AR1 modeling).

Matches MATLAB's `fit_seasonal_daily_noise.m`.

# References
Original MATLAB implementation in GEMB/src/fit_simulated_climate_to_data/

# See also
`simulate_seasonal_daily_noise`

# Example
```julia
# Fit model to wind speed data
coeffs = fit_seasonal_daily_noise(dec_year, wind_speed)

# Coeffs can be used with simulate_seasonal_daily_noise to generate
# synthetic data with matching statistics
```
"""
function fit_seasonal_daily_noise(dec_year::Vector{Float64},
                                  y_data::Vector{Float64})

    # Ensure inputs are column vectors (already are in Julia, but for clarity)
    t = dec_year
    y = y_data

    ## 1. Construct the Design Matrix (X)

    # Frequency Constants
    omega_yr = 2π           # Once per year
    omega_day = 2π * 365.25 # Approx 365.25 times per year

    # Basis Functions
    # Col 1: Mean (Intercept)
    # Col 2: Daily Cosine
    # Col 3: Daily Sine
    # Col 4: Yearly Cosine
    # Col 5: Yearly Sine
    X = hcat(
        ones(length(t)),
        cos.(omega_day .* t),
        sin.(omega_day .* t),
        cos.(omega_yr .* t),
        sin.(omega_yr .* t)
    )

    ## 2. Linear Regression (Least Squares)
    # Solve y = X*beta
    beta = X \ y

    ## 3. Residual Analysis (Noise Fitting)

    # Calculate the deterministic model prediction
    y_model = X * beta

    # Calculate residuals (What's left over)
    residuals = y .- y_model

    # Calculate Noise Statistics
    noise_std = std(residuals)

    # Calculate Lag-1 Autocorrelation
    # Correlation between x[t] and x[t-1]
    if length(residuals) > 1
        # Use corrcoef equivalent
        r1 = @view residuals[1:end-1]
        r2 = @view residuals[2:end]

        # Pearson correlation
        mean_r1 = mean(r1)
        mean_r2 = mean(r2)

        n = length(r1)
        cov_val = sum((r1 .- mean_r1) .* (r2 .- mean_r2)) / n
        std_r1 = std(r1, corrected=false)
        std_r2 = std(r2, corrected=false)

        if std_r1 > 0 && std_r2 > 0
            noise_lag1 = cov_val / (std_r1 * std_r2)
        else
            noise_lag1 = 0.0
        end
    else
        noise_lag1 = 0.0
    end

    return (
        beta = beta,
        noise_std = noise_std,
        noise_lag1 = noise_lag1
    )
end
