"""
    fit_precipitation(dec_year, precipitation; wet_threshold=0.1)

Fit a seasonal Markov-Gamma model to hourly precipitation data.

# Arguments
- `dec_year::Vector{Float64}`: Decimal years (e.g., 2021.45)
- `precipitation::Vector{Float64}`: Precipitation amounts [mm]
- `wet_threshold::Float64=0.1`: Minimum amount to consider "wet" [mm]

# Returns
- `NamedTuple` with fields:
  - `P01_harmonics`: 3×1 vector of dry→wet transition probability harmonics
  - `P11_harmonics`: 3×1 vector of wet→wet transition probability harmonics
  - `Alpha_harmonics`: 3×1 vector of Gamma shape parameter harmonics
  - `Beta_harmonics`: 3×1 vector of Gamma scale parameter harmonics
  - `wet_threshold`: Threshold used for wet/dry classification

# Details
Fits a two-state Markov chain (wet/dry) with seasonally-varying transition
probabilities and Gamma-distributed precipitation amounts when wet.

The model:
1. Bins data into 12 monthly bins
2. Calculates Markov transition probabilities (P01: dry→wet, P11: wet→wet)
3. Fits Gamma distribution to wet amounts using method of moments
4. Fits harmonic functions to capture seasonality: Y = C1 + C2*sin(2πt) + C3*cos(2πt)

Matches MATLAB's `fit_precipitation.m`.

# References
Original MATLAB implementation in GEMB/src/fit_simulated_climate_to_data/
"""
function fit_precipitation(dec_year::Vector{Float64},
                          precipitation::Vector{Float64};
                          wet_threshold::Float64=0.1)

    ## 1. Setup

    # Extract fractional year (0 to 1) for seasonality
    t_season = mod.(dec_year, 1.0)

    # Determine Wet/Dry states
    is_wet = precipitation .>= wet_threshold

    ## 2. Calculate Monthly Statistics (12 bins)

    months = ceil.(Int, t_season .* 12.0)
    months[months .== 0] .= 1  # Handle edge case exactly at 0.0

    # Columns: [P01, P11, Alpha, Beta]
    stats = zeros(12, 4)

    for m in 1:12
        idx = months .== m
        n_obs = sum(idx)

        if n_obs < 2
            continue
        end

        data_m = is_wet[idx]
        precip_m = precipitation[idx]

        # --- Markov Transition Probabilities ---
        # Current state vs Next state (shift by 1)
        current = data_m[1:end-1]
        next = data_m[2:end]

        # P01: Probability Dry → Wet (Count 0→1 / Count 0)
        idx_0 = .!current
        if sum(idx_0) > 0
            stats[m, 1] = sum(next[idx_0]) / sum(idx_0)
        end

        # P11: Probability Wet → Wet (Count 1→1 / Count 1)
        idx_1 = current
        if sum(idx_1) > 0
            stats[m, 2] = sum(next[idx_1]) / sum(idx_1)
        end

        # --- Gamma Distribution Parameters (Method of Moments) ---
        wet_amts = precip_m[precip_m .>= wet_threshold]
        if length(wet_amts) > 5
            mu = mean(wet_amts)
            v = var(wet_amts)

            # Gamma params: Alpha (Shape) = mean^2 / var, Beta (Scale) = var / mean
            if v > 0
                stats[m, 3] = mu^2 / v  # Alpha
                stats[m, 4] = v / mu    # Beta
            end
        end
    end

    ## 3. Fit Harmonic Curves
    # Model: Y = C1 + C2*sin(2*pi*t) + C3*cos(2*pi*t)
    # We fit this to the 12 monthly values

    t_bins = collect(0.5:1:11.5) ./ 12.0  # Midpoints of months
    X = hcat(ones(12), sin.(2π .* t_bins), cos.(2π .* t_bins))

    # Solve least squares: B = (X'X)^-1 X'Y
    # Do this for all 4 parameters at once
    B = X \ stats

    # Pack into NamedTuple
    return (
        P01_harmonics = B[:, 1],
        P11_harmonics = B[:, 2],
        Alpha_harmonics = B[:, 3],
        Beta_harmonics = B[:, 4],
        wet_threshold = wet_threshold
    )
end
