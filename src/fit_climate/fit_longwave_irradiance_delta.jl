"""
    fit_longwave_irradiance_delta(delta)

Fit a Gaussian mixture model to longwave radiation residuals (cloud effects).

# Arguments
- `delta::Vector{Float64}`: Longwave radiation residuals [W/m²]
  (observed - clear_sky_estimate)

# Returns
- `NamedTuple` with fields:
  - `mu`: Mean values of Gaussian components [2×1 vector]
  - `sigma`: Standard deviations of Gaussian components [2×1 vector]
  - `P`: Mixture proportions [2×1 vector]
  - `limits`: [min, max] of observed data
  - `prctile_bounds`: [0.5th, 95th] percentiles for soft limiting

# Details
Fits a 2-component Gaussian Mixture Model to represent the distribution
of longwave radiation anomalies caused by clouds. The bimodal distribution
typically captures clear-sky (low delta) and cloudy-sky (high delta) states.

Note: Julia version uses 2 components like MATLAB code (line 32) rather than
the 3 mentioned in comments.

Matches MATLAB's `fit_longwave_irradiance_delta.m`.

# References
Original MATLAB implementation in GEMB/src/fit_simulated_climate_to_data/

# Example
```julia
# Calculate residuals
clear_sky = simulate_longwave_irradiance(temp, vapor)
delta = observed_longwave .- clear_sky

# Fit model
coeffs = fit_longwave_irradiance_delta(delta)
```
"""
function fit_longwave_irradiance_delta(delta::Vector{Float64})
    # Clean data - remove NaNs
    delta_clean = delta[.!isnan.(delta)]

    if length(delta_clean) < 10
        error("Insufficient data for GMM fitting. Need at least 10 points.")
    end

    # Fit a 2-component Gaussian Mixture Model
    # Note: Requires GaussianMixtures.jl or similar
    # For now, provide a simplified version using k-means initialization

    # Sort data to get initial estimates
    sorted_delta = sort(delta_clean)
    n = length(sorted_delta)

    # Split into two groups for initial estimates
    mid = div(n, 2)
    group1 = sorted_delta[1:mid]
    group2 = sorted_delta[mid+1:end]

    # Initial estimates
    mu1_init = mean(group1)
    mu2_init = mean(group2)
    sigma1_init = std(group1)
    sigma2_init = std(group2)
    p1_init = 0.5
    p2_init = 0.5

    # EM algorithm for GMM fitting
    mu = [mu1_init, mu2_init]
    sigma = [sigma1_init, sigma2_init]
    P = [p1_init, p2_init]

    # Run EM iterations
    max_iter = 100
    tol = 1e-6

    # Responsibility matrix is fully overwritten each E-step, so it is allocated once
    # and reused across iterations rather than per iteration.
    n_samples = length(delta_clean)
    responsibilities = zeros(n_samples, 2)

    for iter in 1:max_iter
        # E-step: Calculate responsibilities. Accumulated in scalars and written
        # once per element — slicing `responsibilities[i, :]` to normalize would
        # allocate two row copies per sample per iteration.
        # `P[k] * (pdf)` keeps the baseline's multiplication grouping — writing it as
        # `(P[k] * norm) * exp(...)` reassociates and shifts the fitted coefficients
        # in the last few ulps.
        @inbounds for i in 1:n_samples
            x = delta_clean[i]
            r1 = P[1] * ((1.0 / (sigma[1] * sqrt(2π))) * exp(-0.5 * ((x - mu[1]) / sigma[1])^2))
            r2 = P[2] * ((1.0 / (sigma[2] * sqrt(2π))) * exp(-0.5 * ((x - mu[2]) / sigma[2])^2))
            row_sum = r1 + r2
            if row_sum > 0
                r1 /= row_sum
                r2 /= row_sum
            end
            responsibilities[i, 1] = r1
            responsibilities[i, 2] = r2
        end

        # M-step: Update parameters
        mu_old = copy(mu)

        for k in 1:2
            r_k = @view responsibilities[:, k]
            N_k = sum(r_k)

            if N_k > 1e-10
                # These two reductions keep their allocating broadcast form on purpose.
                # A single-pass accumulation loop is faster, but any such loop — with or
                # without `@simd` — reassociates the additions away from the pairwise
                # order `sum` uses over an `Array`, and this is MATLAB-parity code where
                # a ~1e-14 drift in the fitted coefficients is a regression. (`sum(f, idx)`
                # is no help: reducing lazily over indices blocks differently from
                # reducing over a materialized array, so it is not bit-identical either.)
                # The cost is bounded — six n-element temporaries per iteration, against
                # the E-step above, which is where this function's time actually goes.
                mu[k] = sum(r_k .* delta_clean) / N_k
                sigma[k] = sqrt(sum(r_k .* (delta_clean .- mu[k]) .^ 2) / N_k)
                P[k] = N_k / n_samples
            end
        end

        # Check convergence
        if maximum(abs.(mu .- mu_old)) < tol
            break
        end
    end

    # Add regularization to prevent singular covariance
    sigma .= max.(sigma, 0.1)

    # Normalize proportions
    P ./= sum(P)

    # Statistical range anchors
    limits = [minimum(delta_clean), maximum(delta_clean)]

    # Percentile bounds (soft limits)
    sorted_clean = sort(delta_clean)
    idx_low = max(1, round(Int, 0.005 * length(sorted_clean)))
    idx_high = min(length(sorted_clean), round(Int, 0.95 * length(sorted_clean)))
    prctile_bounds = [sorted_clean[idx_low], sorted_clean[idx_high]]

    return (
        mu = mu,
        sigma = sigma,
        P = P,
        limits = limits,
        prctile_bounds = prctile_bounds
    )
end
