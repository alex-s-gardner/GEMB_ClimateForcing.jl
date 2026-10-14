# Solar geometry and the shortwave forcing layers derived from it.
#
# A forcing step stamped `t` carries the mean flux over its accumulation interval `(t − Δt, t]`
# (ERA5-Land's hourly fluxes are means over the preceding hour), so every geometric quantity here
# is an average over that interval rather than a value at its stamp:
#
#   solar_zenith_angle          — degrees; `acos` of the irradiance-weighted mean cosine of the
#                                 zenith angle over the sunlit part of the interval
#   shortwave_downward_diffuse  — W/m²; diffuse part of the downward shortwave
#   cloud_optical_thickness     — dimensionless
#
# The weighting matters at low sun: an unweighted mean cosine over an hour that contains sunrise
# averages in the dark minutes, while the irradiance-weighted mean is the cosine the hour's
# irradiance actually arrived at.

"""
Total solar irradiance at the mean Earth–Sun distance, W m⁻² (Kopp & Lean 2011,
*Geophys. Res. Lett.* 38, L01706).
"""
const SOLAR_CONSTANT = 1361.0

# Ceiling on the cloud optical thickness inverted from the clearness index. Transmissivity falls
# toward zero only slowly with τ, so a dark overcast step would otherwise invert to an arbitrarily
# large τ on a few W/m² of difference. Clouds thicker than this are rare, and the Gardner & Sharp
# (2010) albedo term that consumes τ is saturated well below it.
const CLOUD_OPTICAL_THICKNESS_MAX = 100.0

# Solar declination δ (degrees), equation of time E (minutes) and Earth–Sun distance R (au) at UTC
# time `t`, by the low-precision solar coordinates of Meeus (1991, *Astronomical Algorithms*,
# ch. 25) in the form of NOAA's Solar Calculator. Accurate to about 0.01° in δ and a few seconds in
# E for 1950–2050. (The shorter Fourier series of Spencer (1971) departs from these by up to 0.5°
# in δ near the equinoxes.)
function _solar_position(t::DateTime)
    T = (datetime2julian(t) - 2451545.0) / 36525                   # Julian centuries from J2000
    L₀ = mod(280.46646 + T * (36000.76983 + 0.0003032T), 360)      # geometric mean longitude
    M = 357.52911 + T * (35999.05029 - 0.0001537T)                 # mean anomaly
    e = 0.016708634 - T * (0.000042037 + 0.0000001267T)            # orbital eccentricity
    C = sind(M) * (1.914602 - T * (0.004817 + 0.000014T)) +
        sind(2M) * (0.019993 - 0.000101T) + sind(3M) * 0.000289    # equation of center
    Ω = 125.04 - 1934.136T
    λ = L₀ + C - 0.00569 - 0.00478sind(Ω)                          # apparent longitude
    ε = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - 0.001813T))) / 60) / 60 +
        0.00256cosd(Ω)                                             # corrected obliquity
    δ = asind(sind(ε) * sind(λ))
    y = tand(ε / 2)^2
    E = 4rad2deg(y * sind(2L₀) - 2e * sind(M) + 4e * y * sind(M) * cosd(2L₀) -
                 0.5y^2 * sind(4L₀) - 1.25e^2 * sind(2M))
    R = 1.000001018 * (1 - e^2) / (1 + e * cosd(M + C))
    return (; δ, E, R)
end

"""
    solar_declination(t::DateTime) -> Float64

Solar declination (degrees) at UTC time `t`, from the solar coordinates of Meeus (1991) as NOAA's
Solar Calculator implements them. Accurate to about 0.01° for 1950–2050.

# References
- Meeus, J. (1991). *Astronomical Algorithms*. Willmann-Bell, ch. 25.
- NOAA Global Monitoring Laboratory, Solar Calculator (https://gml.noaa.gov/grad/solcalc/).
"""
solar_declination(t::DateTime) = _solar_position(t).δ

"""
    equation_of_time(t::DateTime) -> Float64

Equation of time (minutes, apparent minus mean solar time) at UTC time `t`. See
[`solar_declination`](@ref) for the source.
"""
equation_of_time(t::DateTime) = _solar_position(t).E

"""
    earth_sun_distance_factor(t::DateTime) -> Float64

`(d₀/d)²`, the ratio of the top-of-atmosphere irradiance at UTC time `t` to that at the mean
Earth–Sun distance, with `d` from the Earth's orbit as in [`solar_declination`](@ref). Ranges
from about 0.967 in early July to 1.034 in early January.
"""
earth_sun_distance_factor(t::DateTime) = 1 / _solar_position(t).R^2

# Hour angle (radians, zero at local solar noon) at UTC time `t` and longitude `lon` (degrees
# east), with the equation of time `eot` (minutes) supplied by the caller.
function _hour_angle(t::DateTime, lon::Real, eot::Real)
    minutes = hour(t) * 60 + minute(t) + (second(t) + millisecond(t) / 1000) / 60
    return deg2rad((minutes + eot + 4lon) / 4 - 180)
end

"""
    cos_solar_zenith(latitude, longitude, t::DateTime) -> Float64

Instantaneous cosine of the solar zenith angle at UTC time `t`, for latitude and longitude in
degrees: `sin φ sin δ + cos φ cos δ cos ω`. Negative when the Sun is below the horizon.
"""
function cos_solar_zenith(latitude::Real, longitude::Real, t::DateTime)
    sun = _solar_position(t)
    ω = _hour_angle(t, longitude, sun.E)
    return sind(latitude) * sind(sun.δ) + cosd(latitude) * cosd(sun.δ) * cos(ω)
end

# ∫ μ dω and ∫ μ² dω over the sunlit part of the hour-angle interval [ω₁, ω₂], with
# μ(ω) = a + b·cos ω. The Sun is up where μ > 0: every ω for polar day, none for polar night, and
# otherwise the windows |ω − 2πk| < h₀ with cos h₀ = −a/b.
function _sunlit_moments(a::Float64, b::Float64, ω₁::Float64, ω₂::Float64)
    c = b > 0 ? -a / b : (a > 0 ? -Inf : Inf)
    c >= 1 && return (0.0, 0.0)
    h₀ = c <= -1 ? Float64(π) : acos(c)
    I₁ = I₂ = 0.0
    for k in floor(Int, (ω₁ - h₀) / 2π):ceil(Int, (ω₂ + h₀) / 2π)
        lo, hi = max(ω₁, 2π * k - h₀), min(ω₂, 2π * k + h₀)
        hi > lo || continue
        L, Δsin = hi - lo, sin(hi) - sin(lo)
        I₁ += a * L + b * Δsin
        I₂ += a^2 * L + 2a * b * Δsin + b^2 * (L / 2 + (sin(2hi) - sin(2lo)) / 4)
    end
    return (I₁, I₂)
end

"""
    solar_geometry(latitude, longitude, times; step = Hour(1)) -> NamedTuple

Interval-averaged solar geometry for forcing steps stamped `times` (UTC), each the mean over its
accumulation interval `(t − step, t]`, at `latitude` and `longitude` in degrees. Returns two
vectors on `times`' indices:

- `solar_zenith_angle` (degrees): `acos` of the irradiance-weighted mean `cos θ` over the sunlit
  part of the interval, `∫μ² dω / ∫μ dω` with `μ = cos θ ≥ 0`. A step with no sunlit part reports
  90 exactly.
- `shortwave_toa` (W m⁻²): mean top-of-atmosphere irradiance on a horizontal surface over the
  interval, `S₀·(d₀/d)²·⟨max(cos θ, 0)⟩`, with `S₀ =` [`SOLAR_CONSTANT`](@ref). Zero at night.

Both integrals are evaluated in closed form over the hour angle, so a step containing sunrise or
sunset is exact rather than sampled. Declination, the equation of time and `(d₀/d)²` are taken at
the interval's midpoint; their change over a step of a few hours is negligible.
`step` must be positive and at most a day.
"""
function solar_geometry(latitude::Real, longitude::Real, times::AbstractVector{DateTime};
                        step::Dates.Period = Hour(1))
    Δt = Dates.value(Millisecond(step)) / 3_600_000           # hours
    0 < Δt <= 24 || throw(ArgumentError("step must be positive and at most a day (got $(step))"))
    φ = deg2rad(Float64(latitude))
    Δω = deg2rad(15Δt)
    half = Millisecond(round(Int, Dates.value(Millisecond(step)) / 2))
    sza = similar(times, Float64)
    toa = similar(times, Float64)
    for i in eachindex(times)
        sun = _solar_position(times[i] - half)
        a, b = sin(φ) * sind(sun.δ), cos(φ) * cosd(sun.δ)
        ω₂ = _hour_angle(times[i], longitude, sun.E)
        I₁, I₂ = _sunlit_moments(a, b, ω₂ - Δω, ω₂)
        if I₁ > 0
            sza[i] = acosd(clamp(I₂ / I₁, 0.0, 1.0))
            toa[i] = SOLAR_CONSTANT * I₁ / (Δω * sun.R^2)
        else
            sza[i] = 90.0
            toa[i] = 0.0
        end
    end
    return (; solar_zenith_angle = sza, shortwave_toa = toa)
end

"""
    erbs_diffuse_fraction(clearness_index) -> Float64

Diffuse fraction of global horizontal shortwave from the clearness index `k_t`, by the
piecewise correlation of Erbs, Klein & Duffie (1982):

    f = 1 − 0.09 k_t                                                       k_t ≤ 0.22
    f = 0.9511 − 0.1604 k_t + 4.388 k_t² − 16.638 k_t³ + 12.336 k_t⁴       0.22 < k_t ≤ 0.80
    f = 0.165                                                              k_t > 0.80

Defined for `k_t ≥ 0`; a negative `k_t` throws.

# References
- Erbs, D. G., Klein, S. A. & Duffie, J. A. (1982). Estimation of the diffuse radiation fraction
  for hourly, daily and monthly-average global radiation. *Solar Energy* 28(4), 293–302.
"""
function erbs_diffuse_fraction(kt::Real)
    kt >= 0 || throw(ArgumentError("clearness index must be non-negative (got $(kt))"))
    kt <= 0.22 && return 1 - 0.09kt
    kt <= 0.80 && return 0.9511 - 0.1604kt + 4.388kt^2 - 16.638kt^3 + 12.336kt^4
    return 0.165
end

"""
    cloud_transmissivity(τ) -> Float64

Shortwave transmissivity of a non-absorbing cloud layer of optical thickness `τ`,
`1 / (1 + 0.13 τ)`: the complement of the cloud reflectivity `R = 0.13 τ / (1 + 0.13 τ)` of
Lacis & Hansen (1974), a two-stream result for conservative scattering with asymmetry parameter
0.85 (the coefficient is `√3(1 − g)/2`). Absorption within the cloud, and the extra downward
irradiance that multiple reflection between cloud base and a bright snow surface adds, are both
outside this relation.

# References
- Lacis, A. A. & Hansen, J. E. (1974). A parameterization for the absorption of solar radiation
  in the Earth's atmosphere. *J. Atmos. Sci.* 31, 118–133.
"""
cloud_transmissivity(τ::Real) = 1 / (1 + 0.13τ)

"""
    cloud_optical_thickness(transmissivity) -> Float64

Cloud optical thickness from a cloud transmissivity in `(0, 1]`, inverting
[`cloud_transmissivity`](@ref): `τ = (1/t − 1)/0.13`, capped at
`CLOUD_OPTICAL_THICKNESS_MAX` = $(CLOUD_OPTICAL_THICKNESS_MAX). `t = 1` (clear sky) gives `τ = 0`;
`t` outside `(0, 1]` throws.
"""
function cloud_optical_thickness(t::Real)
    0 < t <= 1 || throw(ArgumentError("cloud transmissivity must be in (0, 1] (got $(t))"))
    return min((1 / t - 1) / 0.13, CLOUD_OPTICAL_THICKNESS_MAX)
end

# Clear-sky clearness index at surface elevation `z` (m): the FAO-56 clear-sky solar radiation
# R_so = (0.75 + 2×10⁻⁵ z)·R_a (Allen et al. 1998, *FAO Irrigation and Drainage Paper 56*,
# Eq. 37), as a fraction of the top-of-atmosphere irradiance R_a.
_clear_sky_clearness(z::Real) = 0.75 + 2e-5 * z

"""
    _partition_shortwave(shortwave_downward, shortwave_toa, elevation) -> (diffuse, τ)

Diffuse shortwave (W m⁻²) and cloud optical thickness for one step from the step's downward
shortwave, its top-of-atmosphere irradiance and the surface elevation (m):

- clearness index `k_t = SW/TOA`, clamped to `[0, 1]`: at low sun, timing and refraction let
  reanalysis shortwave exceed the geometric top-of-atmosphere value, and such a step is treated
  as clear;
- diffuse `= f·SW` with `f` from [`erbs_diffuse_fraction`](@ref);
- cloud transmissivity `= min(k_t / k_clear, 1)`, `k_clear` the FAO-56 clear-sky clearness at the
  elevation, and `τ` from [`cloud_optical_thickness`](@ref); `k_t = 0` under daylight gives the
  cap.

With no top-of-atmosphere irradiance (night), the clearness index is undefined: all the
shortwave, normally zero, is taken as diffuse and `τ = 0`. `NaN` shortwave gives `NaN` for both.
"""
function _partition_shortwave(sw::Real, toa::Real, elevation::Real)
    isnan(sw) && return (NaN, NaN)
    toa > 0 || return (Float64(sw), 0.0)
    kt = clamp(sw / toa, 0.0, 1.0)
    transmissivity = min(kt / _clear_sky_clearness(elevation), 1.0)
    τ = transmissivity > 0 ? cloud_optical_thickness(transmissivity) : CLOUD_OPTICAL_THICKNESS_MAX
    return (erbs_diffuse_fraction(kt) * sw, τ)
end

"""
    shortwave_forcing_layers(shortwave_downward, latitude, longitude, times, elevation;
                             step = Hour(1)) -> NamedTuple

The three shortwave layers GEMB reads beside `shortwave_downward`, for steps stamped `times`
(UTC) whose fluxes are means over `(t − step, t]`, at a cell of `latitude`, `longitude` (degrees)
and surface `elevation` (m):

- `solar_zenith_angle` (degrees) from [`solar_geometry`](@ref);
- `shortwave_downward_diffuse` (W m⁻²) and `cloud_optical_thickness` (dimensionless) from the
  clearness index against that geometry's top-of-atmosphere irradiance — see
  [`_partition_shortwave`](@ref).

None of the three depends on air temperature or on elevation adjustments: the adjustments leave
`shortwave_downward` unchanged, so the diffuse part stays consistent with it.

An average of `solar_zenith_angle` over several cells must weight `cos θ` by each cell's
shortwave, and take `acos` of the result, to remain the angle the averaged irradiance arrives at.
"""
function shortwave_forcing_layers(shortwave_downward::AbstractVector, latitude::Real,
                                  longitude::Real, times::AbstractVector{DateTime},
                                  elevation::Real; step::Dates.Period = Hour(1))
    geometry = solar_geometry(latitude, longitude, times; step)
    diffuse = similar(shortwave_downward, Float64)
    τ = similar(shortwave_downward, Float64)
    for i in eachindex(shortwave_downward, geometry.shortwave_toa)
        diffuse[i], τ[i] = _partition_shortwave(shortwave_downward[i], geometry.shortwave_toa[i],
                                                elevation)
    end
    return (; geometry.solar_zenith_angle, shortwave_downward_diffuse = diffuse,
            cloud_optical_thickness = τ)
end
