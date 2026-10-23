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
# large τ on a few W/m² of difference. The Gardner & Sharp (2010) albedo parameterization that
# consumes τ is fitted over τ = 0–30 (Table B1, p. 12) and keeps rising with τ beyond that range,
# so the ceiling is the top of the fitted range.
const CLOUD_OPTICAL_THICKNESS_MAX = 30.0

# Largest solar zenith angle (degrees) at which τ is inverted from the clearness index. See
# `_partition_shortwave_series`.
const CLOUD_OPTICAL_THICKNESS_ZENITH_MAX = 80.0

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
    sin_M = sind(M)
    C = sin_M * (1.914602 - T * (0.004817 + 0.000014T)) +
        sind(2M) * (0.019993 - 0.000101T) + sind(3M) * 0.000289    # equation of center
    Ω = 125.04 - 1934.136T
    λ = L₀ + C - 0.00569 - 0.00478sind(Ω)                          # apparent longitude
    ε = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - 0.001813T))) / 60) / 60 +
        0.00256cosd(Ω)                                             # corrected obliquity
    δ = asind(sind(ε) * sind(λ))
    y = tand(ε / 2)^2
    E = 4rad2deg(y * sind(2L₀) - 2e * sin_M + 4e * y * sin_M * cosd(2L₀) -
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
    sin_φ, cos_φ = sin(φ), cos(φ)
    Δω = deg2rad(15Δt)
    half = Millisecond(round(Int, Dates.value(Millisecond(step)) / 2))
    sza = similar(times, Float64)
    toa = similar(times, Float64)
    for i in eachindex(times)
        sun = _solar_position(times[i] - half)
        a, b = sin_φ * sind(sun.δ), cos_φ * cosd(sun.δ)
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
`CLOUD_OPTICAL_THICKNESS_MAX` = $(CLOUD_OPTICAL_THICKNESS_MAX), the top of the τ range over which
Gardner & Sharp (2010, Table B1, p. 12) fitted the albedo parameterization that consumes it.
`t = 1` (clear sky) gives `τ = 0`; `t` outside `(0, 1]` throws.

# References
- Gardner, A. S. & Sharp, M. J. (2010). A review of snow and ice albedo and the development of a
  new physically based broadband albedo parameterization. *J. Geophys. Res.* 115, F01009.
"""
function cloud_optical_thickness(t::Real)
    0 < t <= 1 || throw(ArgumentError("cloud transmissivity must be in (0, 1] (got $(t))"))
    return min((1 / t - 1) / 0.13, CLOUD_OPTICAL_THICKNESS_MAX)
end

"""
    _clear_sky_clearness(cos_zenith, pressure, vapor_pressure) -> Float64

Clear-sky shortwave as a fraction of the top-of-atmosphere irradiance, `R_so/R_a = K_B + K_D`,
by FAO-56 Annex 3 eqs. 3-17 to 3-20 (Allen et al. 1998) for clean air (`K_t = 1`):

    K_B = 0.98·exp[−0.00146·P/sin φ − 0.091·(W/sin φ)^0.25]          (3-18)
    W   = 0.14·e_a·P + 2.1                                           (3-19)
    K_D = 0.35 − 0.33·K_B   for K_B ≥ 0.15,   0.18 + 0.82·K_B otherwise   (3-20)

with `sin φ = cos_zenith` the sine of the sun's elevation, `P` the surface pressure and `e_a` the
vapor pressure, both in kPa in the equations and passed here in Pa, and `W` the precipitable water
in mm. `K_B` is the direct-beam and `K_D` the diffuse clearness. A `cos_zenith` outside `(0, 1]`,
a non-positive pressure or a negative vapor pressure throws.

# References
- Allen, R. G., Pereira, L. S., Raes, D. & Smith, M. (1998). *Crop evapotranspiration:
  guidelines for computing crop water requirements.* FAO Irrigation and Drainage Paper 56,
  Annex 3, eqs. 3-17 to 3-20.
"""
function _clear_sky_clearness(cos_zenith::Real, pressure::Real, vapor_pressure::Real)
    0 < cos_zenith <= 1 ||
        throw(ArgumentError("cosine of the solar zenith angle must be in (0, 1] (got $(cos_zenith))"))
    pressure > 0 || throw(ArgumentError("surface pressure must be positive (got $(pressure) Pa)"))
    vapor_pressure >= 0 ||
        throw(ArgumentError("vapor pressure must be non-negative (got $(vapor_pressure) Pa)"))
    P, e_a = pressure / 1000, vapor_pressure / 1000                       # kPa
    W = 0.14 * e_a * P + 2.1                                              # mm
    K_B = 0.98 * exp(-0.00146 * P / cos_zenith - 0.091 * (W / cos_zenith)^0.25)
    K_D = K_B >= 0.15 ? 0.35 - 0.33 * K_B : 0.18 + 0.82 * K_B
    return K_B + K_D
end

# Clearness index `k_t = SW/TOA` of a sunlit step, clamped to `[0, 1]`: at low sun, timing and
# refraction let reanalysis shortwave exceed the geometric top-of-atmosphere value, and such a step
# is treated as clear.
_clearness_index(sw::Real, toa::Real) = clamp(sw / toa, 0.0, 1.0)

# Diffuse shortwave (W m⁻²) of one step: `f·SW` with `f` from `erbs_diffuse_fraction(k_t)`. With
# no top-of-atmosphere irradiance (night) the clearness index is undefined and all the shortwave,
# normally zero, is taken as diffuse. `NaN` shortwave gives `NaN`.
function _diffuse_shortwave(sw::Real, toa::Real)
    isnan(sw) && return NaN
    toa > 0 || return Float64(sw)
    return erbs_diffuse_fraction(_clearness_index(sw, toa)) * sw
end

"""
    _inverted_optical_thickness(shortwave_downward, shortwave_toa, cos_zenith, pressure,
                                vapor_pressure) -> Float64

Cloud optical thickness of one sunlit step (`shortwave_toa > 0`): the cloud transmissivity
`min(k_t / k_clear, 1)`, with `k_t` the clamped clearness index `SW/TOA` and `k_clear` the FAO-56
clear-sky clearness of [`_clear_sky_clearness`](@ref) at the step's sun angle, surface pressure
and vapor pressure (Pa), inverted by [`cloud_optical_thickness`](@ref). `k_t = 0` gives
`CLOUD_OPTICAL_THICKNESS_MAX`. `NaN` shortwave gives `NaN`; a step with no top-of-atmosphere
irradiance throws.
"""
function _inverted_optical_thickness(sw::Real, toa::Real, cos_zenith::Real, pressure::Real,
                                     vapor_pressure::Real)
    isnan(sw) && return NaN
    toa > 0 || throw(ArgumentError("τ is inverted only on a sunlit step (got TOA = $(toa) W/m²)"))
    k_clear = _clear_sky_clearness(cos_zenith, pressure, vapor_pressure)
    transmissivity = min(_clearness_index(sw, toa) / k_clear, 1.0)
    return transmissivity > 0 ? cloud_optical_thickness(transmissivity) : CLOUD_OPTICAL_THICKNESS_MAX
end

"""
    _partition_shortwave_series(shortwave_downward, shortwave_toa, solar_zenith_angle,
                                pressure_air, vapor_pressure) -> (diffuse, τ)

Diffuse shortwave (W m⁻²) and cloud optical thickness on a time series of steps, from each step's
downward shortwave, top-of-atmosphere irradiance, solar zenith angle (degrees), surface pressure
and vapor pressure (Pa). All five must share their indices; the results do too.

- Diffuse: `f·SW` with `f` from [`erbs_diffuse_fraction`](@ref) of the clearness index
  `k_t = SW/TOA` (clamped to `[0, 1]`), which does not involve the clear-sky model; at night
  (`TOA = 0`) all the shortwave is diffuse.
- `τ` at night is 0.
- `τ` at a sunlit step with zenith angle at most `CLOUD_OPTICAL_THICKNESS_ZENITH_MAX` =
  $(CLOUD_OPTICAL_THICKNESS_ZENITH_MAX)° is inverted from the step's own clearness index against
  the FAO-56 Annex 3 clear sky at `sin φ = cos θ` — see [`_inverted_optical_thickness`](@ref).
- `τ` at a sunlit step with a larger zenith angle is not inverted. It takes the τ of the latest
  earlier inverted step in the same daylight period, or, where the period has none yet (dawn),
  of the first later one; a period with no step at or below the threshold (polar twilight) gets
  `τ = 0`. A daylight period is a maximal run of consecutive steps with `TOA > 0`, so under polar
  day the run spans several days and the low-sun hours around midnight carry the previous
  evening's τ.

The inversion is restricted to a sun well above the horizon for three reasons. The FAO-56 clear
sky takes the atmospheric path as the plane-parallel air mass `1/sin φ` (Annex 3, eq. 3-18),
which grows without bound as the sun approaches the horizon, so the clear-sky reference is least
reliable there. `τ` is ill-conditioned when `SW` and `TOA` are both small: a few W m⁻² of timing
or refraction error moves it across its whole range. And the Gardner & Sharp (2010) albedo
parameterization that consumes `τ` is fitted for zenith angles of at most 85° (Table B1, p. 12).

`NaN` shortwave gives `NaN` for both at that step, and a low-sun step that takes its τ from a
`NaN` step is `NaN`.

# References
- Allen, R. G., Pereira, L. S., Raes, D. & Smith, M. (1998). FAO Irrigation and Drainage Paper
  56, Annex 3.
- Gardner, A. S. & Sharp, M. J. (2010). *J. Geophys. Res.* 115, F01009.
"""
function _partition_shortwave_series(sw::AbstractVector, toa::AbstractVector,
                                     sza::AbstractVector, pressure::AbstractVector,
                                     vapor_pressure::AbstractVector)
    indices = eachindex(sw, toa, sza, pressure, vapor_pressure)
    diffuse = similar(sw, Float64)
    τ = similar(sw, Float64)
    held = NaN                       # τ of the latest inverted step in the current period
    have_held = false
    dawn = eltype(indices)[]         # low-sun steps of the current period before its first inversion
    for i in indices
        diffuse[i] = _diffuse_shortwave(sw[i], toa[i])
        if !(toa[i] > 0)
            τ[i] = 0.0                                  # night ends the daylight period
            have_held = false
            empty!(dawn)                                # a period with no inversion keeps τ = 0
        elseif sza[i] <= CLOUD_OPTICAL_THICKNESS_ZENITH_MAX
            τ[i] = _inverted_optical_thickness(sw[i], toa[i], cosd(sza[i]), pressure[i],
                                               vapor_pressure[i])
            for j in dawn
                τ[j] = τ[i]
            end
            empty!(dawn)
            held, have_held = τ[i], true
        elseif isnan(sw[i])
            τ[i] = NaN
        elseif have_held
            τ[i] = held
        else
            τ[i] = 0.0                                  # replaced if the period has an inversion
            push!(dawn, i)
        end
    end
    return (diffuse, τ)
end

"""
    shortwave_forcing_layers(shortwave_downward, pressure_air, vapor_pressure, latitude,
                             longitude, times; step = Hour(1)) -> NamedTuple

The three shortwave layers GEMB reads beside `shortwave_downward`, for steps stamped `times`
(UTC) whose fluxes are means over `(t − step, t]`, at a cell of `latitude`, `longitude`
(degrees), from the cell's downward shortwave (W m⁻²), surface pressure and vapor pressure (Pa)
on those steps:

- `solar_zenith_angle` (degrees) from [`solar_geometry`](@ref);
- `shortwave_downward_diffuse` (W m⁻²) and `cloud_optical_thickness` (dimensionless) from the
  clearness index against that geometry's top-of-atmosphere irradiance, with τ inverted against
  the FAO-56 Annex 3 clear sky at the step's sun angle, pressure and vapor pressure, and held
  across low-sun steps — see [`_partition_shortwave_series`](@ref).

None of the three depends on air temperature, and the forcing adjustments do not re-derive them:
the adjustments leave `shortwave_downward` unchanged, so the diffuse part stays consistent with
it, and τ remains the cloud inferred from that shortwave at the pressure and vapor pressure it
was derived with.

An average of `solar_zenith_angle` over several cells must weight `cos θ` by each cell's
shortwave, and take `acos` of the result, to remain the angle the averaged irradiance arrives at.
"""
function shortwave_forcing_layers(shortwave_downward::AbstractVector, pressure_air::AbstractVector,
                                  vapor_pressure::AbstractVector, latitude::Real, longitude::Real,
                                  times::AbstractVector{DateTime}; step::Dates.Period = Hour(1))
    geometry = solar_geometry(latitude, longitude, times; step)
    diffuse, τ = _partition_shortwave_series(shortwave_downward, geometry.shortwave_toa,
                                             geometry.solar_zenith_angle, pressure_air,
                                             vapor_pressure)
    return (; geometry.solar_zenith_angle, shortwave_downward_diffuse = diffuse,
            cloud_optical_thickness = τ)
end
