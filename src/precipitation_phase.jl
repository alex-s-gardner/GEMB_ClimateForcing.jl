"""
    JENNINGS2018_BIVARIATE

Coefficients `(α, β, γ)` of the bivariate binary logistic precipitation-phase model of Jennings et
al. (2018), optimized on 17.8 million Northern Hemisphere station observations (1978–2007):

    p(snow) = 1 / (1 + exp(α + β·Tₛ + γ·RH))

with `Tₛ` the air temperature [°C] and `RH` the relative humidity [%]. Values from their
Supplementary Table 2: α = −10.04 (0.68), β = 1.41 (0.05), γ = 0.09 (0.007), standard deviations in
parentheses. The probability is 0.5 at 0.74 °C at 100 % relative humidity and at 3.93 °C at 50 %.

# References
- Jennings, K. S., Winchell, T. S., Livneh, B. and Molotch, N. P. (2018). Spatial variation of the
  rain–snow temperature threshold across the Northern Hemisphere. *Nat. Commun.* 9, 1148.
  doi:10.1038/s41467-018-03629-7
"""
const JENNINGS2018_BIVARIATE = (α = -10.04, β = 1.41, γ = 0.09)

"""
    snowfall_fraction(temperature_air, vapor_pressure) -> Float64

Fraction [-] of precipitation falling as snow at air temperature `temperature_air` [K] and vapor
pressure `vapor_pressure` [Pa], from the bivariate logistic model of Jennings et al. (2018)
([`JENNINGS2018_BIVARIATE`](@ref)). The model gives the probability of snow; it is used here as the
snow fraction of the step's precipitation, which spreads the transition over the several degrees
the observations show rather than switching at one temperature.

Relative humidity is taken with respect to water at every temperature
([`saturation_vapor_pressure`](@ref)), the convention of the station humidity the model was fitted
to, and is capped at 100 %.

Accepts scalars; broadcast over arrays.
"""
function snowfall_fraction(temperature_air::Real, vapor_pressure::Real)
    (; α, β, γ) = JENNINGS2018_BIVARIATE
    T = Float64(temperature_air)
    rh = min(100 * Float64(vapor_pressure) / saturation_vapor_pressure(T), 100.0)
    return 1 / (1 + exp(α + β * (T - 273.15) + γ * rh))
end

"""
    snowfall_fraction_layer(stack::DimStack) -> DimArray

[`snowfall_fraction`](@ref) at every step of a forcing stack, from its `temperature_air` and
`vapor_pressure`, on the stack's time dimension.
"""
snowfall_fraction_layer(stack::DimStack) =
    snowfall_fraction.(stack[:temperature_air], stack[:vapor_pressure])

"""
    with_snowfall_fraction(stack::DimStack) -> DimStack

`stack` with a `snowfall_fraction` layer ([`snowfall_fraction_layer`](@ref)) added, or replaced, which
GEMB partitions each step's precipitation into snow and rain by instead of a temperature threshold.
"""
function with_snowfall_fraction(stack::DimStack)
    layers = merge(NamedTuple(stack), (; snowfall_fraction = snowfall_fraction_layer(stack)))
    return DimStack(layers; metadata = metadata(stack))
end
