# Finding 10 (review 2026-10-13): ΔT = 0 and Δz = 0 were not identities. Vapor pressure came from
# the dewpoint by one formula (Magnus over water), relative humidity was formed against another
# (Buck over ice below 0 °C), and the adjustments clamped RH to 1. Saturated cold air, and warm air
# made near-saturated by the glacier decoupling, therefore lost vapor and longwave even under an
# identity adjustment.
#
# Run (from the repository root):
#   julia --project=. scripts/audit/humidity_identity.jl                   # this branch
#   julia --project=<checkout of 352d7db> scripts/audit/humidity_identity.jl
#
# Expected: at 352d7db the dewpoint-saturated steps report RH 1.22, 1.10 and 1.05 at −20, −10 and
# −5 °C, and lose 20–27 Pa and 0.9–2.5 W/m² under Δz = 0 and ΔT = 0; the decoupled warm step
# loses 316 Pa and 6.9 W/m². On this branch RH(Td = T) = 1 and every Δ is exactly 0.
include(joinpath(@__DIR__, "_common.jl"))

T = [253.15, 263.15, 268.15, 273.65, 283.15, 293.15]
e = G.dewpoint_to_vapor_pressure(T)                      # dewpoint = air temperature: saturated
println("e(Td = T) / eₛ(T) used for RH: ", round.(e ./ G.saturation_vapor_pressure.(T); digits = 3))
s = forcing_stack(; T, e, LW = 260.0)
for (name, adj) in (("climate_adjust_for_elevation(s, 0)", G.climate_adjust_for_elevation(s, 0.0)),
                    ("temperature_adjust(s, 0)", G.temperature_adjust(s, 0.0)))
    println(rpad(name, 36), "Δe = ", round.(parent(adj[:vapor_pressure]) .- e; digits = 2),
            " Pa,  ΔLW = ", round.(parent(adj[:longwave_downward]) .- 260; digits = 3), " W/m²")
end
# Warm humid air decoupled with k = 0.6, then a ΔT = 0 perturbation.
T2 = fill(285.15, 3)
s2 = forcing_stack(; T = T2, e = 0.95 .* G.saturation_vapor_pressure.(T2), LW = 330.0)
dec = G.climate_adjust_for_glacier(s2, 0.6)
again = G.temperature_adjust(dec, 0.0)
println("after decoupling, ΔT = 0:            Δe = ",
        round(parent(again[:vapor_pressure])[1] - parent(dec[:vapor_pressure])[1]; digits = 2),
        " Pa,  ΔLW = ",
        round(parent(again[:longwave_downward])[1] - parent(dec[:longwave_downward])[1]; digits = 3),
        " W/m²")
