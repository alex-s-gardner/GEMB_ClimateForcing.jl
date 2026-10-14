# Finding 18 (review 2026-10-13): the glacier decoupling recomputed downward longwave from the
# cooled 2 m temperature, about −4 W/m² per K of cooling, although the decoupling factor is
# measured from 2 m temperature alone and ERA5-Land longwave does not depend on the glacier mask.
#
# Run (from the repository root):
#   julia --project=. scripts/audit/decoupling_longwave.jl                 # this branch
#   julia --project=<checkout of 352d7db> scripts/audit/decoupling_longwave.jl
#
# Expected: at 352d7db ΔLW is −4.2 to −4.6 W/m² per K of cooling and Δe = 0; on this branch
# ΔLW = 0 and, the input being at 99 % relative humidity, vapor pressure falls to saturation of
# the cooled air (−32 to −673 Pa).
include(joinpath(@__DIR__, "_common.jl"))

T = [275.15, 280.15, 285.15, 290.15]
s = forcing_stack(; T, e = 0.99 .* G.saturation_vapor_pressure.(T), LW = 320.0)
out = G.climate_adjust_for_glacier(s, 0.6)
ΔT = parent(out[:temperature_air]) .- T
ΔLW = parent(out[:longwave_downward]) .- 320.0
Δe = parent(out[:vapor_pressure]) .- parent(s[:vapor_pressure])
println("ΔT  = ", round.(ΔT; digits = 2), " K")
println("ΔLW = ", round.(ΔLW; digits = 2), " W/m²  (per K of cooling: ",
        round.(ΔLW ./ ΔT; digits = 2), ")")
println("Δe  = ", round.(Δe; digits = 1), " Pa")
