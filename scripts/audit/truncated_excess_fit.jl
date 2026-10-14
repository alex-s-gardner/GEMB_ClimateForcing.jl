# Finding 11 (review 2026-10-13): the decoupling fit treated the truncated warm excess
# max(T − 273.15, 0) as linear in elevation, so where the melting level falls inside a region its
# cells at or below melting bent the fitted surface and biased k low; the lapse fit inherited it.
#
# Synthetic region: 20 cells at 500–2500 m, k = 0.8, ambient excess E(z) = 14 − 0.006 z, so the
# cells above 2333 m are below melting.
#
# Run (from the repository root):
#   julia --project=. scripts/audit/truncated_excess_fit.jl                # this branch
#   julia --project=<checkout of 352d7db> scripts/audit/truncated_excess_fit.jl
#
# Expected: at 352d7db k = 0.658 and the slope 3.98 K/km (its target was the on-glacier 4.80); on
# this branch k = 0.800 and the slope is the ambient 6.00 K/km.
include(joinpath(@__DIR__, "_common.jl"))

k = 0.8
E_amb(z) = 14 - 0.006z
z = collect(range(500.0, 2500.0; length = 20))
glm = region_glm(z)
T = [E_amb(zi) > 0 ? T_REF + (1 - g * (1 - k)) * E_amb(zi) : T_REF + E_amb(zi) for (zi, g) in zip(z, glm)]
acc = region_sums(z, glm, T)
d = G.derive_decoupling_factor(acc)
Γ = G.derive_lapse_rate(acc, d).lapse_rate[1]
println("cells below melting: ", count(<=(T_REF), T), " of ", length(T))
println("fitted k(z̄) = ", round(d.decoupling_factor[1]; digits = 3), " (truth 0.800)")
println("fitted lapse rate = ", round(Γ; digits = 2), " K/km (ambient 6.00, on-glacier 4.80)")
