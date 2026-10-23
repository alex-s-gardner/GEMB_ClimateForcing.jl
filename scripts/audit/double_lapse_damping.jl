# Finding 2 (review 2026-10-13): the fitted lapse rate was the on-glacier slope k·Γ, and the
# interval forcing applied it to a donor's temperature and then decoupled again, so the vertical
# term was damped twice. A band far above its donors came out too warm, more so the higher it sat.
#
# Synthetic region: 20 cells at 500–2500 m, glacier fraction rising with elevation, k = 0.8,
# ambient excess E(z) = 20 − 0.006 z (all cells above melting). The forcing at a band at 3000 m
# should be 273.15 + k·E(3000) = 274.75 K from every donor.
#
# Run (from the repository root):
#   julia --project=. scripts/audit/double_lapse_damping.jl                # this branch
#   julia --project=<checkout of 352d7db> scripts/audit/double_lapse_damping.jl
#
# Expected: at 352d7db the fitted slope is 4.8 K/km (k times the ambient 6) and the band error
# grows with height above the donor, from +0.01 K (donor at 2500 m) to +2.39 K (donor at 500 m);
# on this branch the slope is 6.0 K/km and every donor's error is below 1e-12 K.
include(joinpath(@__DIR__, "_common.jl"))

k, z_band = 0.8, 3000.0
E_amb(z) = 20 - 0.006z
z = collect(range(500.0, 2500.0; length = 20))
glm = region_glm(z)
T = T_REF .+ (1 .- glm .* (1 - k)) .* E_amb.(z)
acc = region_sums(z, glm, T)
d = G.derive_decoupling_factor(acc)
Γ = G.derive_lapse_rate(acc, d).lapse_rate[1]
k_band = G.decoupling_factor_at_elevation(d, z_band)[1]
truth = T_REF + k * E_amb(z_band)
println("fitted k(z̄) = ", round(d.decoupling_factor[1]; digits = 4), ", lapse rate = ",
        round(Γ; digits = 3), " K/km (ambient 6.0, k·ambient 4.8)")
for i in (1, 5, 10, 15, 20)
    fd = forcing_stack(; T = [T[i]], e = 0.6 * G.saturation_vapor_pressure(T[i]), elevation = z[i])
    band = if AMBIENT_FIRST
        k_donor = G.decoupling_factor_at_elevation(d, z[i])[1]
        G._cell_forcing_at_interval(fd, z_band - z[i], [k_band], [Γ], glm[i], [k_donor])
    else
        k_eff = 1 - (1 - k_band) * (1 - glm[i])   # GEMB_GlacierSims' per-cell weighting
        G._cell_forcing_at_interval(fd, z_band - z[i], [k_eff], [Γ])
    end
    err = band[:temperature_air][1] - truth
    println("donor z = ", round(Int, z[i]), " m, glm = ", round(glm[i]; digits = 2),
            ": band T − truth = ", round(err; sigdigits = 3), " K")
end
