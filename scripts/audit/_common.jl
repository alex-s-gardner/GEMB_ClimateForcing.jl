# Shared by the audit scripts: a one-cell forcing stack and the cross-cell sums the downscaling
# fits read, written so the same script runs against both versions of the package. The sums carry
# the field names of both accumulator layouts (all-cell `g, w, zg, …` and the `_warm` set); each
# version reads the fields it knows.
using GEMB_ClimateForcing, DimensionalData, Dates, Statistics
const G = GEMB_ClimateForcing
const T_REF = 273.15

# A forcing stack of `n` steps from per-layer values (scalars are filled).
function forcing_stack(; T, e = 300.0, P = 80000.0, LW = 280.0, SW = 500.0, n = length(T),
                       elevation = 1000.0, start = DateTime(2020, 7, 1, 12))
    ti = Ti(start .+ Hour.(0:(n - 1)))
    layer(v) = DimArray(v isa Real ? fill(Float64(v), n) : Float64.(v), ti)
    return DimStack((temperature_air = layer(T), pressure_air = layer(P), vapor_pressure = layer(e),
                     wind_speed = layer(3.0), precipitation = layer(0.0),
                     shortwave_downward = layer(SW), longwave_downward = layer(LW));
                    metadata = Dict{String,Any}("elevation" => elevation))
end

# One-timestep accumulator over cells with elevations `z`, glacier fractions `glm` and
# temperatures `T`, carrying both layouts' fields.
function region_sums(z, glm, T)
    E = max.(T .- T_REF, 0.0)
    w = glm .* z
    warm = E .> 0
    s(v) = [sum(v)]
    sw(v) = [sum(v[warm])]
    S = (s1 = [length(z)], z = s(z), zz = s(z .^ 2), T = s(T), zT = s(z .* T),
         g = s(glm), w = s(w), zg = s(z .* glm), zw = s(z .* w), gg = s(glm .^ 2), gw = s(glm .* w),
         ww = s(w .^ 2),
         s1_warm = [count(warm)], z_warm = sw(z), g_warm = sw(glm), w_warm = sw(w),
         zz_warm = sw(z .^ 2), zg_warm = sw(z .* glm), zw_warm = sw(z .* w), gg_warm = sw(glm .^ 2),
         gw_warm = sw(glm .* w), ww_warm = sw(w .^ 2),
         E = s(E), zE = s(z .* E), gE = s(glm .* E), wE = s(w .* E), EE = s(E .^ 2))
    return (n = 1, time = [DateTime(2020, 7, 1, 12)], mean_elevation = mean(z), sums = S)
end

# Glacier fraction rising with elevation, plus a fixed perturbation.
region_glm(z) = [clamp(clamp((zi - 500) / 2000, 0, 1)^2 + 0.05 * sin(3i), 0, 1) for (i, zi) in enumerate(z)]

# Which version is loaded: the six-argument interval forcing exists only on fix/forcing-audit.
const AMBIENT_FIRST = hasmethod(G._cell_forcing_at_interval, NTuple{6,Any})
println("GEMB_ClimateForcing at ", pkgdir(GEMB_ClimateForcing), AMBIENT_FIRST ?
        " (ambient-lapse API)" : " (original API)")
