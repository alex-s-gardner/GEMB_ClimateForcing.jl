# Latent finding (review 2026-10-13): validate_climate_forcing_units passed a stack containing NaN,
# because `minimum`/`maximum` propagate NaN and every range comparison against NaN is false.
#
# Run (from the repository root):
#   julia --project=. scripts/audit/nan_validation.jl                      # this branch
#   julia --project=<checkout of 352d7db> scripts/audit/nan_validation.jl
#
# Expected: at 352d7db every layer reports "passed"; on this branch every layer reports an
# ArgumentError naming it as non-finite.
include(joinpath(@__DIR__, "_common.jl"))

for name in (:temperature_air, :pressure_air, :vapor_pressure, :wind_speed, :precipitation,
             :shortwave_downward, :longwave_downward)
    s = forcing_stack(; T = [270.0, 271.0, 272.0], e = 400.0, SW = 100.0)
    v = copy(parent(s[name])); v[2] = NaN
    bad = rebuild(s; data = merge(map(parent, NamedTuple(s)), (; name => v)))
    result = try
        G.validate_climate_forcing_units(bad)
        "passed"
    catch err
        err isa ArgumentError ? "ArgumentError: " * first(split(err.msg, '\n'; limit = 3)[2:end]) :
                                rethrow()
    end
    println(rpad(name, 20), result)
end
