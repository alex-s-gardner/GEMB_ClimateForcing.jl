#!/usr/bin/env julia

"""
Analytic (offline) tests for direct temperature / precipitation perturbations of
climate forcing (`temperature_adjust`, `precipitation_adjust`).
No network access required.
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData
using Statistics

using GEMB_ClimateForcing: saturation_vapor_pressure,
    konzelmann_clear_sky_emissivity, _SIGMA_SB

# Physically reasonable high-elevation cold-surface stack. Temperature varies by
# month so the perturbations are exercised on a non-constant series.
function make_perturb_stack(; n_months::Int=12)
    times = [DateTime(2020, m, 15, 12) for m in 1:n_months]
    time_dim = Ti(times)
    npts = length(times)
    T  = [255.0 + 10.0 * sinpi((m - 1) / 12) for m in 1:npts]  # K
    P  = fill(80000.0, npts)            # Pa
    e  = 0.7 .* saturation_vapor_pressure.(T)                  # Pa, 70 % RH
    ws = fill(6.0, npts)                # m/s
    pr = fill(0.2, npts)                # kg/m²/hr
    sw = fill(150.0, npts)              # W/m²
    lw = fill(220.0, npts)              # W/m²
    return DimStack((
        temperature_air = DimArray(T, (time_dim,); metadata=Dict("units" => "K")),
        pressure_air = DimArray(P, (time_dim,); metadata=Dict("units" => "Pa")),
        vapor_pressure = DimArray(e, (time_dim,); metadata=Dict("units" => "Pa")),
        wind_speed = DimArray(ws, (time_dim,); metadata=Dict("units" => "m/s")),
        precipitation = DimArray(pr, (time_dim,); metadata=Dict("units" => "kg/m²/hr")),
        shortwave_downward = DimArray(sw, (time_dim,); metadata=Dict("units" => "W/m²")),
        longwave_downward = DimArray(lw, (time_dim,); metadata=Dict("units" => "W/m²")),
    ); metadata=Dict("elevation" => 1500.0, "dataset" => "test"))
end

const _FORCING_VARS = (:temperature_air, :pressure_air, :vapor_pressure, :wind_speed,
    :precipitation, :shortwave_downward, :longwave_downward)

@testset "Climate Adjustment (temperature / precipitation)" begin

    @testset "Identity perturbations" begin
        stack = make_perturb_stack()

        # ΔT = 0 reproduces every variable bit for bit.
        adj = temperature_adjust(stack, 0.0)
        for var in _FORCING_VARS
            @test parent(adj[var]) == parent(stack[var])
        end

        # scaling = 1 reproduces every variable exactly.
        adj = precipitation_adjust(stack, 1.0)
        for var in _FORCING_VARS
            @test all(parent(adj[var]) .≈ parent(stack[var]))
        end
    end

    @testset "Identity on saturated input" begin
        # Saturated with respect to the dewpoint's own curve (dewpoint = air temperature), cold
        # enough to be supersaturated over ice. The identity adjustments must return vapor pressure
        # and longwave bit for bit — nothing in them may clamp humidity.
        T = 253.15
        e = dewpoint_to_vapor_pressure(T)
        ti = Ti([DateTime(2020, 1, 15, 12)])
        layer(v) = DimArray([v], (ti,))
        s = DimStack((temperature_air = layer(T), pressure_air = layer(70000.0),
                      vapor_pressure = layer(e), wind_speed = layer(5.0),
                      precipitation = layer(0.1), shortwave_downward = layer(0.0),
                      longwave_downward = layer(200.0)); metadata = Dict("elevation" => 2000.0))
        for adj in (climate_adjust_for_elevation(s, 0.0), temperature_adjust(s, 0.0))
            @test parent(adj[:vapor_pressure]) == [e]
            @test parent(adj[:longwave_downward]) == [200.0]
        end
    end

    @testset "ΔT = 0 after decoupling is the identity" begin
        # Decoupling cools the air at fixed vapor pressure, raising relative humidity; a following
        # ΔT = 0 must still change nothing.
        n = 6
        ti = Ti(DateTime(2020, 7, 1) .+ Hour.(0:(n - 1)))
        layer(v) = DimArray(fill(v, n), (ti,))
        T = 285.15
        s = DimStack((temperature_air = layer(T), pressure_air = layer(85000.0),
                      vapor_pressure = layer(0.95 * saturation_vapor_pressure(T)),
                      wind_speed = layer(3.0), precipitation = layer(0.0),
                      shortwave_downward = layer(400.0), longwave_downward = layer(330.0));
                     metadata = Dict("elevation" => 1500.0))
        decoupled = climate_adjust_for_glacier(s, 0.6)
        again = temperature_adjust(decoupled, 0.0)
        for var in _FORCING_VARS
            @test parent(again[var]) == parent(decoupled[var])
        end
    end

    @testset "Decoupling caps vapor at saturation of the cooled air" begin
        # Decoupling cools air toward 273.15 K: T′ = 273.15 + k·(T − 273.15). Water-saturated input
        # cooled to T′ is capped at the carry curve's saturation, over ice below 273.16 K and over
        # water above. IFS eq. 7.5 constants written out independently.
        es_w(T) = 611.21 * exp(17.502 * (T - 273.16) / (T - 32.19))
        es_i(T) = 611.21 * exp(22.587 * (T - 273.16) / (T + 0.7))
        ti = Ti([DateTime(2020, 7, 1, 12)])
        layer(v) = DimArray([v], (ti,))
        stack_at(T) = DimStack((temperature_air = layer(T), pressure_air = layer(85000.0),
                                vapor_pressure = layer(es_w(T)), wind_speed = layer(3.0),
                                precipitation = layer(0.0), shortwave_downward = layer(400.0),
                                longwave_downward = layer(330.0)); metadata = Dict("elevation" => 1500.0))
        # T = 274.15 K, k = 0.005 ⇒ T′ = 273.155 K, below 273.16 K: capped at ice saturation.
        cold = climate_adjust_for_glacier(stack_at(274.15), 0.005)
        T′ = parent(cold[:temperature_air])[1]
        @test T′ ≈ 273.155 rtol = 1e-12
        @test parent(cold[:vapor_pressure])[1] ≈ es_i(T′) rtol = 1e-12
        # T = 283.15 K, k = 0.6 ⇒ T′ = 279.15 K: capped at water saturation.
        warm = climate_adjust_for_glacier(stack_at(283.15), 0.6)
        T′w = parent(warm[:temperature_air])[1]
        @test T′w ≈ 279.15 rtol = 1e-12
        @test parent(warm[:vapor_pressure])[1] ≈ es_w(T′w) rtol = 1e-12
    end

    @testset "Temperature offset is exact" begin
        stack = make_perturb_stack()
        for ΔT in (-5.0, 0.5, 2.0, 8.0)
            adj = temperature_adjust(stack, ΔT)
            @test all(parent(adj[:temperature_air]) .≈ parent(stack[:temperature_air]) .+ ΔT)
        end
        # Integer argument is accepted and promoted.
        adj = temperature_adjust(stack, 2)
        @test all(parent(adj[:temperature_air]) .≈ parent(stack[:temperature_air]) .+ 2.0)
    end

    @testset "temperature_adjust leaves independent variables untouched" begin
        stack = make_perturb_stack()
        adj = temperature_adjust(stack, 3.0)
        for var in (:pressure_air, :wind_speed, :precipitation, :shortwave_downward)
            @test parent(adj[var]) == parent(stack[var])
        end
    end

    @testset "Relative humidity conserved under warming" begin
        stack = make_perturb_stack()
        T = parent(stack[:temperature_air])
        e = parent(stack[:vapor_pressure])
        # Relative humidity with respect to ice below 273.16 K and water above.
        es_ref = GEMB_ClimateForcing.saturation_vapor_pressure_reference
        RH = e ./ es_ref.(T)

        for ΔT in (-4.0, 2.0, 6.0)
            adj = temperature_adjust(stack, ΔT)
            T′ = parent(adj[:temperature_air])
            e′ = parent(adj[:vapor_pressure])
            RH′ = e′ ./ es_ref.(T′)
            @test all(RH′ .≈ RH)
        end
    end

    @testset "Vapor pressure follows the saturation curve" begin
        stack = make_perturb_stack()
        e = parent(stack[:vapor_pressure])
        @test all(parent(temperature_adjust(stack, 3.0)[:vapor_pressure]) .> e)
        @test all(parent(temperature_adjust(stack, -3.0)[:vapor_pressure]) .< e)
    end

    @testset "Longwave response and preserved cloud increment" begin
        stack = make_perturb_stack()
        T = parent(stack[:temperature_air])
        e = parent(stack[:vapor_pressure])
        LW = parent(stack[:longwave_downward])

        # Cloud/aerosol emissivity increment diagnosed from the input.
        Δε = LW ./ (_SIGMA_SB .* T .^ 4) .- konzelmann_clear_sky_emissivity.(e, T)

        ΔT = 2.0
        adj = temperature_adjust(stack, ΔT)
        T′ = parent(adj[:temperature_air])
        e′ = parent(adj[:vapor_pressure])
        LW′ = parent(adj[:longwave_downward])

        # Warming increases downwelling longwave.
        @test all(LW′ .> LW)

        # Independently recomputed increment matches the input increment.
        Δε′ = LW′ ./ (_SIGMA_SB .* T′ .^ 4) .- konzelmann_clear_sky_emissivity.(e′, T′)
        @test all(Δε′ .≈ Δε)

        # Cooling decreases it.
        @test all(parent(temperature_adjust(stack, -2.0)[:longwave_downward]) .< LW)
    end

    @testset "Precipitation scaling is exact" begin
        stack = make_perturb_stack()
        pr = parent(stack[:precipitation])
        for s in (0.0, 0.85, 1.15, 3.0)
            adj = precipitation_adjust(stack, s)
            @test all(parent(adj[:precipitation]) .≈ pr .* s)
        end
        @test all(parent(precipitation_adjust(stack, 0.0)[:precipitation]) .== 0.0)
    end

    @testset "precipitation_adjust leaves all other variables untouched" begin
        stack = make_perturb_stack()
        adj = precipitation_adjust(stack, 1.4)
        for var in _FORCING_VARS
            var === :precipitation && continue
            @test parent(adj[var]) == parent(stack[var])
        end
    end

    @testset "Invalid scaling rejected" begin
        stack = make_perturb_stack()
        @test_throws ArgumentError precipitation_adjust(stack, -0.1)
        @test_throws ArgumentError precipitation_adjust(stack, -1.0)
    end

    @testset "Layer metadata preserved" begin
        stack = make_perturb_stack()
        adj = precipitation_adjust(temperature_adjust(stack, 2.0), 1.2)
        for var in _FORCING_VARS
            @test metadata(adj[var])["units"] == metadata(stack[var])["units"]
        end
    end

    @testset "Metadata bookkeeping" begin
        stack = make_perturb_stack()

        adj = temperature_adjust(stack, 2.0)
        @test metadata(adj)["delta_temperature"] == 2.0
        @test metadata(adj)["temperature_offset"] == 2.0
        @test metadata(adj)["temperature_air_mean"] ≈
              mean(parent(stack[:temperature_air])) + 2.0
        # Untouched source keys survive.
        @test metadata(adj)["elevation"] == 1500.0
        @test metadata(adj)["dataset"] == "test"

        # Temperature offsets accumulate additively.
        adj2 = temperature_adjust(adj, 1.5)
        @test metadata(adj2)["delta_temperature"] == 1.5
        @test metadata(adj2)["temperature_offset"] == 3.5

        # Precipitation scalings accumulate multiplicatively.
        padj = precipitation_adjust(stack, 1.1)
        @test metadata(padj)["precipitation_scaling"] ≈ 1.1
        @test metadata(padj)["precipitation_mean"] ≈
              mean(parent(stack[:precipitation])) * 1.1 * 8760.0
        padj2 = precipitation_adjust(padj, 1.1)
        @test metadata(padj2)["precipitation_scaling"] ≈ 1.21

        # Neither function touches the elevation bookkeeping.
        @test !haskey(metadata(adj2), "elevation_offset")
        @test !haskey(metadata(padj2), "elevation_offset")
        @test metadata(padj2)["elevation"] == 1500.0
    end

    @testset "Composition of the two perturbations" begin
        stack = make_perturb_stack()
        combined = precipitation_adjust(temperature_adjust(stack, 2.0), 0.85)

        @test all(parent(combined[:temperature_air]) .≈
                  parent(stack[:temperature_air]) .+ 2.0)
        @test all(parent(combined[:precipitation]) .≈ parent(stack[:precipitation]) .* 0.85)
        @test metadata(combined)["temperature_offset"] == 2.0
        @test metadata(combined)["precipitation_scaling"] ≈ 0.85

        # Order does not matter — the two are independent.
        reversed = temperature_adjust(precipitation_adjust(stack, 0.85), 2.0)
        for var in _FORCING_VARS
            @test all(parent(reversed[var]) .≈ parent(combined[var]))
        end
    end

    @testset "Composition with elevation adjustment" begin
        stack = make_perturb_stack()
        adj = temperature_adjust(climate_adjust_for_elevation(stack, 200.0), 2.0)
        @test metadata(adj)["elevation_offset"] == 200.0
        @test metadata(adj)["temperature_offset"] == 2.0
        @test metadata(adj)["elevation"] == 1700.0
    end

    @testset "Out-of-range perturbation fails validation" begin
        stack = make_perturb_stack()
        # Pushes temperature below the 180 K floor.
        @test_throws "temperature_air" temperature_adjust(stack, -100.0)
    end

    @testset "A precipitation scaling is not a unit error" begin
        stack = make_perturb_stack()
        # The 100 kg/m²/hr ceiling exists to catch metres read as kg/m² in the ingested rate. A
        # chosen factor is the experiment, so the product is not re-tested against it — otherwise a
        # wet enough cell cannot be run at a large scaling at all.
        adj = precipitation_adjust(stack, 1.0e4)
        @test maximum(parent(adj[:precipitation])) > 100.0
        @test metadata(adj)["precipitation_scaling"] == 1.0e4

        # Composing with a temperature perturbation re-validates, and must not resurrect the ceiling.
        @test precipitation_adjust(temperature_adjust(stack, -1.0), 1.0e4) isa DimStack
        @test temperature_adjust(precipitation_adjust(stack, 1.0e4), -1.0) isa DimStack
    end

end

println("\n✓ All climate adjustment tests passed!")
