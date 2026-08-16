#!/usr/bin/env julia

"""
Analytic (offline) tests for the ambient → on-glacier correction
(`climate_adjust_for_glacier`). Uses the vendored Shaw et al. (2025) lookup table for
the ID/position forms, so no network access is required.
"""

using Test
using GEMB_ClimateForcing
using Dates
using DimensionalData
using Statistics

using GEMB_ClimateForcing: saturation_vapor_pressure,
    konzelmann_clear_sky_emissivity, _SIGMA_SB

# Melt-season stack straddling 0 °C, so both branches of the gated correction are
# exercised. Positioned on Haut Glacier d'Arolla so the metadata lookup resolves.
function make_glacier_stack(; npts::Int=48)
    times = DateTime(2020, 7, 15, 0) .+ Hour.(0:(npts - 1))
    time_dim = Ti(times)
    # ±6 K about freezing over a diurnal cycle.
    T = [273.15 + 6.0 * sinpi(2 * (i - 1) / 24) for i in 1:npts]
    e = 0.7 .* saturation_vapor_pressure.(T)
    return DimStack((
        temperature_air = DimArray(T, (time_dim,); metadata=Dict("units" => "K")),
        pressure_air = DimArray(fill(72000.0, npts), (time_dim,); metadata=Dict("units" => "Pa")),
        vapor_pressure = DimArray(e, (time_dim,); metadata=Dict("units" => "Pa")),
        wind_speed = DimArray(fill(3.0, npts), (time_dim,); metadata=Dict("units" => "m/s")),
        precipitation = DimArray(fill(0.1, npts), (time_dim,); metadata=Dict("units" => "kg/m²/hr")),
        shortwave_downward = DimArray(fill(250.0, npts), (time_dim,); metadata=Dict("units" => "W/m²")),
        longwave_downward = DimArray(fill(260.0, npts), (time_dim,); metadata=Dict("units" => "W/m²")),
    ); metadata=Dict("latitude" => 45.97, "longitude" => 7.53, "elevation" => 2900.0,
                     "dataset" => "test"))
end

# Variables the glacier correction must never touch.
const _GLACIER_UNCHANGED = (:pressure_air, :vapor_pressure, :wind_speed, :precipitation,
    :shortwave_downward)

@testset "Glacier Adjustment (ambient → on-glacier)" begin

    @testset "k = 1 is the identity" begin
        stack = make_glacier_stack()
        adj = climate_adjust_for_glacier(stack, 1.0)

        # Temperature is bit-exact, not merely approximate: the increment form
        # multiplies by (k - 1) = 0.
        @test parent(adj[:temperature_air]) == parent(stack[:temperature_air])
        for var in _GLACIER_UNCHANGED
            @test parent(adj[var]) == parent(stack[var])
        end
        # Longwave round-trips through the emissivity diagnosis to float tolerance.
        @test all(isapprox.(parent(adj[:longwave_downward]),
                            parent(stack[:longwave_downward]); rtol=1e-12))
        # Identity holds in the ungated branch too.
        adj_ub = climate_adjust_for_glacier(stack, 1.0; apply_below_freezing=true)
        @test parent(adj_ub[:temperature_air]) == parent(stack[:temperature_air])
    end

    @testset "cooling is gated at the melting point" begin
        stack = make_glacier_stack()
        k = 0.8
        adj = climate_adjust_for_glacier(stack, k)

        T = parent(stack[:temperature_air])
        T′ = parent(adj[:temperature_air])
        warm = T .> 273.15

        @test any(warm) && any(.!warm)          # the fixture exercises both branches
        # Below/at freezing: untouched by default.
        @test T′[.!warm] == T[.!warm]
        # Above freezing: cooled, by exactly (1-k)·(T - T_ref).
        @test all(T′[warm] .< T[warm])
        @test all(isapprox.(T[warm] .- T′[warm], (1 - k) .* (T[warm] .- 273.15)))
        # Cooling never overshoots past the melting point (k > 0).
        @test all(T′[warm] .> 273.15)
    end

    @testset "below-freezing branch" begin
        stack = make_glacier_stack()
        T = parent(stack[:temperature_air])
        cold = T .< 273.15

        adj = climate_adjust_for_glacier(stack, 0.8; apply_below_freezing=true)
        T′ = parent(adj[:temperature_air])

        # The unphysical consequence this kwarg opts into: k·T warms sub-freezing air.
        # Pinned deliberately, since it is why the default is false.
        @test all(T′[cold] .> T[cold])
        @test all(isapprox.(T′, 273.15 .+ 0.8 .* (T .- 273.15)))
    end

    @testset "reference_temperature shifts the pivot" begin
        stack = make_glacier_stack()
        T = parent(stack[:temperature_air])

        # Pivoting about the series maximum leaves everything untouched under gating,
        # since no step exceeds the threshold.
        adj = climate_adjust_for_glacier(stack, 0.5; reference_temperature=maximum(T))
        @test parent(adj[:temperature_air]) == T

        # A lower pivot cools every step.
        adj = climate_adjust_for_glacier(stack, 0.5; reference_temperature=minimum(T) - 1)
        @test all(parent(adj[:temperature_air]) .< T)
    end

    @testset "only temperature and longwave change" begin
        stack = make_glacier_stack()
        adj = climate_adjust_for_glacier(stack, 0.7)

        for var in _GLACIER_UNCHANGED
            @test parent(adj[var]) == parent(stack[var])
        end

        # Vapor pressure is left alone *on purpose* — SM10 Eq. 4 pivots e about
        # 6.11 hPa, opposite in sign to constant-RH propagation. Contrast with
        # temperature_adjust, which does propagate at constant RH.
        @test parent(adj[:vapor_pressure]) == parent(stack[:vapor_pressure])
        ta = temperature_adjust(stack, -1.0)
        @test parent(ta[:vapor_pressure]) != parent(stack[:vapor_pressure])

        # Longwave falls wherever the air was cooled, and is unchanged elsewhere.
        T = parent(stack[:temperature_air])
        warm = T .> 273.15
        LW = parent(stack[:longwave_downward])
        LW′ = parent(adj[:longwave_downward])
        @test all(LW′[warm] .< LW[warm])
        @test all(isapprox.(LW′[.!warm], LW[.!warm]; rtol=1e-12))

        # And it matches an explicit Konzelmann recomputation at unchanged e.
        e = parent(stack[:vapor_pressure])
        T′ = parent(adj[:temperature_air])
        Δε = @. LW / (_SIGMA_SB * T^4) - konzelmann_clear_sky_emissivity(e, T)
        expected = @. (konzelmann_clear_sky_emissivity(e, T′) + Δε) * _SIGMA_SB * T′^4
        @test all(isapprox.(LW′, expected))
    end

    @testset "metadata" begin
        stack = make_glacier_stack()
        adj = climate_adjust_for_glacier(stack, 0.75)
        m = metadata(adj)

        @test m["decoupling_factor"] == 0.75
        @test m["glacier_decoupling_factor"] == 0.75
        @test m["glacier_decoupling_reference_temperature"] == 273.15
        @test m["glacier_decoupling_below_freezing"] == false
        @test m["temperature_air_mean"] ≈ mean(parent(adj[:temperature_air]))
        # Source metadata is carried through.
        @test m["elevation"] == 2900.0
        @test m["dataset"] == "test"
        # No lookup happened, so no provenance keys.
        @test !haskey(m, "glacier_decoupling_rgi_id")

        # Repeated calls compose multiplicatively.
        twice = climate_adjust_for_glacier(climate_adjust_for_glacier(stack, 0.9), 0.9)
        @test metadata(twice)["glacier_decoupling_factor"] ≈ 0.81
        @test metadata(twice)["decoupling_factor"] == 0.9      # most recent call only
    end

    @testset "lookup forms" begin
        stack = make_glacier_stack()
        expected = glacier_decoupling("RGI60-11.02810")

        by_id = climate_adjust_for_glacier(stack; rgi_id="RGI60-11.02810")
        # The fixture sits on Arolla, so the implicit centroid lookup finds the same
        # glacier as the explicit ID.
        by_pos = climate_adjust_for_glacier(stack)

        for adj in (by_id, by_pos)
            m = metadata(adj)
            @test m["decoupling_factor"] == expected.k
            @test m["glacier_decoupling_rgi_id"] == "RGI60-11.02810"
            @test m["glacier_decoupling_k_lower"] == expected.k_lower
            @test occursin("Shaw", m["glacier_decoupling_source"])
            @test haskey(m, "glacier_decoupling_match_distance")
        end
        @test metadata(by_id)["glacier_decoupling_match_distance"] == 0.0
        @test metadata(by_pos)["glacier_decoupling_match_distance"] > 0.0

        # Same k reached by either route, so the fields agree.
        @test parent(by_id[:temperature_air]) == parent(by_pos[:temperature_air])
        # And it agrees with passing k explicitly.
        explicit = climate_adjust_for_glacier(stack, expected.k)
        @test parent(explicit[:temperature_air]) == parent(by_id[:temperature_air])
    end

    @testset "error handling" begin
        stack = make_glacier_stack()

        # k outside (0, 1].
        @test_throws ArgumentError climate_adjust_for_glacier(stack, 1.2)
        @test_throws ArgumentError climate_adjust_for_glacier(stack, 0.0)
        @test_throws ArgumentError climate_adjust_for_glacier(stack, -0.5)
        # k = 1 is allowed (the identity).
        @test climate_adjust_for_glacier(stack, 1.0) isa DimStack

        # Uncovered RGI regions surface the table's explanatory error.
        @test_throws ErrorException climate_adjust_for_glacier(stack; rgi_id="RGI60-19.00001")
        @test_throws ErrorException climate_adjust_for_glacier(stack; rgi_id="RGI60-05.00001")

        # Implicit lookup needs coordinates in the metadata.
        bare = rebuild(stack; metadata=Dict("dataset" => "test"))
        @test_throws ArgumentError climate_adjust_for_glacier(bare)

        # A location far from any glacier is rejected by the distance bound.
        ocean = rebuild(stack; metadata=Dict("latitude" => 0.0, "longitude" => 0.0,
                                            "dataset" => "test"))
        @test_throws ErrorException climate_adjust_for_glacier(ocean)
    end

    @testset "composes with the elevation correction" begin
        stack = make_glacier_stack()

        # Documented order: elevation first, glacier second. Both metadata records
        # survive the composition.
        at_z = climate_adjust_for_elevation(stack, -200.0)
        on_glacier = climate_adjust_for_glacier(at_z, 0.8)
        m = metadata(on_glacier)
        @test m["elevation_offset"] == -200.0
        @test m["glacier_decoupling_factor"] == 0.8
        @test m["elevation"] == 2700.0

        # The two orders differ (the glacier correction is nonlinear in T through the
        # gate), which is why the order is specified rather than left to the caller.
        other_order = climate_adjust_for_elevation(climate_adjust_for_glacier(stack, 0.8), -200.0)
        @test parent(on_glacier[:temperature_air]) != parent(other_order[:temperature_air])
    end
end

println("\n✓ All glacier adjustment tests passed!")
