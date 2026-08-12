using Test
using GEMB_ClimateForcing
using Rasters

# Independent reference implementation of the List (1968) geopotential → geometric-height
# conversion, with constants re-typed fresh so a coefficient/sign error in the package is
# caught rather than mirrored (same convention as test_elevation_physics.jl).
const REF_G0 = 9.80665            # standard gravity, m/s²
ref_gravity(ϕdeg) = (ϕ = deg2rad(ϕdeg); 9.80616 * (1 - 0.0026373 * cos(2ϕ) + 0.0000059 * cos(2ϕ)^2))
ref_radius(ϕdeg)  = (ϕ = deg2rad(ϕdeg); 6378137.0 / (1.006803 - 0.006706 * sin(ϕ)^2))
function ref_geometric(Φ, ϕdeg)
    H = Φ / REF_G0
    g = ref_gravity(ϕdeg)
    R = ref_radius(ϕdeg)
    return R * H / ((g / REF_G0) * R - H)
end

@testset "Geoid / geopotential2height" begin

    @testset "Surface gravity (List 1968) vs published values" begin
        # Equator ≈ 9.780, pole ≈ 9.832 m/s² (standard normal-gravity endpoints).
        @test ref_gravity(0.0)  ≈ 9.7803 atol=3e-3
        @test ref_gravity(90.0) ≈ 9.8322 atol=3e-3
        # Gravity increases monotonically from equator to pole.
        @test ref_gravity(0.0) < ref_gravity(45.0) < ref_gravity(90.0)
    end

    @testset "Geometric height (no geoid) matches independent reference" begin
        # Use the package's internal geometric-height step (steps 1-4, before adding N).
        gh = GEMB_ClimateForcing._geopotential_to_geometric
        # Summit-ish geopotential (~3216 m of geopotential height).
        Φ = 3216.0 * REF_G0
        for ϕ in (-80.0, -30.0, 0.0, 45.0, 72.58, 89.0)
            @test gh(Φ, ϕ) ≈ ref_geometric(Φ, ϕ) atol=1e-6
        end
        # Zero geopotential → zero height at every latitude.
        @test gh(0.0, 45.0) == 0.0
        # Geometric height slightly exceeds naive Φ/g₀ (gravity < g₀ over most of the globe,
        # and radius correction lifts it); check it is within a sensible band.
        @test gh(Φ, 45.0) ≈ 3216.0 rtol=0.01
    end

    @testset "geopotential2height with supplied geoid raster (no network)" begin
        # Build a tiny synthetic geoid grid with a known undulation so we can test the
        # end-to-end arithmetic h = Z + N without streaming from the CDN.
        lons = X(-180.0:10.0:180.0)
        lats = Y(-90.0:10.0:90.0)
        Ngrid = fill(30.0, length(lons), length(lats))   # constant N = 30 m
        geoid = Raster(Ngrid, (lons, lats))

        Φ = 3216.0 * REF_G0
        ϕ, λ = 72.58, -38.46
        N = geoid_undulation(ϕ, λ; geoid_raster=geoid)
        @test N ≈ 30.0 atol=1e-9

        h = geopotential2height(Φ, ϕ, λ; geoid_raster=geoid)
        @test h ≈ ref_geometric(Φ, ϕ) + 30.0 atol=1e-6

        # Longitude in 0–360°E convention is wrapped to the grid's −180…180.
        @test geoid_undulation(ϕ, 321.54; geoid_raster=geoid) ≈
              geoid_undulation(ϕ, -38.46; geoid_raster=geoid) atol=1e-9

        # Array broadcasting: scalars and vectors both work, shapes preserved.
        Φs = [Φ, 2Φ, 0.0]
        ϕs = [10.0, -45.0, 80.0]
        λs = [0.0, 100.0, 250.0]
        hs = geopotential2height(Φs, ϕs, λs; geoid_raster=geoid)
        @test length(hs) == 3
        @test hs[3] ≈ 30.0 atol=1e-6   # zero geopotential → just N
        for i in 1:3
            @test hs[i] ≈ ref_geometric(Φs[i], ϕs[i]) + 30.0 atol=1e-6
        end
    end

    @testset "height_reference = :orthometric (no network)" begin
        Φ = 3216.0 * REF_G0
        ϕ, λ = 72.58, -38.46
        # :orthometric returns the geometric height above the geoid directly (no N, no geoid
        # streamed), so it equals the independent reference geometric height.
        @test geopotential2height(Φ, ϕ, λ; height_reference=:orthometric) ≈ ref_geometric(Φ, ϕ) atol=1e-6
        # Ellipsoidal (with a supplied geoid) exceeds orthometric by exactly N.
        lons = X(-180.0:10.0:180.0); lats = Y(-90.0:10.0:90.0)
        geoid = Raster(fill(30.0, length(lons), length(lats)), (lons, lats))
        h_ell = geopotential2height(Φ, ϕ, λ; height_reference=:wgs84, geoid_raster=geoid)
        h_ort = geopotential2height(Φ, ϕ, λ; height_reference=:orthometric)
        @test h_ell - h_ort ≈ 30.0 atol=1e-6
        # Array broadcasting on the orthometric path.
        hs = geopotential2height([Φ, 0.0], [10.0, 80.0], [0.0, 250.0]; height_reference=:orthometric)
        @test hs[2] == 0.0
    end

    @testset "Argument validation (offline)" begin
        # geopotential2height with an unsupported ellipsoid errors before any network use.
        @test_throws ArgumentError geopotential2height(1000.0, 45.0, 0.0; ellipsoid=:grs80)
        # Unsupported height_reference.
        @test_throws ArgumentError geopotential2height(1000.0, 45.0, 0.0; height_reference=:navd88)
        # Unknown geoid model.
        @test_throws ArgumentError GEMB_ClimateForcing._load_geoid(:not_a_geoid)
        # Registry advertises the two streamable grids.
        @test haskey(GEOID_MODELS, :egm96)
        @test haskey(GEOID_MODELS, :egm2008)
    end

    # ------------------------------------------------------------------
    # Integration: real geoid grid streamed via /vsicurl/ (opt-in, no token).
    # ------------------------------------------------------------------
    run_integration = get(ENV, "GEMB_TEST_GEOID", "") in ("1", "true", "TRUE")

    if run_integration
        println("\n" * "="^70)
        println("Running geoid streaming tests (network; EGM96 ~2.6 MB via /vsicurl/)")
        println("="^70)

        @testset "Stream EGM96 and sample undulation" begin
            g = GEMB_ClimateForcing._load_geoid(:egm96; verbose=false)
            @test g isa Rasters.AbstractRaster
            @test hasdim(g, X)
            @test hasdim(g, Y)

            # Global undulation lies within the known EGM96 range (~ -107 … +86 m).
            N_summit = geoid_undulation(72.58, -38.46; geoid=:egm96)
            @test -110.0 < N_summit < 90.0

            # Ellipsoidal height = geometric + N, consistent with the sampled undulation.
            Φ = 3216.0 * REF_G0
            h = geopotential2height(Φ, 72.58, -38.46; geoid=:egm96)
            @test h ≈ ref_geometric(Φ, 72.58) + N_summit atol=1.0
            # ... and differs from the naive Φ/g₀ elevation by roughly N (tens of metres).
            @test abs(h - Φ / REF_G0) > 5.0
        end
    else
        println("\n⚠️  Skipping geoid streaming tests (set GEMB_TEST_GEOID=1 to enable)")
    end
end
