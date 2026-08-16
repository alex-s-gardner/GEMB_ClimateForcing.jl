# Tests for the Shaw et al. (2025) per-glacier decoupling lookup table.
#
# Fully offline: the table is vendored in data/shaw2025_glacier_decoupling.csv.gz.

@testset "Glacier decoupling (Shaw et al. 2025)" begin

    @testset "table loads and is self-consistent" begin
        tbl = glacier_decoupling_table()

        # Row count of the published v1 dataset.
        @test length(tbl) == 186_792
        @test length(tbl.rgi_id) == length(tbl.k) == length(tbl.k_lower)
        @test length(tbl.lat) == length(tbl.lon) == length(tbl.ele) == length(tbl.region)

        # RGI ids are unique and well-formed.
        @test length(unique(tbl.rgi_id)) == length(tbl)
        @test all(id -> occursin(r"^RGI60-\d{2}\.\d{5}$", id), tbl.rgi_id)

        # Sorted: the ID lookup binary-searches this column, so the invariant is load-bearing.
        @test issorted(tbl.rgi_id)

        # k is clamped upstream to [0.2, 1.0]; the published global mean is ~0.83.
        @test all(0.2 .<= tbl.k .<= 1.0)
        @test isapprox(sum(tbl.k) / length(tbl.k), 0.833, atol=0.005)

        # Coordinate ranges, -180..180 longitude convention (unlike the ERA5-Land grid).
        @test all(-90 .<= tbl.lat .<= 90)
        @test all(-180 .<= tbl.lon .<= 180)
        @test minimum(tbl.lon) < 0      # confirms the -180..180 convention

        # Returns the cached object on the second call.
        @test glacier_decoupling_table() === tbl
    end

    @testset "region coverage" begin
        tbl = glacier_decoupling_table()
        regions = Set(tbl.region)

        # Regions 5 (Greenland periphery) and 19 (Antarctic) are excluded upstream, as is
        # region 20. This is a real limitation for ice-sheet forcing, so pin it.
        @test 5 ∉ regions
        @test 19 ∉ regions
        @test 1 ∈ regions && 11 ∈ regions && 18 ∈ regions
        @test minimum(tbl.region) == 1 && maximum(tbl.region) == 18
    end

    @testset "lookup by RGI id" begin
        tbl = glacier_decoupling_table()

        r = glacier_decoupling(tbl, "RGI60-01.00001")
        @test r.rgi_id == "RGI60-01.00001"
        @test r.region == 1
        @test r.distance == 0.0
        @test isapprox(r.k, 0.8436, atol=1e-4)
        @test isapprox(r.lat, 63.689, atol=1e-4)
        @test isapprox(r.lon, -146.823, atol=1e-4)
        @test isapprox(r.ele, 2385.0, atol=0.5)

        # k_lower is a lower CI *bound*, not a half-width, so it sits below k.
        @test r.k_lower < r.k

        # The convenience form (no explicit table) agrees.
        @test glacier_decoupling("RGI60-01.00001").k == r.k

        # Every row is retrievable by its own id, and round-trips to the same row.
        for i in (1, 1000, 100_000, length(tbl))
            row = glacier_decoupling(tbl, tbl.rgi_id[i])
            @test row.k == tbl.k[i]
            @test row.region == tbl.region[i]
        end
    end

    @testset "lookup by position" begin
        tbl = glacier_decoupling_table()

        # Haut Glacier d'Arolla, Switzerland — one of the study glaciers.
        r = glacier_decoupling(tbl, 45.97, 7.53)
        @test r.region == 11
        @test r.distance < 5.0
        @test 0.2 <= r.k <= 1.0

        # Looking up a glacier's own centroid returns that glacier at ~zero distance.
        i = 12_345
        r2 = glacier_decoupling(tbl, tbl.lat[i], tbl.lon[i])
        @test r2.rgi_id == tbl.rgi_id[i]
        @test r2.distance < 1e-6

        # Longitude wrapping: -146.823 and +213.177 are the same meridian.
        a = glacier_decoupling(tbl, 63.689, -146.823)
        b = glacier_decoupling(tbl, 63.689, 213.177)
        @test a.rgi_id == b.rgi_id
    end

    @testset "error handling" begin
        tbl = glacier_decoupling_table()

        # Unknown id. Includes ids sorting before and after the whole table, which the
        # binary search must handle without going out of bounds.
        @test_throws ErrorException glacier_decoupling(tbl, "RGI60-11.99999")
        @test_throws ErrorException glacier_decoupling(tbl, "not-an-rgi-id")
        @test_throws ErrorException glacier_decoupling(tbl, "AAAA")
        @test_throws ErrorException glacier_decoupling(tbl, "zzzz")
        @test_throws ErrorException glacier_decoupling(tbl, "")

        # Uncovered regions get a message naming the right region, not a bare "not found".
        @test_throws "region 05" glacier_decoupling(tbl, "RGI60-05.00001")
        @test_throws "Greenland" glacier_decoupling(tbl, "RGI60-05.00001")
        @test_throws "region 19" glacier_decoupling(tbl, "RGI60-19.00001")
        @test_throws "Antarctic" glacier_decoupling(tbl, "RGI60-19.00001")

        # ...and must not name both regions regardless of which one matched.
        err19 = try
            glacier_decoupling(tbl, "RGI60-19.00001")
        catch e
            sprint(showerror, e)
        end
        @test !occursin("Greenland", err19)

        # Out-of-range latitude.
        @test_throws ArgumentError glacier_decoupling(tbl, 91.0, 0.0)
        @test_throws ArgumentError glacier_decoupling(tbl, -91.0, 0.0)

        # A point in the open ocean has no glacier within max_distance.
        @test_throws ErrorException glacier_decoupling(tbl, 0.0, 0.0)
        # ...but succeeds when the bound is relaxed.
        @test glacier_decoupling(tbl, 0.0, 0.0; max_distance=1e5) isa NamedTuple
    end
end
