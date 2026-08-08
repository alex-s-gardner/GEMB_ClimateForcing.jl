"""
verify_chunk_alignment.jl

Empirically confirms that `era5_land_chunk_map` assigns chunk IDs that agree with the
actual Zarr chunk files fetched from the live ARCO store.

Strategy
--------
1. Open the store with a DirectoryStore cache so fetched chunks land as files on disk.
2. Read `metadata.chunks` and note the claimed (lon_chunk_size, lat_chunk_size) tuple.
3. Fetch single elements at several deliberately chosen (lon, lat) positions that span
   known chunk boundaries:
     - two points inside the same chunk
     - two points in adjacent lon chunks
     - two points in adjacent lat chunks
4. After each fetch, inspect which new chunk files appeared in the cache directory.
   Zarr v2 chunk filenames encode the chunk indices: `<array>/<lon_cid>.<lat_cid>.<time_cid>`
   (following the same empirical on-disk dimension order as the data reads).
5. Compare the filename-derived chunk indices against what `(i-1) ÷ chunk_size` predicts.
   A mismatch means either `chunks[1]`/`chunks[2]` are swapped, or the formula is wrong.

Run with:
    julia --project=. test/verify_chunk_alignment.jl

Requires CDS_API_KEY to be set.
"""

using Zarr
using Dates

# ── load the package so AuthenticatedHTTPStore etc. are available ─────────────
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GEMB_ClimateForcing: era5_land_url
# We need the internal store type — reach into the module
using GEMB_ClimateForcing

token = get(ENV, "CDS_API_KEY", nothing)
if isnothing(token) || isempty(strip(token))
    error("CDS_API_KEY environment variable is not set. " *
          "Get a free key at https://cds.climate.copernicus.eu/")
end

const VARIABLE_GROUP = "sfc-2m-temperature"
const VAR_NAME       = "t2m"

# Use a temp dir so we start with an empty cache each run
cache_dir = mktempdir(; prefix="zarr_chunk_verify_")
println("Cache directory: ", cache_dir)

# ── open store with caching ───────────────────────────────────────────────────
for chunk_strategy in (:geo, :time)
    println("\n", "="^70)
    println("Verifying chunk_strategy = :$(chunk_strategy)")
    println("="^70)

    url   = era5_land_url(VARIABLE_GROUP, chunk_strategy)
    # Access internal constructor via the module
    store = GEMB_ClimateForcing.AuthenticatedHTTPStore(url; token=token)
    sub_cache = joinpath(cache_dir, string(chunk_strategy))
    cstore = Zarr.CachingStore(store, Zarr.DirectoryStore(sub_cache))

    zg = Zarr.zopen(cstore, consolidated=true, fill_as_missing=false)

    lat_values = zg["latitude"][:]
    lon_values = zg["longitude"][:]
    time_values_raw = zg["time"][:]

    n_lon = length(lon_values)
    n_lat = length(lat_values)

    zarr_arr = zg[VAR_NAME]
    chunks   = zarr_arr.metadata.chunks  # claimed NTuple{3,Int}
    println("metadata.chunks = ", chunks)
    println("  claimed → chunks[1]=$(chunks[1]) (lon?), chunks[2]=$(chunks[2]) (lat?), chunks[3]=$(chunks[3]) (time?)")
    println("Grid: n_lon=$(n_lon), n_lat=$(n_lat)")

    # ── choose test indices ───────────────────────────────────────────────────
    # Pick a lon and lat that are NOT at chunk boundaries so we can jump over one.
    # Use offset 10 into the first chunk to stay away from edges.
    lon_i_A = 10                        # in chunk 0
    lon_i_B = 10                        # same lon chunk
    lon_i_C = chunks[1] + 10            # next lon chunk (chunk 1)

    lat_j_A = 10                        # in chunk 0
    lat_j_B = chunks[2] + 10            # next lat chunk (chunk 1)
    lat_j_C = 10                        # same lat chunk as A

    # Use time index 1 (first time step) — we just need one element
    t_idx = 1

    # Predicted chunk IDs from our formula
    n_lat_chunks = cld(n_lat, chunks[2])

    function predicted_chunk_id(i, j)
        lon_cid = (i - 1) ÷ chunks[1]
        lat_cid = (j - 1) ÷ chunks[2]
        return (lon_cid, lat_cid, lon_cid * n_lat_chunks + lat_cid)
    end

    println("\nTest points (1-based grid indices):")
    println("  A: lon_i=$(lon_i_A), lat_j=$(lat_j_A) → predicted chunk $(predicted_chunk_id(lon_i_A, lat_j_A))")
    println("  B: lon_i=$(lon_i_B), lat_j=$(lat_j_B) → predicted chunk $(predicted_chunk_id(lon_i_B, lat_j_B))  [same lon chunk, next lat chunk]")
    println("  C: lon_i=$(lon_i_C), lat_j=$(lat_j_C) → predicted chunk $(predicted_chunk_id(lon_i_C, lat_j_C))  [next lon chunk, same lat chunk]")

    # ── fetch elements and record which chunk files appear ────────────────────
    # Clear any pre-cached t2m chunks (consolidated metadata may have already been fetched)
    var_cache = joinpath(sub_cache, VAR_NAME)

    function cached_chunks_now()
        isdir(var_cache) || return Set{String}()
        return Set(readdir(var_cache))
    end

    println("\nFetching element A (lon=$(lon_i_A), lat=$(lat_j_A), t=1)...")
    before = cached_chunks_now()
    _ = zarr_arr[lon_i_A, lat_j_A, t_idx]
    after_A = cached_chunks_now()
    new_A = setdiff(after_A, before)
    println("  New chunk files: ", isempty(new_A) ? "(none — may have been pre-cached)" : join(sort(collect(new_A)), ", "))

    println("\nFetching element B (lon=$(lon_i_B), lat=$(lat_j_B), t=1) — different lat chunk...")
    before = cached_chunks_now()
    _ = zarr_arr[lon_i_B, lat_j_B, t_idx]
    after_B = cached_chunks_now()
    new_B = setdiff(after_B, before)
    println("  New chunk files: ", isempty(new_B) ? "(none — already cached)" : join(sort(collect(new_B)), ", "))

    println("\nFetching element C (lon=$(lon_i_C), lat=$(lat_j_C), t=1) — different lon chunk...")
    before = cached_chunks_now()
    _ = zarr_arr[lon_i_C, lat_j_C, t_idx]
    after_C = cached_chunks_now()
    new_C = setdiff(after_C, before)
    println("  New chunk files: ", isempty(new_C) ? "(none — already cached)" : join(sort(collect(new_C)), ", "))

    println("\nAll cached chunk files for $(VAR_NAME): ", sort(collect(cached_chunks_now())))

    # ── parse chunk filenames and verify ─────────────────────────────────────
    # Zarr v2 filenames: "<d0>.<d1>.<d2>"
    # We check whether the dimension that changes between A→B matches the lat axis
    # and the dimension that changes between A→C matches the lon axis.
    println("\n--- Alignment verification ---")

    all_chunks = sort(collect(cached_chunks_now()))
    if isempty(all_chunks)
        println("  WARNING: no chunk files found in cache — cannot verify from filenames.")
        println("  (The CachingStore may not have written to disk, or the store is read-only.)")
    else
        # Parse all chunk filenames into tuples of ints
        parsed = [Tuple(parse(Int, x) for x in split(f, ".")) for f in all_chunks if occursin(".", f)]
        println("  Cached chunk index tuples: ", parsed)

        # Group by what changed
        # A's chunk: the one fetched at (lon_i_A, lat_j_A)
        a_lon_cid = (lon_i_A - 1) ÷ chunks[1]
        a_lat_cid = (lat_j_A - 1) ÷ chunks[2]
        b_lon_cid = (lon_i_B - 1) ÷ chunks[1]  # same lon chunk as A
        b_lat_cid = (lat_j_B - 1) ÷ chunks[2]  # different lat chunk
        c_lon_cid = (lon_i_C - 1) ÷ chunks[1]  # different lon chunk
        c_lat_cid = (lat_j_C - 1) ÷ chunks[2]  # same lat chunk as A

        println("  Predicted chunk tuples under assumption chunks=(lon,lat,time):")
        println("    A → ($(a_lon_cid), $(a_lat_cid), 0)")
        println("    B → ($(b_lon_cid), $(b_lat_cid), 0)  [lat increments]")
        println("    C → ($(c_lon_cid), $(c_lat_cid), 0)  [lon increments]")

        expected_A = (a_lon_cid, a_lat_cid, 0)
        expected_B = (b_lon_cid, b_lat_cid, 0)
        expected_C = (c_lon_cid, c_lat_cid, 0)

        ok_A = expected_A in parsed
        ok_B = expected_B in parsed
        ok_C = expected_C in parsed

        println("\n  Check A in cache: ", ok_A ? "PASS ✓" : "FAIL ✗  (expected $(expected_A))")
        println("  Check B in cache: ", ok_B ? "PASS ✓" : "FAIL ✗  (expected $(expected_B))")
        println("  Check C in cache: ", ok_C ? "PASS ✓" : "FAIL ✗  (expected $(expected_C))")

        # Also check the swapped interpretation
        if !ok_A || !ok_B || !ok_C
            println("\n  Trying swapped interpretation: chunks=(lat,lon,time)...")
            a_lat_swap = (lon_i_A - 1) ÷ chunks[2]
            a_lon_swap = (lat_j_A - 1) ÷ chunks[1]
            b_lat_swap = (lon_i_B - 1) ÷ chunks[2]
            b_lon_swap = (lat_j_B - 1) ÷ chunks[1]
            c_lat_swap = (lon_i_C - 1) ÷ chunks[2]
            c_lon_swap = (lat_j_C - 1) ÷ chunks[1]

            exp_A_swap = (a_lat_swap, a_lon_swap, 0)
            exp_B_swap = (b_lat_swap, b_lon_swap, 0)
            exp_C_swap = (c_lat_swap, c_lon_swap, 0)

            println("    A → ($(exp_A_swap[1]), $(exp_A_swap[2]), 0) : ", exp_A_swap in parsed ? "PASS ✓" : "FAIL ✗")
            println("    B → ($(exp_B_swap[1]), $(exp_B_swap[2]), 0) : ", exp_B_swap in parsed ? "PASS ✓" : "FAIL ✗")
            println("    C → ($(exp_C_swap[1]), $(exp_C_swap[2]), 0) : ", exp_C_swap in parsed ? "PASS ✓" : "FAIL ✗")
        end

        if ok_A && ok_B && ok_C
            println("\n  RESULT: chunks[1]=lon, chunks[2]=lat confirmed. Chunk ID formula is CORRECT.")
        else
            println("\n  RESULT: Mismatch detected — chunk ID formula may be WRONG.")
            println("  Review the raw chunk filenames above and the two interpretation tables.")
        end
    end
end

println("\nDone. Cache written to: ", cache_dir)
