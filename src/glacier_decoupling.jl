# Per-glacier air temperature decoupling factors from Shaw et al. (2025).
#
# The vendored table (data/shaw2025_glacier_decoupling.csv.gz, built by
# data/make_shaw2025_decoupling.jl) is the *published* k per RGI v6 glacier. We read it
# rather than re-implementing the paper's five-predictor regression, whose Table S2 has
# sign errors in a4/a5 and which does not reproduce the authors' own published k even
# after correction. See docs/on_glacier_temperature_correction.md.

const _DECOUPLING_CSV = joinpath(@__DIR__, "..", "data",
                                 "shaw2025_glacier_decoupling.csv.gz")

# RGI first-order regions the upstream dataset deliberately omits, for a clearer error
# than a bare "not found" when someone asks for an ice-sheet periphery glacier.
const _UNCOVERED_REGIONS = Dict("05" => "Greenland periphery", "19" => "Antarctic and Subantarctic")

"""
    GlacierDecoupling

Lookup table of Shaw et al. (2025) air temperature decoupling factors `k` for RGI v6
glaciers. Build one with [`glacier_decoupling_table`](@ref) and query it with
[`glacier_decoupling`](@ref).

Fields are parallel vectors over glaciers: `rgi_id`, `region`, `lat`, `lon`, `ele`
(`missing` where unavailable), `k`, and `k_lower` (lower bound of the 95 % CI — an
absolute bound, not a half-width, and not clamped, so it may be negative).

`rgi_id` is sorted, which [`glacier_decoupling`](@ref) relies on to look up an identifier
by binary search rather than carrying a 186,792-entry hash table.
"""
struct GlacierDecoupling
    rgi_id::Vector{String}
    region::Vector{Int}
    lat::Vector{Float64}
    lon::Vector{Float64}
    ele::Vector{Union{Missing,Float64}}
    k::Vector{Float64}
    k_lower::Vector{Float64}
end

Base.length(t::GlacierDecoupling) = length(t.k)
function Base.show(io::IO, t::GlacierDecoupling)
    print(io, "GlacierDecoupling(", length(t), " RGI v6 glaciers, regions ",
          minimum(t.region), "–", maximum(t.region), ")")
end

# Cached so repeated calls in a session don't re-parse ~10 MB of CSV.
const _DECOUPLING_CACHE = Ref{Union{Nothing,GlacierDecoupling}}(nothing)

"""
    glacier_decoupling_table() -> GlacierDecoupling

Load (and cache) the Shaw et al. (2025) per-glacier decoupling table.

The table covers 186,792 RGI v6 glaciers. **RGI regions 05 (Greenland periphery) and 19
(Antarctic) are absent** — the upstream dataset excludes ice sheets, ice caps and the
large Antarctic/Greenland periphery glaciers — so lookups there will fail. Use the SM10
flowpath-length scheme for those instead.

# Example
```julia
tbl = glacier_decoupling_table()
glacier_decoupling(tbl, "RGI60-11.00001")
```
"""
function glacier_decoupling_table()
    cached = _DECOUPLING_CACHE[]
    cached === nothing || return cached

    isfile(_DECOUPLING_CSV) ||
        error("decoupling table not found at $_DECOUPLING_CSV; regenerate it with " *
              "data/make_shaw2025_decoupling.jl")

    rgi_id = String[]
    region = Int[]
    lat = Float64[]
    lon = Float64[]
    ele = Union{Missing,Float64}[]
    k = Float64[]
    k_lower = Float64[]
    # Row count of the shipped table; a hint only, so a regenerated file of any size works.
    for v in (rgi_id, region, lat, lon, ele, k, k_lower)
        sizehint!(v, 187_000)
    end

    open(GzipDecompressorStream, _DECOUPLING_CSV) do io
        for line in eachline(io)
            (isempty(line) || startswith(line, '#')) && continue
            startswith(line, "rgi_id,") && continue      # header
            f = split(line, ',')
            length(f) == 7 || error("malformed row in $_DECOUPLING_CSV: $line")
            push!(rgi_id, f[1])
            push!(region, parse(Int, f[2]))
            push!(lat, parse(Float64, f[3]))
            push!(lon, parse(Float64, f[4]))
            push!(ele, isempty(f[5]) ? missing : parse(Float64, f[5]))
            push!(k, parse(Float64, f[6]))
            push!(k_lower, parse(Float64, f[7]))
        end
    end

    # ID lookup binary-searches this column, so a regenerated file that is no longer
    # sorted must fail loudly here rather than silently missing glaciers.
    issorted(rgi_id) ||
        error("$_DECOUPLING_CSV is not sorted by rgi_id; regenerate it with " *
              "data/make_shaw2025_decoupling.jl, which preserves the upstream order")

    tbl = GlacierDecoupling(rgi_id, region, lat, lon, ele, k, k_lower)
    _DECOUPLING_CACHE[] = tbl
    return tbl
end

"""
    glacier_decoupling(table, rgi_id::AbstractString) -> NamedTuple
    glacier_decoupling(table, lat, lon; max_distance=10.0) -> NamedTuple

Look up the Shaw et al. (2025) decoupling factor for a glacier, either by RGI v6
identifier (exact) or by position (nearest glacier centroid).

Returns `(; rgi_id, region, lat, lon, ele, k, k_lower, distance)`, where `k` is the
unitless decoupling factor — on-glacier temperature change per unit ambient change,
clamped upstream to `[0.2, 1.0]` — and `distance` is the great-circle distance in km to
the matched centroid (`0.0` for an ID lookup).

The positional form matches against glacier *mean centroids*, so on a large glacier the
nearest centroid may be some distance from the point of interest; `max_distance` (km)
bounds how far a match may be before an error is raised. Prefer the RGI-ID form when the
identifier is known.

!!! warning "Ambient temperature, ablation season"
    `k` multiplies an *ambient* (off-glacier) temperature that has already been adjusted
    to the glacier's elevation — run [`climate_adjust_for_elevation`](@ref) first. The
    underlying regression is calibrated on ablation-season data only.
"""
function glacier_decoupling(t::GlacierDecoupling, rgi_id::AbstractString)
    id = String(rgi_id)
    i = searchsortedfirst(t.rgi_id, id)          # rgi_id is sorted; checked at load
    if i > length(t) || t.rgi_id[i] != id
        m = match(r"^RGI60-(\d{2})\.", id)
        if m !== nothing && haskey(_UNCOVERED_REGIONS, m.captures[1])
            region = m.captures[1]
            error("no decoupling factor for $id: RGI region $region " *
                  "($(_UNCOVERED_REGIONS[region])) is not covered by the Shaw et al. " *
                  "(2025) dataset, which excludes ice sheets, ice caps and large " *
                  "periphery glaciers")
        end
        error("RGI id $id not found in the Shaw et al. (2025) decoupling table")
    end
    return _decoupling_row(t, i, 0.0)
end

function glacier_decoupling(t::GlacierDecoupling, lat::Real, lon::Real;
                            max_distance::Real=10.0)
    -90 <= lat <= 90 || throw(ArgumentError("latitude must be in [-90, 90], got $lat"))
    lon = _wrap_longitude(lon)
    best, best_d = _nearest_centroid(t, lat, lon, max_distance)
    best_d <= max_distance ||
        error("nearest glacier centroid in the Shaw et al. (2025) table is " *
              "$(round(best_d, digits=1)) km from ($lat, $lon), beyond max_distance=" *
              "$max_distance km. Note RGI regions 05 and 19 are not covered.")
    return _decoupling_row(t, best, best_d)
end

glacier_decoupling(rgi_id::AbstractString) =
    glacier_decoupling(glacier_decoupling_table(), rgi_id)
glacier_decoupling(lat::Real, lon::Real; kwargs...) =
    glacier_decoupling(glacier_decoupling_table(), lat, lon; kwargs...)

function _decoupling_row(t::GlacierDecoupling, i::Int, distance::Float64)
    return (rgi_id=t.rgi_id[i], region=t.region[i], lat=t.lat[i], lon=t.lon[i],
            ele=t.ele[i], k=t.k[i], k_lower=t.k_lower[i], distance=distance)
end

# Index of the nearest glacier centroid to (lat, lon), and its distance in km.
#
# A latitude band cheaply rejects most of the 186,792 rows before the haversine: the
# great-circle distance is at least 110.57 km per degree of latitude, so a glacier more
# than max_distance/110.5° away in latitude cannot be within max_distance km. That makes
# the band exact for any match the caller would accept, ~150x faster than the full scan.
#
# It is *not* exact for a match the caller will reject — the nearest glacier inside the
# band can be much farther than the nearest glacier overall. So when the band yields
# nothing acceptable, rescan everything, which costs a full pass only on the error path
# and keeps the reported distance truthful.
function _nearest_centroid(t::GlacierDecoupling, lat::Real, lon::Real, max_distance::Real)
    scan(rows) = begin
        best, best_d = 0, Inf
        @inbounds for i in rows
            d = _haversine_km(lat, lon, t.lat[i], t.lon[i])
            if d < best_d
                best_d, best = d, i
            end
        end
        best, best_d
    end

    Δlat = max_distance / 110.5
    band = Iterators.filter(i -> abs(t.lat[i] - lat) <= Δlat, eachindex(t.k))
    best, best_d = scan(band)
    best_d <= max_distance && return best, best_d

    return scan(eachindex(t.k))
end

function _haversine_km(lat1::Real, lon1::Real, lat2::Real, lon2::Real)
    R = 6371.0088                       # mean Earth radius, km
    φ1, φ2 = deg2rad(lat1), deg2rad(lat2)
    Δφ = φ2 - φ1
    Δλ = deg2rad(lon2 - lon1)
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    return 2R * asin(min(1.0, sqrt(a)))
end
