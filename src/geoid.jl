"""
Geoid undulation and geopotential → **ellipsoidal-height** conversion.

ERA5-Land geopotential `:z` (m² s⁻²) is conventionally turned into elevation with the
crude constant-gravity division `z / 9.80665`, which yields an approximate *geopotential
height* referenced (roughly) to mean sea level — not a geometric height, and not height
above the WGS84 ellipsoid. [`geopotential2height`](@ref) does the full conversion:

1. geopotential → geometric (orthometric) height above the geoid, using latitude-dependent
   gravity and effective Earth radius (List 1968, as in NCL's `gp2gmh`), then
2. adds the geoid undulation `N` so the result is height above the **WGS84 ellipsoid**.

The geoid model is **streamed** from the open PROJ CDN (`cdn.proj.org`, published by
OSGeo/PROJ-data from US NGA) as a Cloud-Optimized GeoTIFF read over HTTP byte ranges via
GDAL `/vsicurl/` — nothing is bundled with the package. The same streaming machinery and
the macOS CA-bundle fix (`_configure_gdal_http`) are shared with the Copernicus DEM loader
(`datasets/copernicus_dem.jl`).

Available geoid grids ([`GEOID_MODELS`](@ref)):
- `:egm96`   — EGM96, 15′ grid  (`us_nga_egm96_15.tif`, ~2.6 MB) — default.
- `:egm2008` — EGM2008, 2.5′ grid (`us_nga_egm08_25.tif`, ~77 MB) — higher resolution.

Each grid stores `N`, the height of the geoid above the WGS84 ellipsoid, i.e. the value
added to an orthometric height to obtain an ellipsoidal height. Both use −180…180°E
longitude and are referenced to WGS 84 (EPSG:4979).
"""

using Rasters
using ArchGDAL          # loads Rasters' GDAL backend (RastersArchGDALExt)
using DimensionalData

# Base URL of the OSGeo PROJ datum-grid CDN. Anonymous HTTPS, no authentication.
const _PROJ_CDN_BASE = "https://cdn.proj.org"

# Streamable geoid grids, keyed by model symbol → (GeoTIFF filename, resolution, description).
# N in these grids is geoid height above the WGS84 ellipsoid (WGS 84 → EGMxx height).
const GEOID_MODELS = Dict{Symbol,NamedTuple{(:file, :resolution, :description),Tuple{String,String,String}}}(
    :egm96   => (file="us_nga_egm96_15.tif", resolution="15 arc-min",
                 description="EGM96 geoid undulation (WGS84→EGM96 height)"),
    :egm2008 => (file="us_nga_egm08_25.tif", resolution="2.5 arc-min",
                 description="EGM2008 geoid undulation (WGS84→EGM2008 height)"),
)

# --- Standard gravity and List (1968) latitude-dependent gravity / effective radius ---
# These are the closed forms used by NCL's `gp2gmh`. The equatorial constant 6378137 m is
# the WGS84 semi-major axis. The WGS84 Somigliana normal-gravity form
#   γ(φ) = γ_e·(1 + k·sin²φ) / √(1 − e²·sin²φ),  γ_e = 9.7803253359 m/s²,
#   k = 0.00193185265241,  e² = 0.00669437999013
# is the geodetic basis; the List-1968 series below is the standard numerical implementation.
const _G0            = 9.80665      # standard gravity, m/s² (geopotential-height definition)
const _GRAV_EQUATOR  = 9.80616      # List 1968 surface-gravity coefficient, m/s²
const _GRAV_C1       = 0.0026373    # List 1968 cos(2φ) coefficient
const _GRAV_C2       = 0.0000059    # List 1968 cos²(2φ) coefficient
const _EARTH_A       = 6378137.0    # WGS84 semi-major axis, m
const _RADIUS_B0     = 1.006803     # List 1968 effective-radius constant
const _RADIUS_B1     = 0.006706     # List 1968 sin²φ coefficient

"""
    _surface_gravity(latitude) -> gravity (m/s²)

Latitude-dependent surface gravity (List 1968), `latitude` in degrees. Broadcast-friendly.
"""
function _surface_gravity(latitude)
    ϕ = @. deg2rad(Float64(latitude))
    return @. _GRAV_EQUATOR * (1 - _GRAV_C1 * cos(2ϕ) + _GRAV_C2 * cos(2ϕ)^2)
end

"""
    _effective_earth_radius(latitude) -> radius (m)

Latitude-dependent effective Earth radius (List 1968), `latitude` in degrees.
Broadcast-friendly.
"""
function _effective_earth_radius(latitude)
    ϕ = @. deg2rad(Float64(latitude))
    return @. _EARTH_A / (_RADIUS_B0 - _RADIUS_B1 * sin(ϕ)^2)
end

"""
    _load_geoid(model; cache_path=nothing, force_download=false, verbose=true) -> Raster

Open a geoid grid as a **lazy** `Raster` streamed from the PROJ CDN via GDAL `/vsicurl/`
(HTTP byte-range reads; the grid is never downloaded in full). `model` is a key of
[`GEOID_MODELS`](@ref) (`:egm96`, `:egm2008`). `cache_path`/`force_download` are accepted
for API symmetry with the other invariant loaders but are unused (the grid is streamed,
not cached to disk).
"""
function _load_geoid(model::Symbol; cache_path=nothing, force_download::Bool=false, verbose::Bool=true)
    haskey(GEOID_MODELS, model) ||
        throw(ArgumentError("Unknown geoid model $(repr(model)). " *
                            "Available: $(join(sort(collect(keys(GEOID_MODELS))), ", "))"))
    _configure_gdal_http()   # shared with the Copernicus DEM loader (CA bundle + retries)
    info = GEOID_MODELS[model]
    url = "/vsicurl/$(_PROJ_CDN_BASE)/$(info.file)"
    verbose && @info "Streaming $(model) geoid ($(info.description), $(info.resolution))" source="$(_PROJ_CDN_BASE)/$(info.file)"
    return Raster(url; lazy=true)
end

"""
    geoid_undulation(latitude, longitude; geoid=:egm96, geoid_raster=nothing,
                     cache_path=nothing, force_download=false, verbose=true)

Geoid undulation `N` (m) — height of the geoid above the WGS84 ellipsoid — at the given
`latitude`/`longitude` (degrees). Positive where the geoid is above the ellipsoid.

`latitude`/`longitude` may be scalars or arrays (evaluated element-wise). Longitudes are
normalized to the grid's −180…180°E convention, so the ERA5-Land 0–360°E convention is
accepted transparently.

# Keyword arguments
- `geoid::Symbol=:egm96`: geoid model, a key of [`GEOID_MODELS`](@ref).
- `geoid_raster=nothing`: a pre-loaded geoid `Raster` to sample (skips streaming). When
  `nothing`, the grid is streamed once and read into memory for sampling.
- `cache_path`, `force_download`, `verbose`: forwarded to [`_load_geoid`](@ref).

# Returns
- `N` (m) as a scalar or array matching the input shape.

# References
- Lemoine et al. (1998), *The Development of the Joint NASA GSFC and NIMA Geopotential
  Model EGM96*, NASA/TP-1998-206861.
- Pavlis et al. (2012), *The development and evaluation of EGM2008*, JGR 117, B04406.
- Grids: OSGeo/PROJ-data (US NGA), served from https://cdn.proj.org.
"""
function geoid_undulation(latitude, longitude; geoid::Symbol=:egm96, geoid_raster=nothing,
                          cache_path=nothing, force_download::Bool=false, verbose::Bool=true)
    # Load and materialize the (small) geoid grid once; sampling many points against an
    # in-memory raster avoids a byte-range round-trip per point.
    g = geoid_raster === nothing ?
        read(_load_geoid(geoid; cache_path=cache_path, force_download=force_download, verbose=verbose)) :
        geoid_raster
    return _sample_geoid.(Ref(g), Float64.(latitude), _wrap_longitude.(Float64.(longitude)))
end

# Nearest-neighbour sample of the geoid grid at one (lat, lon) point.
function _sample_geoid(g, lat, lon)
    return Float64(g[X(Near(lon)), Y(Near(lat))])
end

"""
    geopotential2height(geopotential, latitude, longitude; height_reference=:wgs84,
                        ellipsoid=:wgs84, geoid=:egm96, cache_path=nothing,
                        force_download=false, verbose=true)
    geopotential2height(z::AbstractRaster; kwargs...)

Convert geopotential Φ (m² s⁻²) to geometric height (m) — above the WGS84 ellipsoid, or
above the geoid, per `height_reference`.

The conversion accounts for the latitude variation of gravity and Earth radius. The
`height_reference` keyword selects the output vertical datum:

- **ellipsoidal** (`:wgs84`, default) — adds the geoid undulation `N` streamed from the
  PROJ CDN, giving a true height above the WGS84 ellipsoid.
- **orthometric / geoid-referenced** (`:orthometric`) — returns the geometric height above
  the geoid directly (no `N` added, no network access). Geopotential height is inherently a
  height above mean sea level, so this is the natural frame for comparison with sea-level /
  geoid-referenced elevation products such as the Copernicus GLO-30 DEM (published as heights
  above the EGM2008 geoid). Note this cannot resolve which geoid model applies — ERA5's own
  datum is neither EGM96 nor EGM2008 — so only the single `:orthometric` option is offered
  rather than model-specific aliases that would imply more precision than exists.

Either way this is more accurate than the crude constant-gravity `Φ / 9.80665`.

Steps (`ϕ` = latitude, `λ` = longitude):
1. Geopotential height `H = Φ / g₀`, `g₀ = 9.80665` m s⁻².
2. Surface gravity `g(ϕ) = 9.80616·(1 − 0.0026373·cos2ϕ + 0.0000059·cos²2ϕ)` (List 1968).
3. Effective radius `R(ϕ) = 6378137 / (1.006803 − 0.006706·sin²ϕ)` (List 1968).
4. Geometric (orthometric) height above the geoid `Z = R·H / ((g/g₀)·R − H)`.
5. If `height_reference=:wgs84`: ellipsoidal height `h = Z + N(ϕ, λ)`;
   if `height_reference=:orthometric`: `h = Z`.

# Arguments
- `geopotential`: geopotential Φ in m² s⁻² (scalar or array). E.g. ERA5-Land `:z`.
- `latitude`, `longitude`: degrees (scalar or array; must match `geopotential`'s shape when
  arrays). Longitude may be −180…180 or 0–360°E.

For the `AbstractRaster` method, pass a geopotential raster (e.g.
`climate_model_invariant(parameter=:z)`); latitude/longitude are taken from its `Y`/`X`
dimensions and a `Raster` of heights (m) is returned.

# Keyword arguments
- `height_reference::Symbol=:wgs84`: output vertical datum. `:wgs84` (ellipsoidal) adds the
  geoid undulation; `:orthometric` returns the height above the geoid (≈ mean sea level)
  directly. Only a single geoid-referenced option is offered because this conversion cannot
  resolve which geoid model applies.
- `ellipsoid::Symbol=:wgs84`: reference ellipsoid for the ellipsoidal case. Only `:wgs84` is
  currently supported.
- `geoid::Symbol=:egm96`: geoid model used to stream `N` for the ellipsoidal case, a key of
  [`GEOID_MODELS`](@ref).
- `geoid_raster=nothing`: a pre-loaded geoid `Raster` to sample (skips streaming).
- `cache_path`, `force_download`, `verbose`: forwarded to the geoid loader.

# Returns
- Height (m) as a scalar/array matching the inputs, or a `Raster` for the raster method:
  ellipsoidal for `height_reference=:wgs84`, orthometric otherwise.

# Examples
```julia
# Summit, Greenland (72.58°N, −38.46°E), z ≈ 3216 m orographic:
h = geopotential2height(3216 * 9.80665, 72.58, -38.46)   # ellipsoidal ≈ orographic + N

# Orthometric height (≈ MSL; comparable to sea-level DEMs, no geoid streamed):
h_msl = geopotential2height(3216 * 9.80665, 72.58, -38.46; height_reference=:orthometric)

# Directly from ERA5-Land geopotential:
z = climate_model_invariant(parameter=:z)                # lazy Raster (m² s⁻²)
h = geopotential2height(read(z[X = 335 .. 347, Y = 63 .. 67]))   # ellipsoidal height, m
```

# References
- List, R.J. (1968), *Smithsonian Meteorological Tables*, 6th ed. (gravity & radius series;
  as implemented in NCL's `gp2gmh`).
- EGM96 / EGM2008 geoid grids: OSGeo/PROJ-data (US NGA), https://cdn.proj.org.
"""
function geopotential2height(geopotential, latitude, longitude; height_reference::Symbol=:wgs84,
                             ellipsoid::Symbol=:wgs84, geoid::Symbol=:egm96, geoid_raster=nothing,
                             cache_path=nothing, force_download::Bool=false, verbose::Bool=true)
    Z = _geopotential_to_geometric(geopotential, latitude)
    if height_reference === :orthometric
        # Geopotential height is already referenced to the geoid (≈ MSL); no N added.
        return Z
    elseif height_reference === :wgs84
        ellipsoid === :wgs84 ||
            throw(ArgumentError("Unsupported ellipsoid $(repr(ellipsoid)); only :wgs84 is supported."))
        N = geoid_undulation(latitude, longitude; geoid=geoid, geoid_raster=geoid_raster,
                             cache_path=cache_path, force_download=force_download, verbose=verbose)
        return @. Z + N
    else
        throw(ArgumentError("Unsupported height_reference $(repr(height_reference)); use " *
                            ":wgs84 (ellipsoidal) or :orthometric (above geoid ≈ MSL)."))
    end
end

# Map a `climate_forcing` dataset symbol to the invariant-registry model symbol.
_invariant_model(model::Symbol) = model === :era5land ? :era5_land : model

"""
    surface_elevation(model, lat, lon; height_reference=:orthometric,
                      cache_path=nothing, verbose=false) -> Float64

Surface elevation (m) of `model`'s grid cell nearest `(lat, lon)`, derived from the model's
geopotential invariant (`:z`). Samples the cached geopotential raster once at the nearest cell and
converts it with [`geopotential2height`](@ref).

`height_reference=:orthometric` (default) returns the height above the geoid (≈ mean sea level),
the natural frame for comparison with sea-level DEMs; `:wgs84` returns the ellipsoidal height.

`model` accepts either the invariant symbol (`:era5_land`) or the `climate_forcing` dataset symbol
(`:era5land`). `cache_path` is forwarded to [`climate_model_invariant`](@ref) for the geopotential
file. The geopotential invariant is downloaded once and cached, then opened lazily, so a point
sample reads at most one chunk from local disk.
"""
function surface_elevation(model::Symbol, lat::Real, lon::Real;
                           height_reference::Symbol=:orthometric,
                           cache_path=nothing, verbose::Bool=false)
    z = climate_model_invariant(; model=_invariant_model(model), parameter=:z,
                                cache_path=cache_path, verbose=verbose)  # lazy Raster, native 0–359.9°E
    # The geopotential grid is native 0–360°E; map a negative longitude into that range before
    # the nearest-cell lookup. `geopotential2height(:orthometric)` needs no geoid sample, so the
    # longitude sign is otherwise immaterial here.
    Φ = Float64(z[X(Near(lon < 0 ? lon + 360 : lon)), Y(Near(lat))])
    return geopotential2height(Φ, Float64(lat), Float64(lon); height_reference=height_reference)
end

function geopotential2height(z::AbstractRaster; kwargs...)
    hasdim(z, X) && hasdim(z, Y) ||
        throw(ArgumentError("geopotential2height(::AbstractRaster) needs both X and Y dims."))
    # Reshape lon/lat so each varies along its own dimension regardless of axis order,
    # then broadcast over the full grid (extra dims, e.g. a singleton Ti, broadcast trivially).
    lon_grid = _axis_vector(z, X, lookup(z, X))
    lat_grid = _axis_vector(z, Y, lookup(z, Y))
    h = geopotential2height(parent(z), lat_grid, lon_grid; kwargs...)
    return rebuild(z; data=h)
end

# Reshape a coordinate vector to a 1×…×n×…×1 array that varies only along dimension `D`
# of raster `z`, so it broadcasts correctly against `parent(z)` for any dimension order.
function _axis_vector(z::AbstractRaster, D, coords)
    shape = ntuple(i -> i == dimnum(z, D) ? length(coords) : 1, ndims(z))
    return reshape(collect(coords), shape)
end

"""
    _geopotential_to_geometric(geopotential, latitude) -> Z (m)

Geopotential Φ (m² s⁻²) → geometric (orthometric) height above the geoid (m), latitude in
degrees. Steps 1–4 of [`geopotential2height`](@ref); no geoid undulation. Broadcast-friendly,
with `Float64` promotion inside the fused kernel.
"""
function _geopotential_to_geometric(geopotential, latitude)
    H = @. Float64(geopotential) / _G0            # geopotential height (m)
    g = _surface_gravity(latitude)                 # m/s²
    R = _effective_earth_radius(latitude)          # m
    return @. R * H / ((g / _G0) * R - H)
end
