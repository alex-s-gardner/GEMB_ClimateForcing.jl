"""
Copernicus GLO-30 (30 m) global Digital Elevation Model — a **lazy** reader for the
open-access copy the Alaska Satellite Facility (ASF) points to, hosted as
Cloud-Optimized GeoTIFFs (COGs) in the AWS Open Data bucket
`copernicus-dem-30m`.

The DEM is tiled into 1°×1° COGs anchored at each tile's south-west corner:

    Copernicus_DSM_COG_10_{N|S}{lat:02d}_00_{E|W}{lon:03d}_00_DEM/…_DEM.tif

Because COGs support HTTP byte-range reads, GDAL's `/vsicurl/` driver fetches only
the requested spatial window (plus internal overviews) — the tiles are **never**
downloaded in full. A given lat/lon box is served by opening the covering tiles
(consulting the published tile index, since ocean tiles do not exist), mosaicking
them with an on-the-fly GDAL VRT, and returning a lazy `Raster` cropped to the box.

Reached through [`climate_model_invariant`](@ref) with `model=:copernicus_dem_30m`.
"""

using Rasters
using ArchGDAL          # loads Rasters' GDAL backend (RastersArchGDALExt)
import NetworkOptions   # Julia's CA bundle, so GDAL's curl can verify TLS
import Downloads

# AWS Open Data bucket (the same collection ASF provides open access to). Anonymous
# HTTPS, no authentication required.
const _COPERNICUS_DEM_30M_BASE = "https://copernicus-dem-30m.s3.amazonaws.com"

# Authoritative list of the ~26,450 tiles that actually exist (ocean tiles are absent).
const _COPERNICUS_DEM_30M_TILELIST_URL = "$(_COPERNICUS_DEM_30M_BASE)/tileList.txt"

# GDAL's bundled libcurl does not know where to find CA certificates on macOS, which
# makes `/vsicurl/` TLS handshakes fail intermittently ("unable to get local issuer
# certificate"). Point it at the CA bundle Julia ships, and enable HTTP retries so a
# transient S3 hiccup is retried rather than surfaced as a hard error. Idempotent;
# called at the start of every load.
function _configure_gdal_http()
    ca = NetworkOptions.ca_roots_path()
    if ca !== nothing && isfile(ca)
        ArchGDAL.setconfigoption("GDAL_HTTP_CAINFO", ca)
        ArchGDAL.setconfigoption("CURL_CA_BUNDLE", ca)
    end
    ArchGDAL.setconfigoption("GDAL_HTTP_MAX_RETRY", "5")
    ArchGDAL.setconfigoption("GDAL_HTTP_RETRY_DELAY", "1")
    # Don't probe the bucket directory on open — we already know the exact object key.
    ArchGDAL.setconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    return nothing
end

"""
    _copernicus_dem_tile_id(lat_sw, lon_sw) -> String

Tile directory/file stem for the 1°×1° tile whose south-west corner is the integer
`(lat_sw, lon_sw)` in degrees. Latitude is `N`/`S` zero-padded to 2 digits, longitude
`E`/`W` zero-padded to 3 digits, each followed by `_00`.

```jldoctest
julia> GEMB_ClimateForcing._copernicus_dem_tile_id(59, -150)
"Copernicus_DSM_COG_10_N59_00_W150_00_DEM"
```
"""
function _copernicus_dem_tile_id(lat_sw::Integer, lon_sw::Integer)
    ns = lat_sw < 0 ? "S" : "N"
    ew = lon_sw < 0 ? "W" : "E"
    latstr = lpad(abs(lat_sw), 2, '0')
    lonstr = lpad(abs(lon_sw), 3, '0')
    return "Copernicus_DSM_COG_10_$(ns)$(latstr)_00_$(ew)$(lonstr)_00_DEM"
end

"""
    _copernicus_dem_tile_url(id) -> String

Full HTTPS URL of a tile's COG GeoTIFF given its stem `id`.
"""
_copernicus_dem_tile_url(id::AbstractString) =
    "$(_COPERNICUS_DEM_30M_BASE)/$(id)/$(id).tif"

"""
    _copernicus_dem_vsicurl(id) -> String

GDAL `/vsicurl/` source string for a tile, so GDAL reads it via HTTP byte ranges.
"""
_copernicus_dem_vsicurl(id::AbstractString) = "/vsicurl/" * _copernicus_dem_tile_url(id)

"""
    _copernicus_dem_tile_dir(cache_path) -> String

Where cached tiles live: a `tiles/` subdirectory, kept apart from `tileList.txt` and the
generated `.vrt` files so a tile sweep never has to pattern-match around them.
"""
_copernicus_dem_tile_dir(cache_path::AbstractString) = joinpath(cache_path, "tiles")

"""
    _copernicus_dem_remote_size(url) -> Union{Int,Nothing}

`Content-Length` for `url`, or `nothing` if the server does not report one.

Used only to detect a truncated download. A partial GeoTIFF is the failure mode worth guarding:
GDAL reports it as an unrelated parse error, and once written it would be treated as cached and
reused on every later call.
"""
function _copernicus_dem_remote_size(url::AbstractString)
    try
        r = Downloads.request(url; method = "HEAD", throw = false)
        r isa Downloads.Response || return nothing
        for (k, v) in r.headers
            lowercase(k) == "content-length" && return tryparse(Int, strip(v))
        end
    catch
        return nothing
    end
    return nothing
end

"""
    _copernicus_dem_cache_tile(id; cache_path, force, verbose) -> String

Local path of tile `id`, downloading it once if absent.

This is what makes caching *persistent*. GDAL has no on-disk cache for `/vsicurl/` — its block
cache lives in the process — so repeated point lookups across sessions otherwise re-issue HTTP
range requests for the same bytes. Fetching the whole tile trades one 19–40 MB download for
unlimited local reads.

Downloads to `.part` and renames, so an interrupted fetch cannot leave a short file that later
runs would trust; the length is checked against `Content-Length` when the server reports one.
"""
function _copernicus_dem_cache_tile(id::AbstractString; cache_path::AbstractString,
                                    force::Bool = false, verbose::Bool = true)
    dir = _copernicus_dem_tile_dir(cache_path)
    mkpath(dir)
    local_path = joinpath(dir, id * ".tif")
    (!force && isfile(local_path)) && return local_path

    url = _copernicus_dem_tile_url(id)
    expected = _copernicus_dem_remote_size(url)
    verbose && @info "Downloading Copernicus DEM tile" id size_MB = isnothing(expected) ?
        missing : round(expected / 1e6; digits = 1)
    tmp = local_path * ".part"
    try
        Downloads.download(url, tmp)
        got = filesize(tmp)
        if !isnothing(expected) && got != expected
            error("truncated download: got $(got) bytes, expected $(expected)")
        end
        mv(tmp, local_path; force = true)
    catch e
        isfile(tmp) && rm(tmp; force = true)
        error("Failed to download Copernicus DEM tile $(id) from $(url):\n$(e)")
    end
    return local_path
end

"""
    _copernicus_dem_tile_index(; cache_path, force) -> Set{String}

Download the published `tileList.txt` once (cached under `cache_path`), and return the
set of tile stems that actually exist. Reuses the download-to-`.part`-then-`mv` idiom.
"""
function _copernicus_dem_tile_index(; cache_path::String, force::Bool, verbose::Bool=true)
    mkpath(cache_path)
    local_path = joinpath(cache_path, "tileList.txt")
    if force || !isfile(local_path)
        verbose && @info "Downloading Copernicus DEM tile index" url=_COPERNICUS_DEM_30M_TILELIST_URL
        tmp = local_path * ".part"
        try
            Downloads.download(_COPERNICUS_DEM_30M_TILELIST_URL, tmp)
            mv(tmp, local_path; force=true)
        catch e
            isfile(tmp) && rm(tmp; force=true)
            error("Failed to download Copernicus DEM tile index from " *
                  "$(_COPERNICUS_DEM_30M_TILELIST_URL):\n$e")
        end
    end
    index = Set{String}()
    for line in eachline(local_path)
        s = strip(line)
        isempty(s) || push!(index, String(s))
    end
    return index
end

# Tile geometry (derivable from the tile ID alone — no need to open any tile).
# Every tile has 3600 rows (1 arc-second latitude sampling). Longitude is sub-sampled
# in latitude bands, so the column count drops toward the poles. Verified against real
# tile headers across all bands.
const _COPERNICUS_DEM_ROWS = 3600
const _COPERNICUS_DEM_FINEST_COLS = 3600   # equatorial band; the common (finest) grid

"""
    _copernicus_dem_ncols(lat_sw) -> Int

Number of longitude columns in the 1°×1° tile whose south-west corner latitude is the
integer `lat_sw`. Copernicus GLO-30 sub-samples longitude poleward; the band is keyed by
the tile's equatorward edge magnitude.
"""
function _copernicus_dem_ncols(lat_sw::Integer)
    a = lat_sw < 0 ? abs(lat_sw) - 1 : abs(lat_sw)   # equatorward edge magnitude
    a < 50 && return 3600
    a < 60 && return 2400
    a < 70 && return 1800
    a < 80 && return 1200
    a < 85 && return 720
    return 360
end

"""
    _copernicus_dem_parse_corner(id) -> (lat_sw::Int, lon_sw::Int)

Inverse of [`_copernicus_dem_tile_id`](@ref): the integer SW-corner degrees of a tile.
"""
function _copernicus_dem_parse_corner(id::AbstractString)
    m = match(r"_(N|S)(\d{2})_00_(E|W)(\d{3})_00_DEM$", id)
    m === nothing && throw(ArgumentError("Unrecognized Copernicus DEM tile id: $(repr(id))"))
    lat = parse(Int, m.captures[2]) * (m.captures[1] == "S" ? -1 : 1)
    lon = parse(Int, m.captures[4]) * (m.captures[3] == "W" ? -1 : 1)
    return (lat, lon)
end

"""
    _copernicus_dem_extent(extent) -> NTuple{4,Float64}  # (xmin, xmax, ymin, ymax)

Normalize a user-supplied `extent` to `(xmin, xmax, ymin, ymax)` in degrees. Accepts an
`Extents.Extent` (with `X`/`Y` bounds) or a NamedTuple `(; X=(a,b), Y=(c,d))`.
"""
function _copernicus_dem_extent(extent)
    x = extent.X
    y = extent.Y
    xmin, xmax = float(minimum(x)), float(maximum(x))
    ymin, ymax = float(minimum(y)), float(maximum(y))
    (-180.0 <= xmin && xmax <= 180.0) ||
        throw(ArgumentError("Copernicus DEM longitude must be in [-180, 180]; got X=$(x)"))
    (-90.0 <= ymin && ymax <= 90.0) ||
        throw(ArgumentError("Copernicus DEM latitude must be in [-90, 90]; got Y=$(y)"))
    return (xmin, xmax, ymin, ymax)
end

"""
    _copernicus_dem_tiles_for_extent(xmin, xmax, ymin, ymax, index) -> Vector{String}

Tile stems covering the box, keeping only those present in `index`. Errors if the box
is entirely over ocean (no tiles exist).
"""
function _copernicus_dem_tiles_for_extent(xmin, xmax, ymin, ymax, index::Set{String})
    lon0 = floor(Int, xmin)
    lon1 = floor(Int, xmax)   # SW corner of the tile containing xmax
    lat0 = floor(Int, ymin)
    lat1 = floor(Int, ymax)
    ids = String[]
    for lat_sw in lat0:lat1, lon_sw in lon0:lon1
        id = _copernicus_dem_tile_id(lat_sw, lon_sw)
        id in index && push!(ids, id)
    end
    isempty(ids) && throw(ArgumentError(
        "No Copernicus DEM tiles exist for extent X=($xmin, $xmax), Y=($ymin, $ymax) " *
        "— the region is outside the DEM's land coverage (ocean tiles are not published)."))
    return ids
end

"""
    _copernicus_dem_write_vrt(ids, path)

Author a GDAL VRT XML at `path` that mosaics the given tile `ids` onto a single common
grid, **purely from the tiles' analytically-derived geometry** — no tile is opened. Each
tile is placed on the finest (equatorial, 3600 cols/°) grid; coarser high-latitude tiles
are upsampled by GDAL on read (nearest), which reproduces direct tile reads exactly.
Uncovered cells (absent ocean tiles) read as 0, matching the dataset's "height 0 over
ocean" convention. Sources are `/vsicurl/` URLs, so reads confine to HTTP byte ranges.
"""
function _copernicus_dem_write_vrt(ids::Vector{String}, path::String,
                                   sources::Vector{String} = _copernicus_dem_vsicurl.(ids))
    length(sources) == length(ids) || throw(DimensionMismatch(
        "got $(length(sources)) sources for $(length(ids)) tiles"))
    corners = _copernicus_dem_parse_corner.(ids)
    west  = minimum(c[2] for c in corners)          # min SW-corner longitude
    east  = maximum(c[2] for c in corners) + 1       # max tile east edge
    south = minimum(c[1] for c in corners)           # min SW-corner latitude
    north = maximum(c[1] for c in corners) + 1       # max tile north edge

    cols = _COPERNICUS_DEM_FINEST_COLS
    xsize = (east - west) * cols
    ysize = (north - south) * _COPERNICUS_DEM_ROWS
    dx =  1 / cols
    dy = -1 / _COPERNICUS_DEM_ROWS

    open(path, "w") do io
        println(io, """<VRTDataset rasterXSize="$xsize" rasterYSize="$ysize">""")
        println(io, "  <SRS>EPSG:4326</SRS>")
        println(io, "  <GeoTransform>$(float(west)), $dx, 0, $(float(north)), 0, $dy</GeoTransform>")
        println(io, """  <VRTRasterBand dataType="Float32" band="1">""")
        for (src, (lat_sw, lon_sw)) in zip(sources, corners)
            nc = _copernicus_dem_ncols(lat_sw)
            dxoff = (lon_sw - west) * cols
            dyoff = (north - (lat_sw + 1)) * _COPERNICUS_DEM_ROWS
            println(io, "    <ComplexSource>")
            println(io, """      <SourceFilename relativeToVRT="0">$src</SourceFilename>""")
            println(io, "      <SourceBand>1</SourceBand>")
            println(io, """      <SourceProperties RasterXSize="$nc" RasterYSize="$(_COPERNICUS_DEM_ROWS)" DataType="Float32"/>""")
            println(io, """      <SrcRect xOff="0" yOff="0" xSize="$nc" ySize="$(_COPERNICUS_DEM_ROWS)"/>""")
            println(io, """      <DstRect xOff="$dxoff" yOff="$dyoff" xSize="$cols" ySize="$(_COPERNICUS_DEM_ROWS)"/>""")
            println(io, "    </ComplexSource>")
        end
        println(io, "  </VRTRasterBand>")
        println(io, "</VRTDataset>")
    end
    return path
end

"""
    _load_copernicus_dem_30m(extent; cache_path, force_download) -> Raster

Load the Copernicus GLO-30 DEM as a **lazy** `Raster` (EPSG:4326, `X`/`Y` in
−180…180 / −90…90). `extent === nothing` loads the **full data extent** (a global lazy
mosaic of every published land tile); otherwise only the covering tiles are used.

The mosaic is assembled by authoring a GDAL VRT analytically from the tile list — **no
tiles are opened** to build it — so even the global default is cheap to construct. A
single covering tile is opened directly at native resolution. Reads go through GDAL
`/vsicurl/`, so nothing is fetched until the raster is indexed/cropped/`read`.
"""
function _load_copernicus_dem_30m(extent; cache_path::String, force_download::Bool,
                                  cache_tiles::Bool=false, max_concurrent_downloads::Integer=4,
                                  verbose::Bool=true)
    # Asking for a persistent cache while the location is a temp directory is refused rather
    # than honored: the tiles would be written, reported as cached, and then reaped, so the
    # caller would pay the download repeatedly while believing they had paid it once.
    if cache_tiles && startswith(abspath(cache_path), abspath(tempdir()))
        throw(ArgumentError(
            "cache_tiles=true needs a cache location that survives the OS clearing " *
            "$(tempdir()), but cache_path resolves to $(cache_path). Set " *
            "ENV[\"GEMB_CACHE_PATH\"] to a bulk-storage volume (not a home directory, which " *
            "is commonly a small SSD with a quota), or pass `cache_path` explicitly. Leave " *
            "cache_tiles=false to read over /vsicurl/ without caching anything."))
    end
    _configure_gdal_http()
    index = _copernicus_dem_tile_index(; cache_path=cache_path, force=force_download, verbose=verbose)

    if extent === nothing
        # Full data extent: every published tile.
        ids = sort!(collect(index))
        xmin, xmax, ymin, ymax = -180.0, 180.0, -90.0, 90.0
        verbose && @info "Copernicus DEM: full data extent ($(length(ids)) tiles)"
    else
        xmin, xmax, ymin, ymax = _copernicus_dem_extent(extent)
        ids = _copernicus_dem_tiles_for_extent(xmin, xmax, ymin, ymax, index)
        verbose && @info "Copernicus DEM: $(length(ids)) tile(s) cover the requested extent"
    end

    # `cache_tiles` fetches each covering tile in full and points the sources at local files,
    # so subsequent reads — in this session or any later one — touch no network. Off by default:
    # the /vsicurl/ path reads only the bytes it needs, which is the cheaper choice for a
    # one-off crop and the only tractable one for a continental extent.
    if cache_tiles
        isempty(ids) || verbose && @info "Copernicus DEM: caching $(length(ids)) tile(s) locally" dir=_copernicus_dem_tile_dir(cache_path)
        # Fanned out rather than looped: each tile is an independent 19-40 MB HTTPS GET, so
        # serially the wall time is the sum of the transfers instead of roughly the longest.
        # `_run_concurrent_jobs` is the same primitive `mcd43a3_granules` uses for granule
        # downloads, with no submission stagger — these are bandwidth-bound, not queued.
        sources = Vector{String}(undef, length(ids))
        _run_concurrent_jobs(length(ids), max_concurrent_downloads, 0;
                             verbose=false) do k, _gate
            sources[k] = _copernicus_dem_cache_tile(ids[k]; cache_path=cache_path,
                                                    force=force_download, verbose=verbose)
        end
    else
        sources = _copernicus_dem_vsicurl.(ids)
    end

    if length(ids) == 1
        raster = Raster(sources[1]; lazy=true)
    else
        # Analytical VRT (no tile opens) referencing the sources. Unlike `Rasters.mosaic`, this
        # stays lazy and never materializes the inputs. The mode is part of the filename because
        # the same tile set yields different VRT content local vs remote.
        mkpath(cache_path)
        tag = cache_tiles ? "local" : "vsicurl"
        vrt_path = joinpath(cache_path, "cop30_$(tag)_" * string(hash(ids); base=16) * ".vrt")
        _copernicus_dem_write_vrt(ids, vrt_path, sources)
        raster = Raster(vrt_path; lazy=true)
    end

    # Crop to the requested box (stays lazy) and attach descriptive metadata.
    cropped = view(raster, X=(xmin .. xmax), Y=(ymin .. ymax))
    meta = Dict{String,Any}(
        "source"        => "Copernicus GLO-30 DEM (AWS Open Data / ASF)",
        "base_url"      => _COPERNICUS_DEM_30M_BASE,
        "resolution_m"  => 30,
        "n_tiles"       => length(ids),
        "cached_tiles"  => cache_tiles,
        "extent"        => (X=(xmin, xmax), Y=(ymin, ymax)),
        "units"         => "m",
        "long_name"     => "surface elevation above EGM2008 geoid",
    )
    # Only embed the full tile list when small; a global mosaic has tens of thousands.
    length(ids) <= 64 && (meta["tiles"] = ids)
    return rebuild(cropped; metadata=meta)
end
