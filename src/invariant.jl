"""
Time-invariant (static) climate-model parameters — land–sea mask, geopotential,
vegetation/soil/lake fields — served as ready-made global files rather than in the
time-series ARCO Zarr stores.

For ERA5-Land these are the "auxiliary land invariant parameters" published by
ECMWF, already interpolated to the ERA5-Land 0.1°×0.1° grid, and distributed as
single-field NetCDF4 files attached to the ERA5-Land data-documentation page:
https://confluence.ecmwf.int/spaces/CKB/pages/140385202/ERA5-Land+data+documentation

They are downloaded once, cached on disk, and returned as **lazy** `Raster`s (data
is only read from disk when the raster is indexed/cropped/collected).
"""

using Rasters
using NCDatasets   # loads Rasters' NetCDF backend (RastersNCDatasetsExt)
using DimensionalData
import Downloads

# Base URL of the ERA5-Land invariant-parameter attachments (page id 140385202).
const _ERA5_LAND_INVARIANT_BASE =
    "https://confluence.ecmwf.int/download/attachments/140385202"

# Invariant parameters available for ERA5-Land, keyed by GRIB shortName, mapping to
# (NetCDF filename, paramId, description). NetCDF is preferred so we can read it in
# pure Julia (NCDatasets) without a GRIB decoder. All are on the 0.1°×0.1° grid,
# longitude 0–359.9°E, latitude 90→−90°N (descending), with a singleton time axis.
const ERA5_LAND_INVARIANT_PARAMETERS = Dict{Symbol,NamedTuple{(:file, :paramId, :description),Tuple{String,Int,String}}}(
    :lsm => (file="lsm_1279l4_0.1x0.1.grb_v4_unpack.nc",      paramId=172,    description="Land-sea mask (0–1 land fraction)"),
    :z   => (file="geo_1279l4_0.1x0.1.grib2_v4_unpack.nc",    paramId=129,    description="Geopotential (m² s⁻²); orography = z / 9.80665"),
    :cl  => (file="clake.nc",                                 paramId=26,     description="Lake cover (0–1)"),
    :dl  => (file="dl.nc",                                    paramId=228007, description="Lake total depth (m)"),
    :cvl => (file="cvl.nc",                                   paramId=27,     description="Low vegetation cover (0–1)"),
    :cvh => (file="cvh.nc",                                   paramId=28,     description="High vegetation cover (0–1)"),
    :glm => (file="cicecap_v015_1279_4_regMIR.nc",            paramId=260294, description="Glacier mask"),
    :slt => (file="slt.nc",                                   paramId=43,     description="Soil type"),
    :tvl => (file="tvl.nc",                                   paramId=29,     description="Type of low vegetation"),
    :tvh => (file="tvh.nc",                                   paramId=30,     description="Type of high vegetation"),
)

# Registry of invariant parameters per supported model. Add new models here.
const _INVARIANT_REGISTRY = Dict{Symbol,typeof(ERA5_LAND_INVARIANT_PARAMETERS)}(
    :era5_land => ERA5_LAND_INVARIANT_PARAMETERS,
)

const _INVARIANT_BASE_URL = Dict{Symbol,String}(
    :era5_land => _ERA5_LAND_INVARIANT_BASE,
)

"""
    _default_invariant_cache(model::Symbol) -> String

Default on-disk cache directory for a model's downloaded invariant files.
"""
_default_invariant_cache(model::Symbol) =
    joinpath(tempdir(), "GEMB_ClimateForcing", "invariant", string(model))

"""
    _download_invariant(model, parameter; cache_path, force) -> String

Ensure the NetCDF file for `parameter` is present in `cache_path`, downloading it
once if missing (or if `force`), and return its local path.
"""
function _download_invariant(model::Symbol, parameter::Symbol; cache_path::String, force::Bool)
    params = _INVARIANT_REGISTRY[model]
    info = params[parameter]
    mkpath(cache_path)
    local_path = joinpath(cache_path, info.file)
    url = string(_INVARIANT_BASE_URL[model], "/", info.file)
    if force || !isfile(local_path)
        println("  Downloading $(parameter) ($(info.description))")
        println("    $(url)")
        tmp = local_path * ".part"
        try
            Downloads.download(url, tmp)
            mv(tmp, local_path; force=true)
        catch e
            isfile(tmp) && rm(tmp; force=true)
            error("Failed to download invariant parameter $(parameter) for $(model) from $(url):\n$e")
        end
    else
        println("  Using cached $(parameter): $(local_path)")
    end
    return local_path
end

"""
    _open_invariant_raster(path) -> Raster

Open a single invariant NetCDF file as a lazy 2-D (`X`, `Y`) `Raster`. The stored
files carry a singleton `time` axis; it is dropped so the result is purely spatial.
"""
function _open_invariant_raster(path::String)
    # lazy=true keeps the data disk-backed; only metadata/coords are read now.
    r = Raster(path; lazy=true)
    # Drop the singleton time dimension if present, keeping the raster lazy.
    return hasdim(r, Ti) ? view(r, Ti(1)) : r
end

"""
    climate_model_invariant(; model=:era5_land, parameter=nothing,
                              cache_path=nothing, force_download=false)

Load time-invariant (static) parameters for a climate model as **lazy** `Raster`s.

Invariant fields (land–sea mask, geopotential/orography, vegetation, soil, lake,
glacier) are not part of the time-series ARCO Zarr stores; ECMWF distributes them
as ready-made global files on the ERA5-Land grid. This function downloads the
requested file(s) once, caches them on disk, and opens them lazily — no array data
is read until the returned raster is indexed, cropped, or `read`/`collect`ed.

# Keyword arguments
- `model::Symbol=:era5_land`: source model. Currently only `:era5_land`.
- `parameter::Union{Symbol,Nothing}=nothing`: which invariant to load, by GRIB
  shortName.
  - `nothing` (default) — load **all** available parameters and return a lazy
    `RasterStack`.
  - a single `Symbol` (e.g. `:lsm`, `:z`) — return a single lazy `Raster`.
  Available for `:era5_land`: `:lsm` (land-sea mask), `:z` (geopotential),
  `:cl` (lake cover), `:dl` (lake depth), `:cvl`/`:cvh` (low/high vegetation
  cover), `:tvl`/`:tvh` (low/high vegetation type), `:slt` (soil type),
  `:glm` (glacier mask). See [`ERA5_LAND_INVARIANT_PARAMETERS`](@ref).
- `cache_path::Union{String,Nothing}=nothing`: directory for the downloaded files.
  Defaults to a per-model folder under `tempdir()`.
- `force_download::Bool=false`: re-download even if a cached file exists.

# Returns
- A lazy `Raster` (single `parameter`) or a lazy `RasterStack` (`parameter =
  nothing`), with `X` (longitude) and `Y` (latitude) dimensions and CRS EPSG:4326.

!!! note "Longitude convention"
    The ERA5-Land invariant grid uses **0–359.9°E** longitude (not −180…180) and
    **descending** latitude (90→−90°N), matching the source files. Select longitude
    in that 0–360 convention; the `X = a .. b` / `Y = a .. b` interval selector takes
    `min .. max` regardless of the axis order, so pass `Y = 63 .. 67` for a box near
    Iceland at `X = 335 .. 347` (≈ −25…−13°E). Geopotential `z` is in m² s⁻²; divide
    by `9.80665` for orography (m).

# Examples
```julia
# Land-sea mask as a lazy Raster; nothing is read until you crop/collect it.
lsm = climate_model_invariant(parameter=:lsm)
iceland = read(lsm[X = 335 .. 347, Y = 63 .. 67])   # fractional 0–1 land cover

# Orography (m) from geopotential.
z = climate_model_invariant(parameter=:z)
orography = z ./ 9.80665

# Everything as a lazy stack.
inv = climate_model_invariant()          # RasterStack with :lsm, :z, :cvl, ...
```
"""
function climate_model_invariant(;
    model::Symbol=:era5_land,
    parameter::Union{Symbol,Nothing}=nothing,
    cache_path::Union{String,Nothing}=nothing,
    force_download::Bool=false,
)
    haskey(_INVARIANT_REGISTRY, model) ||
        throw(ArgumentError("Unsupported model $(repr(model)) for invariant parameters. " *
                            "Supported: $(join(sort(collect(keys(_INVARIANT_REGISTRY))), ", "))"))
    params = _INVARIANT_REGISTRY[model]
    cache = isnothing(cache_path) ? _default_invariant_cache(model) : cache_path

    println("Loading $(model) invariant parameter(s) as lazy Raster(s)...")

    if parameter === nothing
        # Load every available parameter into a lazy RasterStack.
        names = sort(collect(keys(params)))
        rasters = map(names) do p
            path = _download_invariant(model, p; cache_path=cache, force=force_download)
            _open_invariant_raster(path)
        end
        stack = RasterStack(NamedTuple{Tuple(names)}(Tuple(rasters)))
        println("  ✓ Loaded $(length(names)) invariant parameters: $(join(names, ", "))")
        return stack
    else
        haskey(params, parameter) ||
            throw(ArgumentError("Unknown invariant parameter $(repr(parameter)) for $(model). " *
                                "Available: $(join(sort(collect(keys(params))), ", "))"))
        path = _download_invariant(model, parameter; cache_path=cache, force=force_download)
        r = _open_invariant_raster(path)
        println("  ✓ Loaded $(parameter) as lazy Raster $(size(r))")
        return r
    end
end
