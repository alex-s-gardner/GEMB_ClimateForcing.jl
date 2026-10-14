"""
    get_cds_api_key()

Get the CDS (Copernicus Climate Data Store) API key from environment variable or .cdsapirc file.

# Returns
- String containing the CDS API key

# Description
Searches for the API key in the following order:
1. Environment variable `CDS_API_KEY`
2. `~/.cdsapirc` file (standard location for CDS API credentials)

# Throws
- `ErrorException` if no API key is found

# Example
```julia
token = get_cds_api_key()
forcing_data = climate_forcing(:era5land, lat, lon; time_range=..., token=token)
```

# See Also
- Get your CDS API key from: https://cds.climate.copernicus.eu/api-how-to
"""
function get_cds_api_key()
    # First check environment variable
    if haskey(ENV, "CDS_API_KEY")
        return ENV["CDS_API_KEY"]
    end

    # Try reading from .cdsapirc file
    cdsapirc = joinpath(homedir(), ".cdsapirc")
    if isfile(cdsapirc)
        for line in readlines(cdsapirc)
            if startswith(line, "key:")
                return String(strip(split(line, ":", limit=2)[2]))
            end
        end
    end

    error("CDS API key not found. Set ENV[\"CDS_API_KEY\"] or create ~/.cdsapirc file")
end

"""
    _retry_transient(f; is_transient, context, attempts, base, max_backoff,
                     deadline=Inf, verbose=true, retry_after=_ -> 0.0,
                     on_retry=(attempt, backoff, err) -> nothing) -> f()

Call `f`, retrying while it fails *transiently* with exponential backoff from `base` to
`max_backoff`. Non-transient errors propagate on the first attempt, unchanged.

The transient/permanent split is supplied by the caller as `is_transient`, because it is
service-specific and getting it wrong is expensive in both directions: retrying a
permanent error produces an endless resubmit loop, while failing fast on a service hiccup
throws away however much work the run had already done.

`deadline` (an absolute `time()` value) bounds the retrying so a caller's own timeout is
still respected — and it, not `attempts`, is what should normally bind. The final failure
is rethrown rather than wrapped, so the message the user sees is the service's own.

`retry_after(err)` lets a server-supplied wait (HTTP `Retry-After`, e.g. on 429) override
the backoff curve upward; it knows better than we do. `on_retry` is called before each
sleep, for logging.

Shared by the CDS Retrieve client (`_cds_retry`) and the NASA Earthdata client
(`_earthdata_retry`), which differ only in their error taxonomy and constants.
"""
function _retry_transient(f; is_transient, context::AbstractString,
                          attempts::Integer, base::Real, max_backoff::Real,
                          deadline::Float64=Inf, verbose::Bool=true,
                          retry_after=_ -> 0.0,
                          on_retry=(attempt, backoff, err) -> nothing)
    attempt = 0
    while true
        attempt += 1
        try
            return f()
        catch err
            (is_transient(err) && attempt < attempts) || rethrow()
            backoff = min(base * 2.0^(attempt - 1), max_backoff)
            requested = retry_after(err)
            requested > 0 && (backoff = max(backoff, requested))
            time() + backoff >= deadline && rethrow()
            verbose && on_retry(attempt, backoff, err)
            sleep(backoff)
        end
    end
end

"""
    saturation_vapor_pressure(T_K::Real) -> Float64

Saturation vapor pressure (Pa) over liquid water at temperature `T_K` (Kelvin), in the form
the ECMWF IFS uses (Tetens' formula with the Buck (1981) constants for water):

    eₛ(T) = 611.21 · exp(17.502 · (T − 273.16) / (T − 32.19))

This is the only saturation curve in the package. ERA5 and ERA5-Land compute their 2 m
dewpoint "with respect to water to be consistent with WMO reporting practice" (IFS Cy45r1
Part IV §3.10.3), so `eₛ` over water at the dewpoint is the vapor pressure the reanalysis
itself implies, and relative humidity is defined against the same curve. Using one curve for
the dewpoint conversion, for relative humidity and for every adjustment is what makes the
identity adjustments (`Δz = 0`, `ΔT = 0`) exact: a ratio `eₛ(T′)/eₛ(T)` of the same function
is exactly 1 at `T′ = T`. Relative humidity against water stays below 1 in sub-freezing air
even where the air is saturated with respect to ice.

# References
- ECMWF (2018). IFS Documentation Cy45r1, Part IV: Physical Processes, Eq. 7.5 and §3.10.3.
- Buck, A. L. (1981). New equations for computing vapor pressure and enhancement factor.
  *J. Appl. Meteorol.* 20, 1527–1532.
"""
saturation_vapor_pressure(T_K::Real) =
    611.21 * exp(17.502 * (Float64(T_K) - 273.16) / (Float64(T_K) - 32.19))

"""
    dewpoint_to_vapor_pressure(T_dewpoint)

Vapor pressure (Pa) from dewpoint temperature (K): `e = eₛ(T_d)`, with `eₛ` the saturation
vapor pressure over water of [`saturation_vapor_pressure`](@ref) — the curve ERA5 and
ERA5-Land define their 2 m dewpoint against (IFS Cy45r1 Part IV §3.10.3).

Accepts a scalar or an array (broadcast); `Float32` input is promoted to `Float64`.
"""
dewpoint_to_vapor_pressure(T_dewpoint) = saturation_vapor_pressure.(T_dewpoint)

"""
    vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)

Relative humidity [%] from vapor pressure [Pa] and air temperature [K], against the
saturation vapor pressure over water of [`saturation_vapor_pressure`](@ref), clamped to
[0, 100]. The clamp makes this a display and fitting convenience; the forcing adjustments
never round-trip through it.
"""
function vapor_pressure_to_relative_humidity(vapor_pressure, temperature_air)
    relative_humidity = (vapor_pressure ./ saturation_vapor_pressure.(temperature_air)) .* 100.0
    return clamp.(relative_humidity, 0.0, 100.0)
end

"""
    relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)

Vapor pressure [Pa] from air temperature [K] and relative humidity [%] against the
saturation vapor pressure over water of [`saturation_vapor_pressure`](@ref).
"""
function relative_humidity_to_vapor_pressure(temperature_air, relative_humidity)
    return saturation_vapor_pressure.(temperature_air) .* (relative_humidity ./ 100.0)
end

# Vapor pressure carried from temperature `T` to `T′` at constant relative humidity:
# `e′ = e·eₛ(T′)/eₛ(T)`. The ratio form needs no relative humidity, so nothing is clamped: an
# input at or below saturation stays at or below saturation at `T′`, and `T′ == T` returns `e`
# bit for bit, since `eₛ(T)/eₛ(T)` is exactly 1.
_vapor_pressure_at(e, T, T′) = e * (saturation_vapor_pressure(T′) / saturation_vapor_pressure(T))

# Vapor pressure after the glacier decoupling cools the air from `T` to `T′` at fixed moisture:
# unchanged, except that cooled air holds at most saturation at its new temperature. Applied only
# where the air was cooled, so an unchanged step returns `e` bit for bit even when its input is
# supersaturated.
_cap_vapor_pressure(e, T, T′) = T′ < T ? min(e, saturation_vapor_pressure(T′)) : e

# ----------------------------------------------------------------------------
# Shared helpers for the forcing-adjustment functions
# (elevation_adjustment.jl, climate_adjustment.jl, glacier_adjustment.jl).
# ----------------------------------------------------------------------------

"""
    _rebuild_forcing(stack, new_metadata; changed_layers...) -> DimStack

Rebuild a forcing `DimStack` with new stack metadata, replacing only the named layers.

Every adjustment function changes a different subset of the seven forcing variables; this
carries the rest through untouched rather than making each call site re-enumerate them.
`rebuild` preserves each layer's own metadata, so broadcast results (which drop it) come
back correctly annotated without threading it by hand.

Values may be `DimArray`s or plain arrays — only the underlying data is used.
"""
function _rebuild_forcing(stack::DimStack, new_metadata; changed_layers...)
    changed = NamedTuple(changed_layers)
    unknown = setdiff(keys(changed), keys(stack))
    isempty(unknown) ||
        throw(ArgumentError("not a layer of this forcing stack: $(join(unknown, ", "))"))
    # Unwrapped after the merge, not before: `rebuild` wraps every layer in the stack's dimensions, and
    # `NamedTuple(stack)` hands back the carried layers as `DimArray`s. Unwrapping only the replacements
    # leaves those to be wrapped a second time, and a doubly wrapped layer's `parent` is a `DimArray`
    # rather than the values a consumer asked for — which reads as an array of numbers everywhere until
    # something dispatches on `Vector`.
    return rebuild(stack; data=map(parent, merge(NamedTuple(stack), changed)),
                   metadata=new_metadata)
end

"""
    _adjust_longwave(LW, e, T, e′, T′)

Recompute downwelling longwave irradiance [W/m²] for adjusted temperature `T′` and vapor
pressure `e′`, preserving the cloud/aerosol emissivity increment diagnosed from the input.

`Δε = LW/(σT⁴) − ε_cs(e, T)` is the departure of the observed irradiance from the
Konzelmann et al. (1994) clear-sky value; holding it fixed means only the clear-sky
response and the `σT⁴` scaling change. Bit-exact identity when `(e′, T′) == (e, T)`.

Δε is preserved *additively* (Glover 1999) and is deliberately left unclamped, so it carries
through unchanged even where the bulk emissivity `LW/(σT⁴)` exceeds 1 — which reanalysis does
produce under a surface inversion with a cloud base warmer than the 2 m temperature. Δε < 0
(ε_cs over-predicting) is likewise carried as-is.
"""
function _adjust_longwave(LW, e, T, e′, T′)
    ε_cs  = konzelmann_clear_sky_emissivity.(e, T)
    ε_cs′ = konzelmann_clear_sky_emissivity.(e′, T′)
    # `(ε_cs′ + Δε)·σT′⁴` rearranged so each term vanishes, or is multiplied by exactly 1, when
    # `(e′, T′) == (e, T)`: the identity adjustment then returns `LW` bit for bit.
    return @. LW * ((T′ / T)^2)^2 + (ε_cs′ - ε_cs) * _SIGMA_SB * T′^4
end

# Normalize longitude to the −180…180°E convention (accepts 0–360°E, or any multiple
# wrap: 540 → −180, −190 → 170). Used by the geoid grid, the glacier lookup table and the
# seam-aware nearest-longitude search.
_wrap_longitude(lon::Real) = mod(lon + 180.0, 360.0) - 180.0

"""
    _gemb_cache_root() -> String

Root directory under which every product caches its downloads.

`ENV["GEMB_CACHE_PATH"]` if set, else the repository's own `data/` directory. Each product
appends its own subdirectory, so one setting moves all of them together:

    data/MCD43A3.061/         # MODIS granules and per-date samples
    data/satellite_albedo/    # C3S ordered timesteps
    data/invariant/<model>/   # ERA5-Land NetCDFs, Copernicus DEM tiles

    export GEMB_CACHE_PATH=/big/volume/gemb_cache   # to put them elsewhere

Defaulting beside the vendored tables and the derived products keeps everything for one
checkout in one place, which is what makes the driver scripts work with no environment set.
The cache subdirectories are gitignored; they are re-derivable bulk, unlike the products next
to them.

**These caches reach hundreds of gigabytes** — the MCD43A3 per-date samples alone are 184 GB —
so the volume holding the checkout has to be sized for it. Point `GEMB_CACHE_PATH` at bulk
storage if it is not. Avoid a home directory, which is commonly a small SSD with a quota;
nothing here writes to `\$HOME`, and the only paths derived from it are credential *reads*
(`~/.cdsapirc`, `~/.netrc`, `~/.edl_token`).

A package installed by `Pkg.add` has a read-only source tree, so the default is unusable there
and this throws rather than failing later inside a download. That is deliberate: the message
names the variable to set, whereas a silent fall back to `tempdir()` would put 184 GB somewhere
the OS reaps and turn a ~30 minute re-fold into a ~1 TB re-download.

Read at call time, not cached in a constant, so setting the variable inside a session takes
effect. Every function still takes an explicit `cache_path` that overrides this, and the
`data/` driver scripts keep their own `RGI7_*` variables for per-run overrides.
"""
function _gemb_cache_root()
    haskey(ENV, "GEMB_CACHE_PATH") && return ENV["GEMB_CACHE_PATH"]
    dir = normpath(joinpath(@__DIR__, "..", "data"))
    (isdir(dir) && _is_writable_dir(dir)) || throw(ErrorException(
        "the default cache location $(dir) is not a writable directory — which is normal for " *
        "a package installed by `Pkg.add`, whose source tree is read only. Set " *
        "ENV[\"GEMB_CACHE_PATH\"] to a writable volume with room for the product you are " *
        "using (the MCD43A3 sample cache alone reaches 184 GB), or pass `cache_path` " *
        "explicitly."))
    return dir
end

# Probe by creating and removing a file rather than reading the mode bits: the directory may sit
# on a read-only mount, or be owned by another user, and neither shows up as a permission bit
# the current process can usefully interpret.
function _is_writable_dir(dir::AbstractString)
    probe = joinpath(dir, ".gemb_write_probe")
    try
        touch(probe)
        rm(probe; force = true)
        return true
    catch
        return false
    end
end

"""
    _assert_bulk_cache(path, why)

Refuse a cache location under `tempdir()`.

[`_gemb_cache_root`](@ref) checks that the *default* location is writable; this checks that a
location the caller chose — via `GEMB_CACHE_PATH`, `RGI7_CACHE_PATH` or an explicit `cache_path`
— is somewhere bulk data can actually survive. The two are complementary, and only the
hundreds-of-gigabytes workflows need the second: a temp directory is reaped by the OS, and
`/tmp` is usually on the root filesystem. `why` names what the caller would lose, since that
differs between a run that downloads a terabyte and one that only re-folds what is already
cached.
"""
function _assert_bulk_cache(path::AbstractString, why::AbstractString)
    startswith(abspath(path), abspath(tempdir())) || return nothing
    error("""
        The cache resolves to $(path), which is under tempdir().
        $(why)
        Set a bulk-storage location on a large volume (NOT a home directory, which is
        commonly a small SSD with a quota):
            export GEMB_CACHE_PATH=/big/volume/gemb_cache   # all products
            export RGI7_CACHE_PATH=/big/volume/MCD43A3.061  # this cache only
        """)
end
