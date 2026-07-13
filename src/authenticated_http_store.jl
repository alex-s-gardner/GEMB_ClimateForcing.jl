"""
AuthenticatedHTTPStore - HTTP store with Bearer token authentication and caching.

Follows the pattern from Zarr.jl's GCStore implementation for authenticating
requests to cloud-optimized Zarr stores over HTTPS. Includes LRU caching to
minimize repeated HTTP requests for the same chunks.
"""

using HTTP
using OpenSSL: OpenSSL
using LRUCache
import Zarr

"""
    AuthenticatedHTTPStore <: Zarr.AbstractStore

HTTP store that supports Bearer token authentication via custom headers with LRU caching.

# Fields
- `url::String`: Base URL of the Zarr store
- `headers::Dict{String,String}`: HTTP headers (includes Authorization)
- `allowed_codes::Set{Int}`: HTTP status codes to treat as "key not found" (default: 404)
- `cache::LRU{String,Union{Vector{UInt8},Nothing}}`: LRU cache for HTTP responses

# Example
```julia
store = AuthenticatedHTTPStore(
    "https://example.com/path/to/data.zarr",
    token="your-bearer-token",
    cache_size=128  # Cache up to 128 chunks
)
arr = Zarr.zopen(store)
```
"""
struct AuthenticatedHTTPStore <: Zarr.AbstractStore
    url::String
    headers::Dict{String,String}
    allowed_codes::Set{Int}
    cache::LRU{String,Union{Vector{UInt8},Nothing}}
end

"""
    AuthenticatedHTTPStore(url::String; token::String, allowed_codes=Set((404,)), cache_size=128)

Create an authenticated HTTP store for accessing Zarr data via HTTPS with caching.

# Arguments
- `url::String`: Base URL of the Zarr store
- `token::String`: Bearer token for authentication
- `allowed_codes`: HTTP status codes to treat as "key not found"
- `cache_size::Int`: Maximum number of chunks to cache (default: 128)

# Returns
- `AuthenticatedHTTPStore` instance with LRU cache

# Notes
- Cache size of 128 chunks ≈ 1-10 MB depending on chunk size
- Cache is shared across all variables in the same store
- Metadata (.zarray, .zattrs, .zgroup) is also cached
"""
function AuthenticatedHTTPStore(url::String; token::String, allowed_codes::Set{Int}=Set((404,)), cache_size::Int=128)
    # Build headers with Bearer token (following GCStore pattern)
    headers = Dict{String,String}(
        "Authorization" => "Bearer $(token)"
    )
    # Create LRU cache for HTTP responses
    cache = LRU{String,Union{Vector{UInt8},Nothing}}(maxsize=cache_size)
    AuthenticatedHTTPStore(url, headers, allowed_codes, cache)
end

"""
    Base.getindex(s::AuthenticatedHTTPStore, k::String)

Retrieve data from the authenticated HTTP store with caching.

Follows the pattern from GCStore: passes headers to HTTP.request() for authentication.
Uses LRU cache to avoid repeated HTTP requests for the same keys.
"""
function Base.getindex(s::AuthenticatedHTTPStore, k::String)
    # Check cache first
    if haskey(s.cache, k)
        return s.cache[k]
    end

    full_url = string(s.url, "/", k)

    # Make authenticated request (following GCStore pattern)
    r = HTTP.request(
        "GET",
        full_url,
        s.headers,  # Pass headers with Bearer token
        status_exception=false,
        socket_type_tls=OpenSSL.SSLStream
    )

    # Handle response
    result = if r.status >= 300
        if r.status in s.allowed_codes
            nothing  # Key not found
        else
            error("HTTP $(r.status) accessing $(full_url): $(String(r.body))")
        end
    else
        r.body
    end

    # Cache the result (including nothing for 404s)
    s.cache[k] = result
    return result
end

# Implement required AbstractStore interface for read-only HTTP stores

"""
    Zarr.storagesize(s::AuthenticatedHTTPStore, p::String)

Return storage size. For HTTP stores, this is typically unknown (return 0).
"""
Zarr.storagesize(s::AuthenticatedHTTPStore, p::String) = 0

"""
    Zarr.subdirs(s::AuthenticatedHTTPStore, p::String)

List subdirectories. HTTP stores typically don't support directory listing.
"""
Zarr.subdirs(s::AuthenticatedHTTPStore, p::String) = String[]

"""
    Zarr.subkeys(s::AuthenticatedHTTPStore, p::String)

List subkeys. HTTP stores typically don't support key listing.
"""
Zarr.subkeys(s::AuthenticatedHTTPStore, p::String) = String[]

"""
    Zarr.isinitialized(s::AuthenticatedHTTPStore, i::String)

Check if a key exists in the store.
"""
function Zarr.isinitialized(s::AuthenticatedHTTPStore, i::String)
    s[i] !== nothing
end

"""
    Base.setindex!(s::AuthenticatedHTTPStore, v, k::String)

HTTP stores are read-only - this operation is not supported.
"""
function Base.setindex!(s::AuthenticatedHTTPStore, v, k::String)
    error("AuthenticatedHTTPStore is read-only")
end
