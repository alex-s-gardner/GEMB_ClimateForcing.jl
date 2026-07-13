"""
AuthenticatedHTTPStore - HTTP store with Bearer token authentication.

Follows the pattern from Zarr.jl's GCStore implementation for authenticating
requests to cloud-optimized Zarr stores over HTTPS.
"""

using HTTP
using OpenSSL: OpenSSL
import Zarr

"""
    AuthenticatedHTTPStore <: Zarr.AbstractStore

HTTP store that supports Bearer token authentication via custom headers.

# Fields
- `url::String`: Base URL of the Zarr store
- `headers::Dict{String,String}`: HTTP headers (includes Authorization)
- `allowed_codes::Set{Int}`: HTTP status codes to treat as "key not found" (default: 404)

# Example
```julia
store = AuthenticatedHTTPStore(
    "https://example.com/path/to/data.zarr",
    token="your-bearer-token"
)
arr = Zarr.zopen(store)
```
"""
struct AuthenticatedHTTPStore <: Zarr.AbstractStore
    url::String
    headers::Dict{String,String}
    allowed_codes::Set{Int}
end

"""
    AuthenticatedHTTPStore(url::String; token::String, allowed_codes=Set((404,)))

Create an authenticated HTTP store for accessing Zarr data via HTTPS.

# Arguments
- `url::String`: Base URL of the Zarr store
- `token::String`: Bearer token for authentication
- `allowed_codes`: HTTP status codes to treat as "key not found"

# Returns
- `AuthenticatedHTTPStore` instance
"""
function AuthenticatedHTTPStore(url::String; token::String, allowed_codes::Set{Int}=Set((404,)))
    headers = Dict{String,String}(
        "Authorization" => "Bearer $(token)"
    )
    AuthenticatedHTTPStore(url, headers, allowed_codes)
end

"""
    Base.getindex(s::AuthenticatedHTTPStore, k::String)

Retrieve data from the authenticated HTTP store.

Follows the pattern from GCStore: passes headers to HTTP.request() for authentication.
"""
function Base.getindex(s::AuthenticatedHTTPStore, k::String)
    full_url = string(s.url, "/", k)

    r = HTTP.request(
        "GET",
        full_url,
        s.headers,
        status_exception=false,
        socket_type_tls=OpenSSL.SSLStream
    )

    if r.status >= 300
        if r.status in s.allowed_codes
            return nothing  # Key not found
        else
            error("HTTP $(r.status) accessing $(full_url): $(String(r.body))")
        end
    else
        return r.body
    end
end

# Implement required AbstractStore interface for read-only HTTP stores

Zarr.storagesize(s::AuthenticatedHTTPStore, p::String) = 0
Zarr.subdirs(s::AuthenticatedHTTPStore, p::String) = String[]
Zarr.subkeys(s::AuthenticatedHTTPStore, p::String) = String[]

function Zarr.isinitialized(s::AuthenticatedHTTPStore, i::String)
    s[i] !== nothing
end

function Base.setindex!(s::AuthenticatedHTTPStore, v, k::String)
    error("AuthenticatedHTTPStore is read-only")
end
