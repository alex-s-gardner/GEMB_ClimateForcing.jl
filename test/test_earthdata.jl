# Offline tests for the NASA Earthdata client (src/earthdata.jl).
#
# Nothing here touches the network: credential resolution runs under `withenv`, the
# transient/permanent split is asserted on constructed errors, and `_cmr_granules` runs
# against a trimmed real UMM-G response (test/fixtures/cmr_mcd43a3_granules.json) injected
# through the `requester` hook `EarthData.jl` provides for exactly this purpose.

using Test
using GEMB_ClimateForcing
using Dates
import JSON3
import HTTP
import Downloads

using GEMB_ClimateForcing: _cmr_granules, _granule_id, _granule_data_url, _granule_bytes,
    _size_unit_factor, _classifying_requester,
    _earthdata_is_transient, _earthdata_retry_after, _earthdata_check_response,
    _earthdata_body_excerpt, _netrc_credentials, EarthdataTransientError,
    _EARTHDATA_RETRY_ATTEMPTS, _EARTHDATA_RETRY_MAX, _EARTHDATA_RETRY_AFTER_MAX

const _CMR_FIXTURE = joinpath(@__DIR__, "fixtures", "cmr_mcd43a3_granules.json")

"""
    fixture_requester(; status=200, headers=[]) -> (requester, recorded)

A stand-in for `HTTP.request` that answers with the CMR fixture and records the calls it
received, so `_cmr_granules` can be driven offline. Mirrors the `recording_requester`
pattern in `EarthData.jl`'s own `test/search.jl`.

No `CMR-Search-After` header is returned, so upstream's pagination loop terminates after one
page — which is what a 3-item response means.
"""
function fixture_requester(; status::Int=200, headers=Pair{String,String}[],
                           body=read(_CMR_FIXTURE, String))
    recorded = []
    requester = function (method, url, hdrs; kwargs...)
        push!(recorded, (; method=String(method), url=String(url), headers=hdrs,
                         kwargs=Dict(kwargs)))
        return HTTP.Response(status, collect(headers); body=body)
    end
    return requester, recorded
end

@testset "NASA Earthdata client" begin

    @testset "Token resolution" begin
        # The environment variable wins, and is stripped.
        withenv("EARTHDATA_TOKEN" => "  tok-from-env  ") do
            @test get_earthdata_token() == "tok-from-env"
        end

        # An empty variable is treated as absent, not as an empty token — otherwise a
        # shell that exports it unset would produce a 401 instead of a setup message.
        mktempdir() do dir
            # Both variables, not just HOME: `homedir()` is libuv's `uv_os_homedir`, which
            # reads USERPROFILE on Windows and HOME elsewhere. Setting only HOME leaves the
            # Windows runner reading the real home directory, where the fallback file does
            # not exist.
            withenv("EARTHDATA_TOKEN" => "", "HOME" => dir, "USERPROFILE" => dir) do
                @test_throws ErrorException get_earthdata_token()
            end

            # Falls back to ~/.edl_token, skipping comments and blank lines.
            write(joinpath(dir, ".edl_token"), "# my token\n\ntok-from-file\n")
            withenv("EARTHDATA_TOKEN" => nothing, "HOME" => dir, "USERPROFILE" => dir) do
                @test get_earthdata_token() == "tok-from-file"
            end
        end

        # The error must name both the mint page and the two-token limit — the single most
        # common first-run failure is trying to create a third token.
        mktempdir() do dir
            withenv("EARTHDATA_TOKEN" => nothing, "HOME" => dir, "USERPROFILE" => dir) do
                msg = try
                    get_earthdata_token()
                    ""
                catch err
                    sprint(showerror, err)
                end
                @test occursin("urs.earthdata.nasa.gov/users/tokens", msg)
                @test occursin("two", msg)
                @test occursin("EARTHDATA_TOKEN", msg)
            end
        end
    end

    @testset ".netrc parsing" begin
        mktempdir() do dir
            path = joinpath(dir, "netrc")
            write(path, """
            machine example.com login wrong password wrong
            machine urs.earthdata.nasa.gov
                login someone
                password s3cret
            """)
            withenv("NETRC" => path) do
                @test _netrc_credentials("urs.earthdata.nasa.gov") == ("someone", "s3cret")
                # A machine with no credentials is an error, not a silent fall-through to
                # the wrong stanza.
                @test_throws ErrorException _netrc_credentials("nasa.example")
            end
            withenv("NETRC" => joinpath(dir, "absent")) do
                @test_throws ErrorException _netrc_credentials("urs.earthdata.nasa.gov")
            end
        end
    end

    @testset "Transient vs permanent" begin
        # Retryable: the service said nothing about the request.
        @test _earthdata_is_transient(EarthdataTransientError("ctx", 503, "down"))
        for status in (500, 502, 503, 504, 408, 429)
            r = HTTP.Response(status, ["Content-Type" => "text/plain"]; body="nope")
            @test_throws EarthdataTransientError _earthdata_check_response(r, "ctx")
        end
        @test _earthdata_is_transient(HTTP.Exceptions.ConnectError("u", ErrorException("x")))
        @test _earthdata_is_transient(EOFError())
        @test _earthdata_is_transient(Base.IOError("reset", -1))

        # NOT retryable. 403 in particular: an unaccepted EULA cannot be waited out, and
        # retrying it is how an endless request loop is born.
        for status in (400, 401, 403, 404, 422)
            r = HTTP.Response(status, String[]; body="nope")
            @test_throws ArgumentError _earthdata_check_response(r, "ctx")
        end
        @test !_earthdata_is_transient(ArgumentError("bad token"))
        @test !_earthdata_is_transient(ErrorException("truncated"))

        # A success returns nothing rather than throwing.
        @test isnothing(_earthdata_check_response(HTTP.Response(200, String[]; body="ok"),
                                                  "ctx"))

        # The 401 and 403 messages must say different things — they have different fixes.
        msg401 = try
            _earthdata_check_response(HTTP.Response(401, String[]; body="x"), "ctx")
        catch err
            sprint(showerror, err)
        end
        msg403 = try
            _earthdata_check_response(HTTP.Response(403, String[]; body="x"), "ctx")
        catch err
            sprint(showerror, err)
        end
        @test occursin("expire", msg401)
        @test occursin("licence", msg403) || occursin("application", msg403)
    end

    @testset "Retry-After and backoff bounds" begin
        mk(v) = HTTP.Response(429, ["Retry-After" => v]; body="slow down")
        @test _earthdata_retry_after(mk("30")) == 30.0
        @test _earthdata_retry_after(mk(" 12.5 ")) == 12.5
        @test _earthdata_retry_after(mk("")) == 0.0
        # The HTTP-date form is legal but unparsed here; it must degrade, not throw.
        @test _earthdata_retry_after(mk("Wed, 21 Oct 2015 07:28:00 GMT")) == 0.0
        # A pathological header cannot park a run for hours.
        @test _earthdata_retry_after(mk("99999")) == _EARTHDATA_RETRY_AFTER_MAX
        @test _earthdata_retry_after(HTTP.Response(429, String[]; body="x")) == 0.0

        # The attempt cap must not be the binding limit — the caller's deadline is. Same
        # assertion the CDS client carries, for the same reason.
        @test _EARTHDATA_RETRY_ATTEMPTS * _EARTHDATA_RETRY_MAX >= 3600
    end

    @testset "CMR granule search (offline, real UMM-G fixture)" begin
        requester, recorded = fixture_requester()
        granules = _cmr_granules(; short_name="MCD43A3", version="061",
                                 temporal="2019-07-01T00:00:00Z,2019-07-01T23:59:59Z",
                                 bounding_box="-51.0,66.0,-49.0,68.0",
                                 verbose=false, requester=requester)

        @test length(granules) == 3
        # One request: 3 items is short of the default page_size, so that is the last page.
        @test length(recorded) == 1
        # Paging is by `page_num`, NOT upstream's `all=true` search-after loop, which sends
        # `page_num` alongside `CMR-Search-After` and gets HTTP 400 `page_num is not allowed
        # with search-after` on page 2. Reproduced live; see the PR 5 note in earthdata.jl.
        @test occursin("page_num=1", recorded[1].kwargs[:body])
        @test !any(p -> lowercase(String(p[1])) == "cmr-search-after", recorded[1].headers)

        # `bounding_box` must reach CMR verbatim. This is the whole reason the shim calls
        # `EarthData.request` rather than `EarthData.granules`, whose keyword whitelist
        # rejects it — if a future refactor routes through `granules`, this fails.
        sent = recorded[1]
        body = sent.kwargs[:body]              # POST, form-encoded
        @test occursin("bounding_box", body)
        @test occursin(HTTP.URIs.escapeuri("-51.0,66.0,-49.0,68.0"), body)
        @test occursin("short_name=MCD43A3", body)
        # Unset keywords are dropped rather than sent as the string "nothing".
        @test !occursin("nothing", body)

        g1 = granules[1]
        @test g1.id == "MCD43A3.A2019175.h16v02.061.2020308035737"

        # `Type == "GET DATA"` picks the file out of EIGHT RelatedUrls. Note the `s3://`
        # copy of the same data is typed `GET DATA VIA DIRECT ACCESS`, and the DOI landing
        # page, the `.cmr.xml` sidecar, the browse imagery and the `s3credentials` endpoint
        # are all `VIEW RELATED INFORMATION` — an unfiltered `urls()` returns all of them.
        @test startswith(g1.url, "https://data.lpdaac.earthdatacloud.nasa.gov/")
        @test endswith(g1.url, ".hdf")
        @test !startswith(g1.url, "s3://")
        @test !occursin("s3credentials", g1.url)
        @test !occursin("doi.org", g1.url)
        @test !occursin(".cmr.xml", g1.url)
        @test !occursin("appeears", g1.url)   # AppEEARS is deliberately never used.

        # Size arrives as `Size` + `SizeUnit` ("MB"), so the unit is read, not assumed.
        @test g1.bytes == round(Int, 53.9037 * 1024^2)
        @test granules[2].bytes == round(Int, 72.2186 * 1024^2)

        # THE 16-DAY TRAP, from a real response: the query above asked for 2019-07-01 alone,
        # yet every granule's nominal date (the `A` field of the id) is in June — the
        # temporal extent spans the product's full 16-day retrieval window. Filtering must
        # use the granule id. Verified live: `CMR-Hits: 32` for this one-day, one-box query.
        @test JSON3.read(read(_CMR_FIXTURE, String)).hits == 32
        for g in granules
            @test occursin(r"\.A2019(175|176)\.", g.id)
        end
        # Two nominal dates and two tiles across three granules, so neither field alone
        # identifies a granule — `mcd43a3_granules` must filter on both.
        @test length(unique(id -> match(r"\.A(\d{7})\.", id)[1], [g.id for g in granules])) == 2
        @test length(unique(id -> match(r"\.(h\d\dv\d\d)\.", id)[1], [g.id for g in granules])) == 2

        # page_size is validated locally against CMR's ceiling rather than deferred to a
        # server error.
        @test_throws ArgumentError _cmr_granules(; short_name="X", version="1",
                                                 page_size=5000, requester=requester)
        @test_throws ArgumentError _cmr_granules(; short_name="X", version="1",
                                                 page_size=0, requester=requester)

        # A FULL page must fetch the next one, or a query over CMR's ceiling truncates
        # silently. With page_size=3 the 3-item fixture looks full, so the loop continues;
        # the second call returns an empty page and ends it.
        pages = 0
        function paging(args...; kwargs...)
            pages += 1
            body = pages == 1 ? read(_CMR_FIXTURE, String) :
                   """{"hits": 3, "took": 1, "items": []}"""
            return HTTP.Response(200, Pair{String,String}[]; body=body)
        end
        @test length(_cmr_granules(; short_name="MCD43A3", version="061", page_size=3,
                                   verbose=false, requester=paging)) == 3
        @test pages == 2
        # Page numbers advance rather than repeating page 1 forever.
        @test length(_cmr_granules(; short_name="MCD43A3", version="061", page_size=3,
                                   verbose=false,
                                   requester=(m, u, h; kwargs...) -> begin
                                       occursin("page_num=1", kwargs[:body]) ||
                                           error("expected page_num=1, got $(kwargs[:body])")
                                       HTTP.Response(200, Pair{String,String}[];
                                                     body="""{"hits":0,"took":1,"items":[]}""")
                                   end)) == 0
    end

    @testset "UMM-G record accessors" begin
        items = JSON3.read(read(_CMR_FIXTURE, String)).items
        # `DataGranule.Identifiers` carries the ProducerGranuleId; for MCD43A3 it equals
        # `GranuleUR`, but the typed identifier is read so the filter does not rely on that.
        @test _granule_id(items[1].umm) == String(items[1].umm.GranuleUR)

        # SizeUnit conversions, and the deliberate refusal to guess an unknown unit: the
        # figure feeds the truncated-download check, so a silent wrong number is worse than
        # an error.
        @test _size_unit_factor("Bytes") == 1
        @test _size_unit_factor("KB") == 1024
        @test _size_unit_factor("MB") == 1024^2
        @test _size_unit_factor("GB") == 1024^3
        @test _size_unit_factor(nothing) == 1024^2   # absent unit: MB, as the DAACs report
        @test_throws ArgumentError _size_unit_factor("furlongs")
    end

    @testset "Transient classification of CMR responses" begin
        # Upstream turns EVERY error status into a plain `ErrorException`, which is
        # (correctly) permanent here — so a 503 would abort the run instead of retrying.
        # `_classifying_requester` restores the split at the one point the status is visible.
        for status in (500, 502, 503, 504, 408, 429)
            requester, _ = fixture_requester(; status=status, body="upstream is down")
            wrapped = _classifying_requester(requester, "ctx")
            err = try
                wrapped("POST", "https://example.test", []; body="")
                nothing
            catch e
                e
            end
            @test err isa EarthdataTransientError
            @test _earthdata_is_transient(err)
            @test err.status == status
        end

        # A `Retry-After` on the 429 is carried through, not discarded.
        requester, _ = fixture_requester(; status=429, headers=["Retry-After" => "42"],
                                         body="slow down")
        err = try
            _classifying_requester(requester, "ctx")("POST", "https://x", []; body="")
        catch e
            e
        end
        @test err.retry_after == 42.0

        # Permanent statuses pass straight through, so upstream's `parse_cmr_error` message
        # (which reads CMR's `errors` array) is what the caller sees.
        for status in (400, 401, 403, 404)
            requester, _ = fixture_requester(; status=status, body="""{"errors":["nope"]}""")
            r = _classifying_requester(requester, "ctx")("POST", "https://x", []; body="")
            @test r.status == status
        end

        # And a permanent CMR error surfaces as a non-retryable exception from the shim.
        requester, _ = fixture_requester(; status=400, body="""{"errors":["bad param"]}""")
        err = try
            _cmr_granules(; short_name="X", version="1", verbose=false, requester=requester)
        catch e
            e
        end
        @test err isa Exception
        @test !_earthdata_is_transient(err)
        @test occursin("bad param", sprint(showerror, err))
    end

    @testset "Error body excerpting" begin
        @test _earthdata_body_excerpt(HTTP.Response(500, String[]; body="")) ==
              "(empty response body)"
        @test _earthdata_body_excerpt(HTTP.Response(500, String[]; body="  boom  ")) == "boom"
        long = HTTP.Response(500, String[]; body=repeat("x", 5000))
        @test occursin("truncated", _earthdata_body_excerpt(long; limit=100))
        @test length(_earthdata_body_excerpt(long; limit=100)) < 200
    end
end
