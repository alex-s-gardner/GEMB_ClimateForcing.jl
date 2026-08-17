"""
Tests for the CDS Retrieve API client (`src/cds_retrieve.jl`).

Offline only: these cover request-header construction, error-body parsing, and the
queue-limit classification that makes concurrent job ordering safe. Nothing here touches
the network — the live-endpoint paths are exercised by the opt-in integration sets in
`test_copernicus_albedo.jl`.
"""

using Test
using GEMB_ClimateForcing
import JSON3
import HTTP

using GEMB_ClimateForcing: _cds_headers, _cds_error_detail, CDSQueueLimitError,
    _CDS_QUEUE_LIMIT_HINT, _CDS_PENDING_STATUSES, _CDS_QUEUE_RETRY_BASE,
    _CDS_QUEUE_RETRY_MAX, _run_concurrent_jobs, CDSTransientError, _cds_is_transient,
    _cds_retry, _cds_check_response, _CDS_TRANSIENT_RETRY_BASE, _CDS_TRANSIENT_RETRY_MAX

@testset "CDS Retrieve client" begin
    @testset "Headers use PRIVATE-TOKEN, not Bearer" begin
        headers = _cds_headers("abc-123")
        names = first.(headers)
        @test "PRIVATE-TOKEN" in names
        # The Retrieve API rejects the Bearer scheme that the ARCO Zarr store uses; these
        # two auth paths share a key but not a header, and conflating them 403s.
        @test !("Authorization" in names)
        @test Dict(headers)["PRIVATE-TOKEN"] == "abc-123"
        @test Dict(headers)["Content-Type"] == "application/json"
    end

    @testset "Error detail extraction" begin
        body = JSON3.write(Dict("title" => "bad request",
                                "detail" => "invalid value for 'year'"))
        detail = _cds_error_detail(body)
        @test occursin("bad request", detail)
        @test occursin("invalid value for 'year'", detail)

        # Non-JSON and empty bodies degrade to something printable rather than throwing.
        @test occursin("<html>", _cds_error_detail("<html>502 Bad Gateway</html>"))
        @test _cds_error_detail("   ") == "(empty response body)"
    end

    @testset "Queue-limit classification" begin
        # The verbatim server message that must be recognised as retryable. Any rewording
        # upstream silently turns concurrent ordering back into hard failures, so this
        # asserts on the real string.
        real_message = "Number queued requests for this dataset is temporarily " *
                       "limited. Please configure your scripts accordingly "
        @test occursin(_CDS_QUEUE_LIMIT_HINT, real_message)

        # And a genuinely bad request must NOT be classified as retryable.
        @test !occursin(_CDS_QUEUE_LIMIT_HINT,
                        "invalid value for 'variable': 'albb_xx'")

        # "rejected" is deliberately not a pending status: it terminates the poll loop so
        # the retry decision is made once, rather than being polled until timeout.
        @test !("rejected" in _CDS_PENDING_STATUSES)
    end

    @testset "CDSQueueLimitError" begin
        err = CDSQueueLimitError("satellite-albedo", "job-1", "queue is full")
        @test err isa Exception
        msg = sprint(showerror, err)
        @test occursin("satellite-albedo", msg)
        @test occursin("job-1", msg)
        @test occursin("queue is full", msg)
        # The message has to point at the fix, since the request itself is fine.
        @test occursin("concurrent", msg)
    end

    @testset "Concurrent job scheduling" begin
        # The point of the whole exercise: N jobs that each spend their time waiting must
        # overlap, so wall-clock tracks the slowest job rather than their sum. Uses sleeps
        # in place of CDS calls, which is faithful because the real work is I/O waiting.
        function instrument(n_jobs, limit, stagger, work; gate=true)
            inflight = Threads.Atomic{Int}(0)
            peak = Threads.Atomic{Int}(0)
            submits = Float64[]
            done = Threads.Atomic{Int}(0)
            lk = ReentrantLock()
            t0 = time()
            _run_concurrent_jobs(n_jobs, limit, stagger; verbose=false) do i, submit_gate
                gate && submit_gate()
                lock(lk) do
                    push!(submits, time() - t0)
                end
                c = Threads.atomic_add!(inflight, 1) + 1
                # `peak` is only ever raised, and every writer holds a monotonic view.
                lock(lk) do
                    peak[] = max(peak[], c)
                end
                sleep(work)
                Threads.atomic_sub!(inflight, 1)
                Threads.atomic_add!(done, 1)
            end
            return (elapsed=time() - t0, peak=peak[], submits=sort(submits), done=done[])
        end

        # Six jobs of 1 s, all allowed in flight: overlapping, so well under 6 s.
        r = instrument(6, 6, 0.05, 1.0)
        @test r.done == 6
        @test r.elapsed < 3.0          # serial would be >= 6 s
        @test r.peak > 1               # genuinely concurrent

        # The stagger spaces submissions without serialising the waiting.
        gaps = diff(r.submits)
        @test all(g -> g >= 0.04, gaps)   # ~0.05 s apart, allowing scheduler slop

        # The first job does not sit out a stagger for nothing: with a 2 s stagger a lone
        # job still finishes promptly.
        r1 = instrument(1, 4, 2.0, 0.05)
        @test r1.done == 1
        @test r1.elapsed < 1.0

        # `max_concurrent` is honoured, and lowering it costs wall-clock.
        r2 = instrument(6, 2, 0.05, 1.0)
        @test r2.done == 6
        @test r2.peak <= 2
        @test r2.elapsed > r.elapsed

        # A job that never opens the gate (e.g. served from cache) is not delayed by it.
        r3 = instrument(4, 4, 5.0, 0.05; gate=false)
        @test r3.done == 4
        @test r3.elapsed < 1.0

        # Failures propagate rather than being swallowed (as a `CompositeException` from
        # `@sync`), and the semaphore is still released — otherwise one bad job would
        # deadlock the rest.
        @test_throws CompositeException _run_concurrent_jobs(3, 2, 0.0; verbose=false) do i, _
            i == 2 && error("job $(i) failed")
            sleep(0.01)
        end
        # Still usable afterwards: no leaked semaphore slot.
        r4 = instrument(3, 2, 0.0, 0.02)
        @test r4.done == 3

        # Degenerate counts are no-ops rather than errors.
        @test _run_concurrent_jobs((i, g) -> error("should not run"), 0, 4, 0.0;
                                   verbose=false) === nothing
    end

    @testset "Transient (5xx) classification" begin
        # A 502 from the fronting proxy says nothing about the request. This exact case
        # killed an 11 h global run: the bare `error()` it used to raise was
        # indistinguishable from a real failure, so two batches of folded work were lost.
        resp502 = HTTP.Response(502, "<html><head><title>502 Bad Gateway</title></head></html>")
        err = try
            _cds_check_response(resp502, "satellite-albedo", "costing")
            nothing
        catch e
            e
        end
        @test err isa CDSTransientError
        @test err.status == 502
        @test _cds_is_transient(err)
        msg = sprint(showerror, err)
        @test occursin("502", msg)
        @test occursin("transient", msg)

        for status in (500, 503, 504)
            e = try
                _cds_check_response(HTTP.Response(status, "upstream down"), "d", "submit")
            catch e
                e
            end
            @test e isa CDSTransientError && _cds_is_transient(e)
        end

        # 4xx must stay non-retryable: a bad request, a missing licence, or an over-limit
        # order fails the same way however many times it is sent. Retrying these is what
        # produced the endless resubmit loop documented at `satellite_albedo`.
        for status in (400, 401, 403, 404)
            e = try
                _cds_check_response(HTTP.Response(status, "nope"), "d", "submit")
            catch e
                e
            end
            @test !_cds_is_transient(e)
        end
        # Nor is a job the server reports as genuinely `failed` (an `ErrorException`).
        @test !_cds_is_transient(ErrorException("CDS job x failed"))
        @test !_cds_is_transient(ArgumentError("bad request"))

        # Connection-level faults are retryable — a dropped socket mid-run is not a verdict
        # on the request.
        @test _cds_is_transient(HTTP.Exceptions.ConnectError("https://x", ErrorException("refused")))
        @test _cds_is_transient(EOFError())
    end

    @testset "_cds_retry" begin
        # Succeeds after transient failures, without propagating them.
        n = 0
        got = _cds_retry(; dataset="d", context="test", verbose=false) do
            n += 1
            n < 3 && throw(CDSTransientError("d", "test", 502, "bad gateway"))
            :ok
        end
        @test got == :ok
        @test n == 3

        # A non-transient error escapes on the first attempt, unwrapped.
        m = 0
        @test_throws ArgumentError _cds_retry(; dataset="d", context="test", verbose=false) do
            m += 1
            throw(ArgumentError("bad request"))
        end
        @test m == 1

        # The attempt cap is honoured and the final failure is rethrown as itself, so the
        # user sees the service's own message rather than a wrapper.
        k = 0
        @test_throws CDSTransientError _cds_retry(; dataset="d", context="test",
                                                  attempts=3, verbose=false) do
            k += 1
            throw(CDSTransientError("d", "test", 500, "boom"))
        end
        @test k == 3

        # A deadline already in the past stops the retrying immediately, so a caller's own
        # timeout still bounds total time rather than being extended by retries.
        j = 0
        @test_throws CDSTransientError _cds_retry(; dataset="d", context="test",
                                                  verbose=false, deadline=time() - 1) do
            j += 1
            throw(CDSTransientError("d", "test", 503, "boom"))
        end
        @test j == 1
    end

    @testset "Transient backoff bounds" begin
        @test 0 < _CDS_TRANSIENT_RETRY_BASE <= _CDS_TRANSIENT_RETRY_MAX
        # Faster than the queue-limit backoff: a proxy hiccup clears in seconds, whereas a
        # queue slot needs another of our jobs to finish.
        @test _CDS_TRANSIENT_RETRY_BASE < _CDS_QUEUE_RETRY_BASE
        backoffs = [min(_CDS_TRANSIENT_RETRY_BASE * 2.0^(i - 1), _CDS_TRANSIENT_RETRY_MAX)
                    for i in 1:8]
        @test issorted(backoffs)
        @test all(<=(_CDS_TRANSIENT_RETRY_MAX), backoffs)
    end

    @testset "Retry backoff bounds" begin
        @test 0 < _CDS_QUEUE_RETRY_BASE <= _CDS_QUEUE_RETRY_MAX
        # A queue slot frees only when another of our jobs finishes (minutes), so polling
        # the submit endpoint every few seconds would just burn requests.
        @test _CDS_QUEUE_RETRY_BASE >= 10
        backoffs = [min(_CDS_QUEUE_RETRY_BASE * 2.0^(i - 1), _CDS_QUEUE_RETRY_MAX)
                    for i in 1:8]
        @test issorted(backoffs)
        @test all(<=(_CDS_QUEUE_RETRY_MAX), backoffs)
    end
end
