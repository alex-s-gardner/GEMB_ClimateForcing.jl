# Latent finding (review 2026-10-13): a longitude in the 0–360 convention, or one near the ±180°
# seam, selected the wrong ERA5-Land cell. The ARCO stores hold longitude −179.9…180.0 ascending
# (with ~1e-11° round-off), and the nearest-index search compared raw differences.
#
# The axis is rebuilt here from that layout, so the script needs no network.
#
# Run (from the repository root):
#   julia --project=. scripts/audit/longitude_selection.jl                 # this branch
#   julia --project=<checkout of 352d7db> scripts/audit/longitude_selection.jl
#
# Expected: at 352d7db lon = 210 and lon = 359.97 both select 180.0 (not −150.0 and 0.0) and
# −179.97 selects −179.9 (not 180.0); on this branch all are correct.
include(joinpath(@__DIR__, "_common.jl"))

lon = collect(range(-179.9, 180.0; length = 3600)) .+ 1e-11
pick(x) = isdefined(G, :_nearest_longitude_index) ? G._nearest_longitude_index(lon, x) :
                                                    G.find_nearest_index(lon, x)
for (request, want) in ((210.0, -150.0), (-179.97, 180.0), (359.97, 0.0), (-38.46, -38.5))
    got = lon[pick(request)]
    println("lon = ", rpad(request, 8), " selects ", rpad(round(got; digits = 4), 9), " (nearest: ",
            want, ")", isapprox(got, want; atol = 1e-6) ? "" : "   WRONG")
end
