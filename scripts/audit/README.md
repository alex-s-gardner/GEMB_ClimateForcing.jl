# Audit scripts

Reproductions of the GEMB_ClimateForcing findings in the GEMB pipeline review of 2026-10-13
(GlobalGlacierAnalysis.jl `notes/gemb_pipeline_review_2026-10-13.md`). Each script builds its
inputs synthetically, needs no network or data, and runs unchanged against the reviewed commit
352d7db and against `fix/forcing-audit`, printing which API it found. Its header states the
finding and the output expected from each version.

| script | finding |
|---|---|
| `double_lapse_damping.jl` | 2: on-glacier lapse rate applied, then decoupled again |
| `truncated_excess_fit.jl` | 11: k and lapse rate fitted to an excess truncated at melting |
| `humidity_identity.jl` | 10: dewpoint and RH on different curves; ΔT = 0 and Δz = 0 not identities |
| `decoupling_longwave.jl` | 18: decoupling lowered downward longwave |
| `longitude_selection.jl` | latent: 0–360 and near-seam longitudes select the wrong cell |
| `nan_validation.jl` | latent: NaN passes `validate_climate_forcing_units` |

`_common.jl` holds the shared stack builder and the cross-cell sums, carrying the field names of
both accumulator layouts.

## Running

From the repository root, on this branch:

    julia --project=. scripts/audit/double_lapse_damping.jl

Against the reviewed commit, from a checkout of it (the scripts are read from this branch):

    git worktree add ../GEMB_ClimateForcing-352d7db 352d7db
    julia --project=../GEMB_ClimateForcing-352d7db -e 'using Pkg; Pkg.instantiate()'
    julia --project=../GEMB_ClimateForcing-352d7db scripts/audit/double_lapse_damping.jl
