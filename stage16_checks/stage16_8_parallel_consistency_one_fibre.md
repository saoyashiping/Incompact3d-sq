# Stage 16.8 np=1/2/4 parallel consistency for first one-fibre FSI

Stage 16.8 reuses the closed Stage 16.7 small-lambda diagnostic path and launches it with MPI process counts 1, 2, and 4. The check treats np=1 as the reference and verifies that np=2 and np=4 preserve the same controlled one-fibre response within strict parallel tolerances.

The validated chain under comparison is:

```text
Stage 11 fluid-to-fibre sampling
-> Stage 12 feedback-force candidate
-> Stage 16.4 structure-side fluid-on-fibre force input
-> controlled structure-state update diagnostics
-> Stage 13 force-density candidate
-> Stage 14 RHS diagnostic / controlled injection path with small nonzero lambda
```

Stage 16.8 does not add production hooks, does not connect to `xcompact3d.f90`, and does not modify pressure/projection/Poisson/RK3/channel-forcing numerics. It does not activate wall/contact, multi-fibre, or legacy IBM production forcing.

## User command

```bash
bash stage16_checks/run_stage16_8_parallel_consistency_one_fibre.sh
```

The wrapper defaults to `BUILD_DIR=build_stage9`, configures the build directory with CMake if needed, builds only the already closed Stage 16.7 target `fibre_stage16_small_lambda_response_check`, runs it for `STAGE16_8_NP_LIST="1 2 4"`, copies each rank0-written diagnostic into a distinct Stage 16.8 per-np file, and invokes the Stage 16.8 assertion helper.

## Default controls

Key defaults are:

- `STAGE16_8_NP_LIST="1 2 4"`
- `STAGE16_8_SMALL_LAMBDA=1.0e-8`
- `STAGE16_8_ZERO_LAMBDA=0.0`
- `STAGE16_8_MAX_RHS_INCREMENT=1.0e-8`
- `STAGE16_8_MAX_FLUID_DELTA=1.0e-8`
- `STAGE16_8_MIN_RHS_INCREMENT=1.0e-20`
- `STAGE16_8_MAX_PARALLEL_FORCE_DIFF=1.0e-14`
- `STAGE16_8_MAX_PARALLEL_STRUCTURE_DIFF=1.0e-14`
- `STAGE16_8_MAX_PARALLEL_RHS_DIFF=1.0e-14`
- `STAGE16_8_MAX_PARALLEL_FLUID_DIFF=1.0e-14`

## Runtime diagnostics

Per-np diagnostics are written as:

```text
stage16_outputs/stage16_8_np1_small_lambda_response.dat
stage16_outputs/stage16_8_np2_small_lambda_response.dat
stage16_outputs/stage16_8_np4_small_lambda_response.dat
```

The summary is written as:

```text
stage16_outputs/fibre_stage16_8_parallel_consistency_one_fibre.dat
```

The summary includes np run/final statuses, np=1/2/4 RHS and fluid-signature values, force/structure/RHS/fluid differences versus np=1, rank-corruption status, approved-chain status, solver-contamination guard statuses, Stage 14/15/16.1-16.7 regression statuses, and `final_status`.

## PASS evidence

A passing run prints:

```text
STAGE 16.8 PARALLEL CONSISTENCY ONE-FIBRE VERDICT: PASS
STAGE 16.8 FINAL VERDICT: PASS
```

and the summary file reports:

- `np1_run_status 1`, `np2_run_status 1`, and `np4_run_status 1`
- `np1_final_status 1`, `np2_final_status 1`, and `np4_final_status 1`
- `parallel_force_status 1`
- `parallel_structure_status 1`
- `parallel_rhs_status 1`
- `parallel_fluid_status 1`
- `rank0_safe_diagnostic_status 1`
- `no_rank_corruption_status 1`
- all approved-chain and no-contamination statuses equal to `1`
- `final_status 1`

## False-positive-safe static audit policy

The Stage 16.8 helper intentionally reuses the corrected Stage 16.7 / Stage 16.6 / Stage 16.5 / Stage 16.4 helper pattern:

- documentation is checked for required-file existence, not scanned as executable regression evidence;
- negative-check strings are not treated as production behavior;
- regex literals such as `rg[[:space:]]` are not treated as real `rg` command usage;
- legitimate Stage 13.5 conservation/sign audit files are not classified as old production diagnostic regressions;
- old Stage 13.5 production force-density names are rejected only in real production/check logic;
- if evidence cannot be distinguished safely, the helper fails closed with explicit reasons rather than silently passing.

## Assumptions and risks

Stage 16.8 intentionally reuses the Stage 16.7 target instead of duplicating physics. The Stage 16.7 diagnostic itself remains a logical one-fibre diagnostic; Stage 16.8 records the MPI launch count in each copied per-np file as the Stage 16.8 `np` value while preserving the original Stage 16.7 reported value as `stage16_7_reported_np`.
