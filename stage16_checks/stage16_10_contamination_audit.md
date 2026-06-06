# Stage 16.10 RHS / IBM / structure contamination audit for first one-fibre FSI

Stage 16.10 is the final contamination audit before Stage 16 total closure. It does not add new physics. It reuses the closed Stage 16.7 small-lambda one-fibre diagnostic as runtime evidence and checks that the validated Stage 16 closed-loop path acts only through the approved chain:

```text
Stage 11 fluid-to-fibre sampling
-> Stage 12 feedback-force candidate
-> Stage 16.4 structure-side fluid-on-fibre force input
-> controlled structure-state diagnostics
-> Stage 13 force-density candidate
-> Stage 14 RHS diagnostic / controlled injection path
```

The audit fails closed if it detects direct RHS injection outside Stage 14, unapproved Stage 14 RHS calls, legacy/unauthorized IBM forcing, pressure/projection/Poisson/RK3/channel-forcing contamination, wall/contact or multi-fibre activation, unapproved bending/tension/full-structure solves, rank-unsafe diagnostics, NaN/Inf evidence, or Stage 14/15/16.1-16.9 regression evidence.

## User command

```bash
bash stage16_checks/run_stage16_10_contamination_audit.sh
```

The wrapper defaults to `BUILD_DIR=build_stage9`, configures the build directory with CMake if needed, builds only the already closed Stage 16.7 target `fibre_stage16_small_lambda_response_check`, runs that diagnostic under `STAGE16_10_NP=2`, copies its rank0-safe output into Stage 16.10 runtime evidence, and invokes the Stage 16.10 assertion helper.

## Default controls

Key defaults are:

- `STAGE16_10_RUN_STAGE16_9=0`
- `STAGE16_10_REQUIRE_STAGE14_CLOSED=1`
- `STAGE16_10_REQUIRE_STAGE15_CLOSED=1`
- `STAGE16_10_REQUIRE_STAGE16_9=1`
- `STAGE16_10_ACCEPT_STAGE16_9_CLOSED_EVIDENCE=1`
- `STAGE16_10_ENABLE=1`
- `STAGE16_10_ONE_FIBRE_FSI_ENABLE=1`
- `STAGE16_10_STRUCTURE_UPDATE_ENABLE=1`
- `STAGE16_10_TWO_WAY_RHS_ENABLE=1`
- `STAGE16_10_DIAGNOSTIC_ONLY=1`
- `STAGE16_10_NP=2`
- `STAGE16_10_NPTS=8`
- `STAGE16_10_FEEDBACK_ALPHA=1.0`
- `STAGE16_10_SMALL_LAMBDA=1.0e-8`
- `STAGE16_10_DT=1.0e-4`
- `STAGE16_10_RHO_TILDE=1.0`
- `STAGE16_10_MAX_FORCE_INPUT=1.0e-6`
- `STAGE16_10_MAX_STRUCTURE_UPDATE=1.0e-12`
- `STAGE16_10_MAX_RHS_INCREMENT=1.0e-8`
- `STAGE16_10_MAX_FLUID_DELTA=1.0e-8`

## Runtime diagnostics

The wrapper writes runtime evidence to:

```text
stage16_outputs/stage16_10_runtime_evidence.dat
```

The required Stage 16.10 summary is:

```text
stage16_outputs/fibre_stage16_10_contamination_audit.dat
```

The summary contains the required closed-loop status fields, approved-chain status, no-direct-RHS / no-unapproved-IBM / no-solver-contamination statuses, no wall/contact/multi-fibre statuses, no unapproved bending/tension/full-structure-solve statuses, rank0/rank-corruption statuses, no-NaN/Inf status, Stage 14/15/16.1-16.9 regression statuses, and `final_status`.

## PASS evidence

A passing run prints:

```text
STAGE 16.10 CONTAMINATION AUDIT VERDICT: PASS
STAGE 16.10 FINAL VERDICT: PASS
```

and the summary file reports:

- `closed_loop_path_status 1`
- `approved_stage11_12_16_4_13_14_chain_status 1`
- `approved_stage12_13_14_chain_status 1`
- `no_direct_rhs_injection_status 1`
- `no_unapproved_stage14_rhs_call_status 1`
- `no_legacy_ibm_forcing_status 1`
- `no_unapproved_production_ibm_forcing_status 1`
- all pressure/projection/Poisson/RK3/channel-forcing contamination statuses equal to `1`
- all wall/contact, multi-fibre, bending, tension, and full-structure solve contamination statuses equal to `1`
- `rank0_safe_diagnostic_status 1`
- `no_rank_corruption_status 1`
- `no_nan_inf_status 1`
- all Stage 14/15/16.1-16.9 regression statuses equal to `1`
- `final_status 1`

## False-positive-safe static audit policy

The Stage 16.10 helper intentionally reuses the corrected Stage 16.9 / Stage 16.8 / Stage 16.7 / Stage 16.6 / Stage 16.5 / Stage 16.4 helper pattern:

- documentation is checked for required-file existence, not scanned as executable regression evidence;
- negative-check strings are not treated as production behavior;
- regex literals such as `rg[[:space:]]` are not treated as real `rg` command usage;
- legitimate Stage 13.5 conservation/sign audit files are not classified as old production diagnostic regressions;
- old Stage 13.5 production force-density names are rejected only in real production/check logic;
- if evidence cannot be distinguished safely, the helper fails closed with explicit reasons rather than silently passing.

## Assumptions and risks

Stage 16.10 intentionally reuses closed Stage 16.7/16.8/16.9 evidence instead of adding a new Fortran physics module. The static audit is intentionally scoped to real source and executable check logic, not markdown or negative-check text, to avoid the known Stage 16.1-16.4 false-positive regressions.


## Stage 16.7 self-regression field handling

Stage 16.10 reuses the Stage 16.7 small-lambda diagnostic output as runtime evidence.
That Stage 16.7 output is produced by Stage 16.7 itself, so it intentionally does not
include a self-referential `stage16_7_regression_status` key. Stage 16.10 must infer
Stage 16.7 closed evidence from `final_status=1`, preserved `stage16_6_regression_status=1`,
and the presence of the Stage 16.7 wrapper/helper/doc/source/check files. Missing the
self-referential key is not a regression and must not produce
`missing_runtime_key_stage16_7_regression_status`, `runtime_stage16_7_regression_status_not_pass`,
or `summary_stage16_7_regression_status_not_pass`.
