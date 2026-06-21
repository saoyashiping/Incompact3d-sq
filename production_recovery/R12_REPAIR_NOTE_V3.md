# Production Recovery R12 Repair Note V3

## Purpose

This repair package fixes the latest R12 failure by restoring the required R10/R11 runtime diagnostics source state and the fixed R11 validation post-processing script, then keeping the corrected R12 validation script.

## Diagnosis

The latest uploaded project regressed to an older R10/R11 state:

1. `src/fibre_prod_main_diagnostics.f90` did not contain `last_status` or `failed_calls`.
2. `src/fibre_prod_runtime_config.f90` did not read `FIBRE_PROD_DIAGNOSTICS_DIR`.
3. `src/fibre_prod_main_hook.f90` wrote only the old fixed R10 diagnostics filename and did not generate lambda-specific audit files for R11.
4. `src/xcompact3d.f90` used the older hook-ready logic.
5. `production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh` ran the six MPI cases but did not copy generated R10 audit/diagnostics into R11 audit files and did not recompute `R11_PASS_FAIL.md` or `R11_MPI_CONSISTENCY_AUDIT.md`.
6. Consequently, R12 found helper/no-fibre/restart evidence PASS but R11 rerun evidence FAIL/BLOCKED.

## Files restored or updated

```text
src/fibre_prod_runtime_config.f90
src/fibre_prod_main_diagnostics.f90
src/fibre_prod_main_hook.f90
src/xcompact3d.f90
production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh
production_recovery/R12_evidence/R12_VALIDATION_COMMAND_FIXED.sh
```

## Scope boundary

The source files above are not new R12 physics. They restore the already-required R10/R11 controlled-hook diagnostic and post-processing fixes that are prerequisites for R12. No RK3, pressure/projection, channel forcing, restart/statistics/visualization semantics, mesh, or new physical model is added.

## Expected result after rerun

R12 should be able to obtain:

```text
R12_STANDALONE_HELPER_AUDIT.md: Result: PASS
R12_R11_RERUN_AUDIT.md: Result: PASS
R12_NO_FIBRE_BASELINE_AUDIT.md: Result: PASS
R12_RESTART_STATS_VISU_AUDIT.md: Result: PASS
R12_FINAL_VALIDATION_MATRIX.md: Result: PASS
R12_FINAL_CLOSURE_REPORT.md: Result: PASS
R12_PASS_FAIL.md: Result: PASS
```
