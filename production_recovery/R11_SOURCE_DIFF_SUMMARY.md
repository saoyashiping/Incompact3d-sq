# Production Recovery R11 Source Diff Summary

## New files

- `production_recovery/R11_PLAN.md`
- `production_recovery/R11_BUILD_LOG.txt`
- `production_recovery/R11_PASS_FAIL.md`
- `production_recovery/R11_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R11_RUN_LOG_hook_check.txt`
- `production_recovery/R11_RUN_LOG_lambda0_np1.txt`
- `production_recovery/R11_RUN_LOG_lambda0_np2.txt`
- `production_recovery/R11_RUN_LOG_lambda0_np4.txt`
- `production_recovery/R11_RUN_LOG_smalllambda_np1.txt`
- `production_recovery/R11_RUN_LOG_smalllambda_np2.txt`
- `production_recovery/R11_RUN_LOG_smalllambda_np4.txt`
- `production_recovery/R11_evidence/README.md`
- `production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh`
- `production_recovery/R11_evidence/R11_MPI_CONSISTENCY_AUDIT.md`
- R11 lambda=0 audit files for np=1/2/4.
- R11 small-lambda audit files for np=1/2/4.
- R11 hook-diagnostics placeholder files for lambda=0 and small-lambda np=1/2/4.

## Modified files

No source files were modified in R11.

## `src/xcompact3d.f90` modified?

No.  R11 did not modify `src/xcompact3d.f90`.

## R10 hook source modified?

No.  R11 did not modify `src/fibre_prod_main_hook.f90`, `src/fibre_prod_rhs_adapter.f90`, `src/fibre_prod_main_diagnostics.f90`, or `src/fibre_prod_runtime_config.f90`.

## RK3 modified?

No.  R11 did not modify RK3.

## Pressure/projection modified?

No.  R11 did not modify pressure/projection.

## Channel forcing modified?

No.  R11 did not modify channel forcing.

## Restart/statistics/visualization modified?

No.  R11 did not modify restart/statistics/visualization semantics.

## Entered R12?

No.  R11 did not enter R12.

## Final result

R11 is BLOCKED in this environment because the validation script could not configure with `mpif90`; no hook-check executable or xcompact3d executable was produced, and no np=1/2/4 MPI run evidence exists.
