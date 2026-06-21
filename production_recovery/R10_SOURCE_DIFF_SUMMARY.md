# Production Recovery R10 Source Diff Summary

## Modified source files

1. `src/fibre_prod_runtime_config.f90`
   - Added `diagnostics_dir` to the runtime config.
   - Added support for environment variable `FIBRE_PROD_DIAGNOSTICS_DIR`.
   - Preserved default no-op behavior when no environment variables are set.

2. `src/fibre_prod_main_diagnostics.f90`
   - Added exact R10 audit writer.
   - Audit files now emit explicit `Result: PASS` or `Result: FAIL`.
   - Lambda=0 PASS requires hook calls, enabled=true, lambda=0, zero modified cells, finite signatures, and no RHS contamination.
   - Small-lambda PASS requires hook calls, enabled=true, lambda>0, modified cells >0, finite signatures, and nonzero controlled response.

3. `src/fibre_prod_main_hook.f90`
   - Finalizer now writes diagnostics and the corresponding R10 audit file into `FIBRE_PROD_DIAGNOSTICS_DIR`.
   - Lambda=0 run writes `R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt`.
   - Small-lambda run writes `R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt`.

## Modified documentation/evidence files

1. `production_recovery/R10_SOURCE_DIFF_SUMMARY.md`
2. `production_recovery/R10_PASS_FAIL.md`
3. `production_recovery/R10_evidence/R10_FIX_NOTE.md`
4. `production_recovery/R10_evidence/R10_VALIDATION_COMMAND_FIXED.sh`
5. Direct standalone gfortran evidence logs under `production_recovery/R10_evidence/`.

## `src/xcompact3d.f90` modified?

No additional change in this fix. The existing R10 hook site remains the controlled hook after `calculate_transeq_rhs(...)`, after the Stage 14 RHS injection candidate, and before `int_time(...)`.

## `src/CMakeLists.txt` modified?

No additional change in this fix.

## Connected to `xcompact3d` executable?

R10 remains connected through the existing controlled main-loop hook.

## Real RHS write?

Only the existing gated diagnostic RHS adapter is used. It writes only when `FIBRE_PROD_ENABLE=1` and `FIBRE_PROD_LAMBDA>0`. Lambda=0 remains a no-write path.

## RK3 modified?

No.

## Pressure/projection modified?

No.

## Channel forcing modified?

No.

## Restart/statistics/visualization modified?

No.

## R11/R12 entered?

No.

## Current result

BLOCKED pending re-run of R10 with the fixed diagnostics command. The previous run evidence shows that the hook check and `xcompact3d` runs can execute, but the audit files were stale and the shell-side PASS check was too permissive. Re-run R10 using `R10_VALIDATION_COMMAND_FIXED.sh` or the corrected command in the assistant response.
