# Production Recovery R10 Source Diff Summary

## Modified source files

1. `src/fibre_prod_main_diagnostics.f90`
   - Added `last_status` and `failed_calls` to the R10 diagnostics state.
   - `fibre_prod_main_diagnostics_record(...)` now accepts an optional `status_code`.
   - Lambda=0 PASS now requires zero failed calls and `last_status=0`.
   - Small-lambda PASS now requires zero failed calls and `last_status=0`.
   - Audit files now report `last_status` and `failed_calls`.

2. `src/fibre_prod_main_hook.f90`
   - `fibre_prod_main_hook_apply(...)` now records diagnostics even if the RHS adapter returns a nonzero status.
   - This prevents silent loss of R10 failure evidence.

3. `src/xcompact3d.f90`
   - Added a separate `fibre_prod_r10_hook_initialized` flag.
   - R10 finalization is now called when the hook was initialized, even if later apply calls were disabled after a nonzero status.
   - Added a short rank-zero message when R10 apply is disabled by a nonzero status.
   - The hook site remains the same: after `calculate_transeq_rhs(...)`, after the Stage 14 RHS candidate, and before `int_time(...)`.

## Modified validation/evidence files

1. `production_recovery/R10_evidence/R10_VALIDATION_COMMAND_FIXED.sh`
   - Uses absolute `FIBRE_PROD_DIAGNOSTICS_DIR`.
   - Removes stale audit files before new validation.
   - Creates explicit FAIL audit files if the finalizer does not generate them.
   - Checks only exact `^Result: PASS` lines.

2. `production_recovery/R10_evidence/R10_FIX_NOTE.md`
3. `production_recovery/R10_SOURCE_DIFF_SUMMARY.md`

## `src/CMakeLists.txt` modified?

No.

## Connected to `xcompact3d` executable?

Yes. R10 remains the first controlled main-loop hook stage.

## Real RHS write?

Only through the existing gated R10 RHS adapter. No write occurs for `FIBRE_PROD_LAMBDA=0`. A small diagnostic RHS contribution is applied only when `FIBRE_PROD_ENABLE=1` and `FIBRE_PROD_LAMBDA>0`.

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

BLOCKED pending real Ubuntu re-run of:

```bash
bash production_recovery/R10_evidence/R10_VALIDATION_COMMAND_FIXED.sh
```

After the re-run, `R10_PASS_FAIL.md` must be written by the script as PASS only if:

1. `R10_FIBRE_PROD_MAIN_HOOK_CHECK PASS` is printed;
2. `R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt` contains `Result: PASS`;
3. `R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt` contains `Result: PASS`.
