# R10 Fix Note — finalized audit reliability repair

## Problem diagnosed

The latest R10 build and both `xcompact3d` np=1 runs completed, but `R10_PASS_FAIL.md` remained FAIL because the audit files were not regenerated. The observed audit files still contained stale `Result: BLOCKED` content.

The main technical weakness was that `src/xcompact3d.f90` used a single flag, `fibre_prod_r10_hook_ready`, for both applying the hook and deciding whether to finalize the hook. If one `fibre_prod_main_hook_apply` call returned a nonzero status, `fibre_prod_r10_hook_ready` was set to `.false.`, and the finalizer was skipped. Therefore no fresh audit file could be written, even though the main solver finished successfully.

## Fix applied

1. Added a separate `fibre_prod_r10_hook_initialized` flag in `src/xcompact3d.f90`.
2. R10 finalization now depends on initialization, not on whether further apply calls remain enabled.
3. If an apply call returns nonzero status, the hook stops applying further RHS changes, but finalization still writes diagnostics and audit files.
4. `fibre_prod_main_hook_apply` now records diagnostics even when the adapter returns a nonzero status.
5. R10 diagnostics now record `last_status` and `failed_calls`.
6. Audit PASS now requires `last_status=0` and `failed_calls=0`.
7. `R10_VALIDATION_COMMAND_FIXED.sh` now uses an absolute `FIBRE_PROD_DIAGNOSTICS_DIR`, removes stale audit files before running, and writes explicit FAIL files if the finalizer fails to generate audit output.
8. The PASS check remains strict: it only accepts exact `^Result: PASS` lines.

## Boundary

- `src/xcompact3d.f90` is modified only to make R10 finalization robust and auditable.
- The R10 hook site remains after `calculate_transeq_rhs(...)` and before `int_time(...)`.
- RK3 coefficients are not modified.
- Pressure/projection order is not modified.
- Channel forcing is not modified.
- Restart/statistics/visualization logic is not modified.
- R11/R12 are not entered.
