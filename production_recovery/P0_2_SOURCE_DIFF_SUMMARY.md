# P0.2 Source Diff Summary

## Files modified

- `src/fibre_prod_rhs_adapter.f90`
- `src/fibre_prod_main_hook_check.f90`
- `production_recovery/P0_2_PASS_FAIL.md`
- `production_recovery/P0_2_evidence/P0_2_VALIDATION_COMMAND.sh`
- `production_recovery/P0_2_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/P0_2_evidence/P0_2_FORCE_BUFFER_TO_RHS_AUDIT.txt`

## Fixes applied

1. The RHS adapter now applies the explicit production scaling
   `lambda_fsi * penalty_beta * force_buffer` when a valid force-density buffer is supplied.

2. The previous P0.2 failure at `ERROR STOP 12` was caused by a mismatch between the P0.2 check and the adapter implementation: the check expected scaled RHS increments, while the adapter added raw `force_x/force_y/force_z` values. This is now fixed.

3. `fibre_prod_main_hook_check.f90` now emits both:
   - `P0_1_FIBRE_PROD_MAIN_HOOK_CHECK PASS`
   - `P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS`

4. The P0.2 validation script now passes the user's real 2decomp-fft/Xcompact3D dependency path to CMake by default:
   `/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4`.

5. The P0.2 validation script no longer trusts a prewritten static PASS file. It rewrites `production_recovery/P0_2_PASS_FAIL.md` to PASS or FAIL based on the actual validation result.

## Preserved safety boundaries

- No pressure/projection/RK3/channel-forcing logic was modified.
- No long production DNS is run in P0.2.
- The old uniform RHS contribution path is not restored.
- `Production-run status` remains `STILL BLOCKED` after P0.2 PASS.
