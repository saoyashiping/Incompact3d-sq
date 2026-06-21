# P0.1 Source Diff Summary

## Modified files

### `src/CMakeLists.txt`

Added the complete production fibre module chain to the `xcompact3d` target while keeping existing standalone `fibre_prod_*_check` targets intact. This makes the main executable build-chain aware of the production modules but does not imply that the main program advances physical fibre state.

### `src/fibre_prod_rhs_adapter.f90`

Replaced the prior small-lambda uniform RHS constant injection with explicit optional Eulerian force-density buffers. Lambda=0 is a strict no-op. Lambda>0 without force buffers returns diagnostic blocking status `13` and leaves RHS unchanged. Lambda>0 with force buffers adds the provided force-density buffers to the matching RHS components and records diagnostic increment quantities.

Uniform RHS injection has been blocked: the old `contribution = lambda * penalty_beta * dt` path is not present, and the adapter no longer adds one global scalar to all `rhs_x` cells.

### `src/fibre_prod_main_hook.f90`

Kept the original compatible main-program call form `call fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status)` and added optional `force_x`, `force_y`, and `force_z` buffers. The hook does not construct pseudo force and does not apply uniform RHS increments; it delegates all RHS loading to the adapter.

### `src/fibre_prod_main_hook_check.f90`

Updated the hook check to cover lambda=0 no-op, lambda>0 missing-buffer blocking, lambda>0 nonzero force-buffer RHS loading, and lambda>0 zero force-buffer no-op behavior. The expected success token is `P0_1_FIBRE_PROD_MAIN_HOOK_CHECK PASS`.

## Production readiness boundary

P0.1 PASS means build-chain integration and RHS safety gating are implemented. It does not mean production DNS-FSI is ready.
