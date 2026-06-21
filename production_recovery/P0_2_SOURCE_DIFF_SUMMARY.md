# P0.2 Source Diff Summary

## Modified files

### `src/fibre_prod_main_hook.f90`

Added the public `fibre_prod_main_hook_apply_force_buffer` API. It validates hook initialization, runtime config, RHS shape, and `fibre_prod_force_buffer_type` allocation/finiteness/shape before delegating to the RHS adapter with `force_buffer%fx`, `force_buffer%fy`, and `force_buffer%fz`.

### `src/fibre_prod_main_hook_check.f90`

Preserved the P0.1 checks and added a buffer-level API smoke check that constructs a small `fibre_prod_force_buffer_type`, calls `fibre_prod_main_hook_apply_force_buffer`, checks scaled RHS increments, and emits `P0_2_FIBRE_PROD_MAIN_HOOK_BUFFER_API_CHECK PASS`.

### `src/fibre_prod_rhs_adapter.f90`

Kept the P0.1 safety gate and updated the explicit force-buffer path to apply `lambda_fsi * penalty_beta` scaling to nonzero force-density buffers. The old uniform RHS contribution path was not restored.

### `src/fibre_prod_force_buffer_rhs_path_check.f90`

Added the P0.2 core micro-check. It initializes a small grid, allocates and resets a production force buffer, generates a buffer through IBM spreading, verifies force conservation, routes the buffer through the main hook into RHS, checks scaled RHS equality and volume-integrated RHS conservation, checks lambda=0 no-op, and checks invalid-buffer blocking.

### `src/CMakeLists.txt`

Added the `fibre_prod_force_buffer_rhs_path_check` target and updated `fibre_prod_main_hook_check` dependencies so both hook-level and buffer-path checks can compile with the buffer-level main-hook API.

## Production readiness boundary

P0.2 PASS means the IBM-generated Eulerian force-density buffer can be passed through the production main hook into DNS RHS in a controlled micro-check. It does not mean production DNS-FSI is ready.
