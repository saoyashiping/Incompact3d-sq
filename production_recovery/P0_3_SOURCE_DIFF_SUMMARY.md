# P0.3 Source Diff Summary

## `src/fibre_prod_runtime_bridge.f90`

Added the runtime bridge container and APIs for RHS-shape initialization, force-buffer reset, finalization, and bridge-mediated no-op application through the production main hook.

## `src/fibre_prod_runtime_bridge_check.f90`

Added a standalone micro-check for shape initialization, force-buffer allocation, initialization no-contamination, lambda=0 no-contamination through bridge/main-hook, lambda>0 zero-buffer no physical response, and invalid-shape guard.

## `src/fibre_prod_main_hook.f90`

Added runtime-config query helpers while preserving P0.1/P0.2 APIs. The main hook still does not construct pseudo force and does not restore uniform RHS contribution.

## `src/xcompact3d.f90`

Added minimal runtime bridge call-path wiring near the existing fibre production RHS hook. The bridge initializes from RHS shape only. Lambda=0/disabled mode may pass through bridge + zero force buffer + main hook; lambda>0 without a physical force buffer remains on the blocked no-buffer path.

## `src/CMakeLists.txt`

Added `fibre_prod_runtime_bridge.f90` to the `xcompact3d` target and added the `fibre_prod_runtime_bridge_check` executable target without removing P0.1/P0.2 targets.

## Production readiness boundary

P0.3 verifies runtime bridge mechanics and lambda=0 no-contamination. It does not mean production DNS-FSI is ready.
