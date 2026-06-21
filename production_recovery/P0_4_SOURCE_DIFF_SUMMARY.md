# P0.4 Source Diff Summary

## `src/fibre_prod_velocity_bridge.f90`

Added one-way velocity sampling APIs and a runtime diagnostic helper. The module allocates sample points, stores interpolated velocity samples, calls the production IBM interpolation kernel, and does not call RHS adapter or force-buffer-to-RHS APIs.

## `src/fibre_prod_velocity_bridge_check.f90`

Added a core P0.4 check for analytic interpolation accuracy, velocity no-contamination, RHS no-contamination, lambda=0 runtime bridge compatibility, lambda>0 without feedback force, and invalid point guarding.

## `src/xcompact3d.f90`

Added default-off one-way velocity sampling call-path wiring guarded by `FIBRE_PROD_VELOCITY_SAMPLE_ENABLE`. The path reads `ux1/uy1/uz1` only and does not write RHS.

## `src/CMakeLists.txt`

Added `fibre_prod_velocity_bridge.f90` to `xcompact3d` and added the `fibre_prod_velocity_bridge_check` executable target while preserving P0.2/P0.3 targets.

## Production readiness boundary

P0.4 verifies one-way velocity sampling mechanics. It does not mean production DNS-FSI is ready.

## P0.4 repair addendum

The original P0.4 attempt failed before dynamic validation because two static integration requirements were incomplete:

- `xcompact3d.f90` called the velocity sampling bridge but did not import `fibre_prod_velocity_bridge` or declare/initialize the velocity sampling guard variables.
- `src/CMakeLists.txt` did not include `fibre_prod_velocity_bridge.f90` in the `xcompact3d` target and did not define the `fibre_prod_velocity_bridge_check` target.

This repair adds only the missing one-way diagnostic integration pieces. It does not modify pressure/projection/RK3/channel-forcing, does not generate feedback force, and does not make the code production DNS-FSI ready.

