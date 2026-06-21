# P0.5 Source Diff Summary

## `src/fibre_prod_state.f90`

Added minimal sampled fluid velocity storage and helper routines to allocate and attach sampled velocity without changing fibre coordinates, structure velocity, acceleration, or forces.

## `src/fibre_prod_velocity_bridge.f90`

Added `fibre_prod_velocity_bridge_sample_points` to sample arbitrary Lagrangian point arrays into caller-provided sampled velocity storage without RHS or force feedback.

## `src/fibre_prod_state_velocity_attachment.f90`

Added the state velocity attachment bridge from state node coordinates to velocity bridge sampling and then into state sampled velocity storage. The module does not call RHS adapter, does not call force-buffer-to-RHS hook, and does not advance structure.

## `src/fibre_prod_state_velocity_attachment_check.f90`

Added the P0.5 core check for state initialization, sampled velocity attachment, analytic velocity consistency, coordinate/structure/velocity/RHS no-contamination, lambda=0 compatibility, lambda>0 without feedback force, force-buffer no-write, and invalid shape guarding.

## `src/xcompact3d.f90`

Added default-off diagnostic-only state velocity attachment call-path wiring guarded by `FIBRE_PROD_STATE_VELOCITY_ATTACH_ENABLE`.

## `src/CMakeLists.txt`

Added `fibre_prod_state_velocity_attachment.f90` to `xcompact3d` and added the `fibre_prod_state_velocity_attachment_check` target.

## Production readiness boundary

P0.5 verifies one-way sampled velocity attachment to production fibre state. It does not mean production DNS-FSI is ready.
