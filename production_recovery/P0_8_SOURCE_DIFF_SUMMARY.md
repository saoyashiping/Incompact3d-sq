# P0.8 Source Diff Summary

## `src/fibre_prod_state.f90`

Added read-only helpers for extracting first-fibre structure coordinates and structure velocity, returning zero velocity when no committed structure velocity is present. These helpers do not modify state.

## `src/fibre_prod_structure_dry_step.f90`

Added the P0.8 dry-step scratch predictor. It computes `dx_trial`, `x_trial`, and `u_trial` from `structure_input_force`, `dt`, and `rho_eff`, records diagnostics, checks boundedness, and keeps all trial data in scratch storage only.

## `src/fibre_prod_structure_dry_step_check.f90`

Added the P0.8 core check for dry-step allocation, predictor formula correctness, finite/bounded trial response, no production state commit, no sampled/hydro/input/velocity/RHS/force-buffer contamination, lambda=0 compatibility, lambda>0 without feedback force, and invalid shape guarding.

## `src/xcompact3d.f90`

Added default-off diagnostic-only dry-step preflight wiring guarded by `FIBRE_PROD_STRUCTURE_DRY_STEP_ENABLE`.

## `src/CMakeLists.txt`

Added `fibre_prod_structure_dry_step.f90` to `xcompact3d` and added the `fibre_prod_structure_dry_step_check` target.

## Production readiness boundary

P0.8 is a scratch dry-step preflight only. It does not commit trial state, does not generate reaction force, does not write Eulerian force buffers, does not modify DNS RHS, and does not mean production DNS-FSI is ready.
