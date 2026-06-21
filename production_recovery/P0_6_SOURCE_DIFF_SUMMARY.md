# P0.6 Source Diff Summary

## `src/fibre_prod_state.f90`

Added structure-side hydrodynamic input candidate storage and helper routines to allocate and attach `hydro_force_candidate(nnode,3)` without changing coordinates, sampled velocity, structure velocity, acceleration, or forces.

## `src/fibre_prod_hydro_input_candidate.f90`

Added the P0.6 hydrodynamic input candidate diagnostic module. It computes `relative_u = sampled_u - structure_u` and `candidate_force = beta_hydro * relative_u`, records diagnostics, and attaches the result to fibre-state diagnostic storage only.

## `src/fibre_prod_hydro_input_candidate_check.f90`

Added the P0.6 core check covering relative velocity correctness, candidate force correctness, attach-to-state behavior, no coordinate/sample/structure/velocity/RHS/force-buffer contamination, lambda=0 compatibility, lambda>0 without feedback force, and invalid shape guarding.

## `src/xcompact3d.f90`

Added default-off diagnostic-only hydro input candidate wiring guarded by `FIBRE_PROD_HYDRO_INPUT_CANDIDATE_ENABLE`.

## `src/CMakeLists.txt`

Added `fibre_prod_hydro_input_candidate.f90` to `xcompact3d` and added the `fibre_prod_hydro_input_candidate_check` target.

## Production readiness boundary

P0.6 is a structure-side diagnostic candidate stage only. It does not write Eulerian force buffers, does not modify DNS RHS, and does not mean production DNS-FSI is ready.
