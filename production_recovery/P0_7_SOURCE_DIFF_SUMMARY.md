# P0.7 Source Diff Summary

## `src/fibre_prod_state.f90`

Added controlled structure-side input storage and helper routines to allocate and attach `structure_input_force(nnode,3)` without changing coordinates, sampled velocity, hydro candidate storage, structure velocity, acceleration, or forces.

## `src/fibre_prod_structure_input_handoff.f90`

Added the P0.7 structure-input handoff module. It copies a finite `hydro_force_candidate(nnode,3)` into a structure-side input slot, records max/sum diagnostics, and attaches the result to fibre-state storage only.

## `src/fibre_prod_structure_input_handoff_check.f90`

Added the P0.7 core check for handoff allocation, handoff correctness, attach-to-state behavior, no coordinate/sample/hydro/structure/velocity/RHS/force-buffer contamination, lambda=0 compatibility, lambda>0 without feedback force, and invalid shape guarding.

## `src/xcompact3d.f90`

Added default-off diagnostic-only structure-input handoff call-path wiring guarded by `FIBRE_PROD_STRUCTURE_INPUT_HANDOFF_ENABLE`.

## `src/CMakeLists.txt`

Added `fibre_prod_structure_input_handoff.f90` to `xcompact3d` and added the `fibre_prod_structure_input_handoff_check` target.

## Production readiness boundary

P0.7 is a structure-side handoff preflight only. It does not advance structure, does not write Eulerian force buffers, does not modify DNS RHS, and does not mean production DNS-FSI is ready.
