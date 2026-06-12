# Stage 19.6 production structure hook no-op invariance

Stage 19.6 defines a diagnostic/helper-level production-structure-hook boundary.  It verifies that when the production-side physical-structure gates are closed, the hook is a strict no-op:

```text
hook_after_state = hook_before_state
```

This is a no-op hook boundary only.  It does not insert a production runtime hook into `xcompact3d.f90` or any DNS path, does not advance or commit structure state, does not overwrite X/V/A, does not call Stage 14 RHS injection, does not spread force to Eulerian RHS, does not modify IBM/DNS-core/projection/Poisson/RK3, and does not modify production restart/statistics/visualization I/O.

## Helper-local model

The Stage 19.6 helper recreates helper-local evidence from the previous Stage 19 boundaries:

1. Stage 19.3-style initialized structure state.
2. Stage 19.4-style bending, tension, helper-fluid-placeholder, and total force candidates.
3. Stage 19.5-style acceleration, next-velocity, and next-position advance candidates.
4. Stage 19.6 no-op hook application with production activation gates closed.

No production Fortran or CMake integration is added in this stage.

## Required no-op invariance

The before and after helper-local hook state must match within tolerance for:

* `X`, `V`, and `A` structure arrays.
* Bending, tension, helper fluid-placeholder, and total force arrays.
* Advance candidates: acceleration, next velocity, and next position.
* Fluid RHS helper marker.
* IBM helper marker.
* DNS-core helper marker.
* Restart/statistics/visualization helper markers.

## Required arrays and markers

The schema includes these before/after fields:

* `X_before`, `V_before`, `A_before`, `F_b_before`, `F_T_before`, `F_h_before`, `F_total_before`.
* `A_advance_candidate_before`, `V_next_candidate_before`, `X_next_candidate_before`.
* `fluid_rhs_before`, `ibm_marker_before`, `dns_core_marker_before`, `restart_io_marker_before`, `statistics_io_marker_before`, `visualization_io_marker_before`.
* Matching `*_after` fields for every before field.
* `owner_rank`, `global_point_id`, and `local_point_id`.

## Required metadata

The schema records:

* `n_fibre`, `n_point`, `component_dim`, `fibre_length`, `ds`, and `dt`.
* `rho_l`, `rho_tilde`, `bending_stiffness`, and `gamma`.
* `init_mode`, `sine_amplitude`, `sine_mode`, `tension_mode`, `tension_value`, and `controlled_force_amplitude`.
* `diagnostic_only`, `force_candidate_only`, `advance_candidate_only`, and `hook_noop_only`.
* `state_valid`, `container_initialized`, `physical_structure_enable`, `hook_enable`, `commit_allowed`, `rhs_spreading_allowed`, `stage14_rhs_injection_allowed`, `fluid_force_input_allowed`, `restart_io_allowed`, `statistics_io_allowed`, and `visualization_io_allowed`.

## Safe defaults

The safe default run uses:

* `n_fibre = 1`, `n_point = 64`, `component_dim = 3`.
* `fibre_length = 1.0`, `dt = 1.0e-4`, and `ds = fibre_length / (n_point - 1)`.
* `rho_l = 1.0`, `rho_tilde = 1.0`, `bending_stiffness = 1.0e-3`, and `gamma = 1.0e-3`.
* `init_mode = small_sine_fibre_zero_velocity`, `sine_amplitude = 1.0e-3`, `sine_mode = 1`.
* `tension_mode = constant`, `tension_value = 0.0`, and `controlled_force_amplitude = 0.0`.
* `diagnostic_only = true`, `force_candidate_only = true`, `advance_candidate_only = true`, and `hook_noop_only = true`.
* `physical_structure_enable = false`, `hook_enable = false`, `commit_allowed = false`, `rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed = false`, `fluid_force_input_allowed = false`, `restart_io_allowed = false`, `statistics_io_allowed = false`, and `visualization_io_allowed = false`.

## Fail-closed policy

The helper fails if any required evidence is missing and cannot be safely accepted, if closed Stage 10-18 or Stage 19.0-19.5 files are modified, if production Fortran/CMake files are modified, if required fields are missing, if default safe values are wrong, if any no-op difference exceeds tolerance, or if production hook/advance/commit/RHS/I/O/MPI/DNS/contact/collision/multifibre activation is detected.

Documentation strings and helper-local marker names are not production activation.  Helper-local arrays and no-op hook checks are not production X/V/A state creation or production runtime state updates.

## Manual command

From the repository root, run:

```bash
bash stage19_checks/run_stage19_6_structure_hook_noop_invariance.sh
```

The wrapper writes:

```text
stage19_outputs/fibre_stage19_6_structure_hook_noop_invariance.dat
```

Expected PASS evidence:

```text
STAGE 19.6 STRUCTURE HOOK NOOP INVARIANCE VERDICT: PASS
STAGE 19.6 FINAL VERDICT: PASS
```
