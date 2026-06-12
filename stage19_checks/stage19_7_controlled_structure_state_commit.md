# Stage 19.7 controlled commit of production structure state

Stage 19.7 defines a diagnostic/helper-local controlled structure-state commit boundary.  It recreates helper-local Stage 19.3 initialization, Stage 19.4 force candidates, Stage 19.5 advance candidates, and Stage 19.6 closed-gate no-op markers, then commits only helper-local structure state when the Stage 19.7 commit gates are explicitly open.

The helper-local commit is:

```text
X_prod_after = X_next_candidate
V_prod_after = V_next_candidate
A_prod_after = A_advance_candidate
```

This is controlled helper-local structure-state commit only.  It is not production runtime hook insertion, not production runtime advance or commit activation, not fluid RHS coupling, not Stage 14 RHS injection, and not production restart/statistics/visualization I/O.

## Required arrays and markers

The controlled-commit schema includes:

* `X_prod_before`, `V_prod_before`, `A_prod_before`.
* `X_prod_after`, `V_prod_after`, `A_prod_after`.
* `F_b_candidate`, `F_T_candidate`, `F_h_candidate`, and `F_total_candidate`.
* `A_advance_candidate`, `V_next_candidate`, `X_next_candidate`, `delta_X_candidate`, and `delta_V_candidate`.
* `fluid_rhs_before`, `fluid_rhs_after`, `ibm_marker_before`, `ibm_marker_after`, `dns_core_marker_before`, `dns_core_marker_after`.
* `restart_io_marker_before`, `restart_io_marker_after`, `statistics_io_marker_before`, `statistics_io_marker_after`, `visualization_io_marker_before`, and `visualization_io_marker_after`.
* `owner_rank`, `global_point_id`, and `local_point_id`.

## Required metadata

The helper records geometry, physical coefficients, gate state, and ownership metadata:

* `n_fibre`, `n_point`, `component_dim`, `fibre_length`, `ds`, and `dt`.
* `rho_l`, `rho_tilde`, `bending_stiffness`, and `gamma`.
* `init_mode`, `sine_amplitude`, `sine_mode`, `tension_mode`, `tension_value`, and `controlled_force_amplitude`.
* `diagnostic_only`, `force_candidate_only`, `advance_candidate_only`, and `controlled_commit_enable`.
* `state_valid`, `container_initialized`, `physical_structure_enable`, `hook_enable`, `commit_allowed`, `rhs_spreading_allowed`, `stage14_rhs_injection_allowed`, `fluid_force_input_allowed`, `restart_io_allowed`, `statistics_io_allowed`, and `visualization_io_allowed`.

## Safe controlled-commit defaults

The default run uses a single helper-local fibre with `n_point = 64`, `component_dim = 3`, `fibre_length = 1.0`, `dt = 1.0e-4`, `rho_l = 1.0`, `rho_tilde = 1.0`, `bending_stiffness = 1.0e-3`, `gamma = 1.0e-3`, `init_mode = small_sine_fibre_zero_velocity`, `sine_amplitude = 1.0e-3`, `sine_mode = 1`, `tension_mode = constant`, `tension_value = 0.0`, and `controlled_force_amplitude = 0.0`.

The helper-local commit gates are intentionally open only for Stage 19.7 evidence: `controlled_commit_enable = true`, `physical_structure_enable = true`, and `commit_allowed = true`.  Runtime coupling gates remain closed: `hook_enable = false`, `rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed = false`, `fluid_force_input_allowed = false`, `restart_io_allowed = false`, `statistics_io_allowed = false`, and `visualization_io_allowed = false`.

## Validation policy

The helper validates that:

* The committed helper-local X/V/A arrays match the next candidate arrays exactly within tolerance.
* The force, acceleration, velocity, position, and delta candidate formulas remain consistent.
* The controlled commit occurs exactly once and remains helper-local.
* Fluid RHS, IBM, DNS-core, restart, statistics, and visualization markers are unchanged.
* Stage 10-18 and Stage 19.0-19.6 files remain unmodified.
* No production Fortran/CMake, production runtime hook, production runtime advance/commit, RHS coupling, Stage 14 RHS injection, MPI, DNS, contact/collision, or production multifibre logic is introduced.

Documentation strings, helper-local arrays, helper-local `physical_structure_enable=1`, and helper-local `commit_allowed=1` are not production runtime activation.

## Manual command

From the repository root, run:

```bash
bash stage19_checks/run_stage19_7_controlled_structure_state_commit.sh
```

The wrapper writes:

```text
stage19_outputs/fibre_stage19_7_controlled_structure_state_commit.dat
```

Expected PASS evidence:

```text
STAGE 19.7 CONTROLLED STRUCTURE STATE COMMIT VERDICT: PASS
STAGE 19.7 FINAL VERDICT: PASS
```
