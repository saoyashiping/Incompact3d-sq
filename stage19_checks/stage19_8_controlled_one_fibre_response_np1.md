# Stage 19.8 controlled one-fibre production response np=1

Stage 19.8 defines a diagnostic/helper-local controlled one-fibre response under `np=1` semantics.  It recreates Stage 19.3 initialization, Stage 19.4 force candidates, Stage 19.5 advance candidates, and Stage 19.7 controlled helper-local commit for a small deterministic response history.

For each helper-local step:

```text
F_total_candidate = F_b_candidate + F_T_candidate + F_h_candidate
A_advance_candidate = F_total_candidate / rho_l
V_next_candidate = V_current + dt * A_advance_candidate
X_next_candidate = X_current + dt * V_current + 0.5 * dt**2 * A_advance_candidate
X_current, V_current, A_current <- X_next_candidate, V_next_candidate, A_advance_candidate
```

This is controlled helper-local one-fibre response only.  It does not insert production runtime hooks, does not activate production advance/commit paths, does not spread forces to Eulerian RHS, does not call Stage 14 RHS injection, does not modify IBM/DNS-core/projection/Poisson/RK3, does not modify production restart/statistics/visualization I/O, and does not run MPI or DNS.

## Required arrays, histories, and markers

The schema includes structure arrays (`X_initial`, `V_initial`, `A_initial`, `X_current`, `V_current`, `A_current`), candidate arrays (`X_next_candidate`, `V_next_candidate`, `A_advance_candidate`, `delta_X_candidate`, `delta_V_candidate`), force arrays (`F_b_candidate`, `F_T_candidate`, `F_h_candidate`, `F_total_candidate`), norm histories (`displacement_norm_history`, `velocity_norm_history`, `acceleration_norm_history`, `force_norm_history`, `bounded_response_history`), fluid/DNS/I/O markers, and ownership metadata (`owner_rank`, `global_point_id`, `local_point_id`).

## Safe defaults and bounds

The default response uses one fibre, `np_semantics = 1`, `n_point = 64`, `n_steps = 5`, `dt = 1.0e-5`, `fibre_length = 1.0`, `rho_l = rho_tilde = 1.0`, `bending_stiffness = gamma = 1.0e-5`, `sine_amplitude = 1.0e-4`, zero tension, and zero helper fluid-force placeholder.

The conservative helper-local bounds are documented and enforced as:

* `max_abs_displacement <= 1.0e-3`
* `max_abs_velocity <= 1.0`
* `max_abs_acceleration <= 1.0e3`
* `max_abs_force <= 1.0e3`

The helper-local response gates are intentionally open only for Stage 19.8 evidence: `controlled_response_enable = true`, `controlled_commit_enable = true`, `physical_structure_enable = true`, and `commit_allowed = true`.  Runtime coupling gates remain closed: `hook_enable = false`, `rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed = false`, `fluid_force_input_allowed = false`, `restart_io_allowed = false`, `statistics_io_allowed = false`, and `visualization_io_allowed = false`.

## Validation policy

The helper validates that the response is finite, bounded, deterministic, and single-fibre under `np=1` semantics without MPI.  It checks formulas at every step, exact helper-local commit count, unchanged fluid/IBM/DNS/I/O markers, preserved Stage 18/19 evidence, no closed-stage file modifications, no production Fortran/CMake changes, no runtime hook/advance/commit/RHS coupling, no Stage 14 RHS call, and no DNS/MPI/contact/collision/multifibre activation.

Documentation strings, helper-local arrays, helper-local `physical_structure_enable=1`, and helper-local `commit_allowed=1` are not production runtime activation.

## Manual command

From the repository root, run:

```bash
bash stage19_checks/run_stage19_8_controlled_one_fibre_response_np1.sh
```

The wrapper writes:

```text
stage19_outputs/fibre_stage19_8_controlled_one_fibre_response_np1.dat
```

Expected PASS evidence:

```text
STAGE 19.8 CONTROLLED ONE-FIBRE RESPONSE NP1 VERDICT: PASS
STAGE 19.8 FINAL VERDICT: PASS
```
