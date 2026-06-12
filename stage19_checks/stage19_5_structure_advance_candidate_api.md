# Stage 19.5 structure advance candidate API

Stage 19.5 is an advance-candidate-boundary diagnostic only. It reconstructs a
helper-local Stage 19.3 initialized state, recomputes Stage 19.4 helper-local
force candidates, and computes candidate-only acceleration, next velocity, next
position, and delta arrays. It does not overwrite production `X_prod`, `V_prod`,
or `A_prod`; does not commit next-state candidates; does not insert production
hooks; does not spread force to Eulerian RHS; does not call Stage 14 RHS
injection; does not modify IBM/DNS-core/projection/Poisson/RK3; does not modify
production restart/statistics/visualization I/O; does not run MPI/DNS; and does
not introduce contact/collision/production multifibre logic.

## Candidate advance equations

The Stage 19.5 boundary validates the helper-local form
`rho_l * X_tt = F_total`:

* `A_advance_candidate = F_total_candidate / rho_l`
* `V_next_candidate = V_prod + dt * A_advance_candidate`
* `X_next_candidate = X_prod + dt * V_prod + 0.5 * dt**2 * A_advance_candidate`
* `delta_V_candidate = V_next_candidate - V_prod`
* `delta_X_candidate = X_next_candidate - X_prod`

With the default `rho_l = rho_tilde = 1.0`, the dimensional and nondimensional
acceleration consistency checks are equivalent.

## Required advance candidate fields

* `X_prod`
* `V_prod`
* `A_prod`
* `F_b_candidate`
* `F_T_candidate`
* `F_h_candidate`
* `F_total_candidate`
* `A_advance_candidate`
* `V_next_candidate`
* `X_next_candidate`
* `delta_X_candidate`
* `delta_V_candidate`
* `owner_rank`
* `global_point_id`
* `local_point_id`
* `n_fibre`
* `n_point`
* `component_dim`
* `fibre_length`
* `ds`
* `dt`
* `rho_l`
* `rho_tilde`
* `bending_stiffness`
* `gamma`
* `init_mode`
* `sine_amplitude`
* `sine_mode`
* `tension_mode`
* `tension_value`
* `controlled_force_amplitude`
* `diagnostic_only`
* `force_candidate_only`
* `advance_candidate_only`
* `state_valid`
* `container_initialized`
* `commit_allowed`
* `rhs_spreading_allowed`
* `stage14_rhs_injection_allowed`

## Safe defaults

* `n_fibre = 1`
* `n_point = 64`
* `component_dim = 3`
* `fibre_length = 1.0`
* `ds = fibre_length / (n_point - 1)`
* `dt = 1.0e-4`
* `rho_l = 1.0`
* `rho_tilde = 1.0`
* `bending_stiffness = 1.0e-3`
* `gamma = 1.0e-3`
* `init_mode = small_sine_fibre_zero_velocity`
* `sine_amplitude = 1.0e-3`
* `sine_mode = 1`
* `tension_mode = constant`
* `tension_value = 0.0`
* `controlled_force_amplitude = 0.0`
* `diagnostic_only = true`
* `force_candidate_only = true`
* `advance_candidate_only = true`
* `state_valid = true`
* `container_initialized = true`
* `commit_allowed = false`
* `rhs_spreading_allowed = false`
* `stage14_rhs_injection_allowed = false`

## Validation and fail-closed policy

Shape checks require all vector arrays to have shape `(N, 3)` and ownership
arrays to have shape `(N,)`. Numeric checks require `n_fibre = 1`, `n_point >=
8`, `component_dim = 3`, positive `fibre_length`, positive `dt`, positive
densities, nonnegative bending/gamma/tension/placeholder amplitudes, finite
arrays, deterministic owner rank, complete `global_point_id` coverage, and no
duplicate global IDs.

Formula checks validate `F_total_candidate = F_b_candidate + F_T_candidate +
F_h_candidate`, acceleration, next velocity, next position, and delta formulas.
The helper keeps all advance arrays helper-local, does not update production
`X_prod`, `V_prod`, or `A_prod`, and records no state commit or production
runtime state update.

If `diagnostic_only = true`, `commit_allowed` must be false. If
`rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed` must be false.
If `advance_candidate_only = true`, no production commit API may be called. If
`force_candidate_only = true`, production force arrays must not be modified.
Controlled helper-local `F_h_candidate` does not imply DNS fluid-force input.

Stage 19.5 preserves the Stage 19.0 source-only Stage 18 closure acceptance fix
and does not require rerunning Stage 18.0 through Stage 18.11.

## Manual command

```bash
stage19_checks/run_stage19_5_structure_advance_candidate_api.sh
```

Expected PASS evidence:

```text
STAGE 19.5 STRUCTURE ADVANCE CANDIDATE API VERDICT: PASS
STAGE 19.5 FINAL VERDICT: PASS
final_status PASS
```
