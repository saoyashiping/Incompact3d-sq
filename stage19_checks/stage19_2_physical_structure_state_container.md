# Stage 19.2 production physical structure state container

Stage 19.2 is a state-container-boundary diagnostic only. It defines the schema
for future production-side physical-structure arrays and ownership metadata, but
it does not insert production hooks, advance state, commit state, call Stage 14
RHS injection, spread force to Eulerian RHS, modify IBM/DNS-core/projection /
Poisson / RK3 channel forcing, modify production restart/statistics/visualization
I/O, run MPI, run DNS, or introduce contact/collision/production multifibre
logic. Local helper arrays are schema validation data, not production runtime
X/V/A state creation.

## Required container fields

* `n_fibre`
* `n_point`
* `component_dim`
* `ds`
* `fibre_length`
* `rho_l`
* `rho_tilde`
* `bending_stiffness`
* `gamma`
* `X_prod`
* `V_prod`
* `A_prod`
* `F_b_prod`
* `F_T_prod`
* `F_h_prod`
* `F_total_prod`
* `X_candidate`
* `V_candidate`
* `A_candidate`
* `owner_rank`
* `global_point_id`
* `local_point_id`
* `state_valid`
* `container_initialized`
* `diagnostic_only`
* `commit_allowed`
* `rhs_spreading_allowed`
* `stage14_rhs_injection_allowed`

## Safe defaults

* `n_fibre = 1`
* `n_point = 64`
* `component_dim = 3`
* `fibre_length = 1.0`
* `ds = fibre_length / (n_point - 1)`
* `rho_l = 1.0`
* `rho_tilde = 1.0`
* `bending_stiffness = 1.0e-3`
* `gamma = 1.0e-3`
* `diagnostic_only = true`
* `commit_allowed = false`
* `rhs_spreading_allowed = false`
* `stage14_rhs_injection_allowed = false`
* `state_valid = false` unless explicitly initialized by the helper
* `container_initialized = false` unless explicitly initialized by the helper

## Shape, numeric, and ownership rules

For a single fibre with `n_point = N` and `component_dim = 3`, `X_prod`,
`V_prod`, `A_prod`, `F_b_prod`, `F_T_prod`, `F_h_prod`, `F_total_prod`,
`X_candidate`, `V_candidate`, and `A_candidate` have shape `(N, 3)`. The
`owner_rank`, `global_point_id`, and `local_point_id` metadata arrays have shape
`(N,)`.

Numeric rules require `n_point >= 2`, `component_dim = 3`, `fibre_length > 0`,
`ds = fibre_length / (n_point - 1)`, `rho_l > 0`, `rho_tilde > 0`,
`bending_stiffness >= 0`, `gamma >= 0`, and finite helper arrays.

Ownership rules require `global_point_id` to cover `0..N-1` exactly once,
`local_point_id` to follow local helper ordering, and `owner_rank` to be
deterministic. Any np emulation is local metadata validation only; Stage 19.2
runs no MPI.

## Fail-closed rules

* If `diagnostic_only = true`, `commit_allowed` must be `false`.
* If `rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed` must be
  `false`.
* If `commit_allowed = false`, Stage 19.2 must not write back to production X/V/A
  state.
* If `state_valid = false`, no advance, commit, or hook is allowed.
* If the Stage 19.1 physical structure config remains default-disabled, the
  Stage 19.2 container remains inactive.
* Unknown or invalid container fields fail the helper.

## Evidence and false-positive policy

Stage 19.2 safely accepts preserved Stage 19.1 config-gate evidence, Stage 19.0
source-only Stage 18 closure acceptance, and Stage 18 closure evidence without
rerunning Stage 18.0 through Stage 18.11. Helper-local stage outputs are not
production I/O. Documentation strings, negative-check names, and state-container
schema text are not production runtime activation.

## Manual command

```bash
stage19_checks/run_stage19_2_physical_structure_state_container.sh
```

Expected PASS evidence:

```text
STAGE 19.2 PHYSICAL STRUCTURE STATE CONTAINER VERDICT: PASS
STAGE 19.2 FINAL VERDICT: PASS
final_status PASS
```
