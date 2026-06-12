# Stage 19.4 bending/tension force candidate API

Stage 19.4 is a force-candidate-boundary diagnostic only. It reconstructs a
helper-local Stage 19.3 initialized state and computes candidate-only bending,
tension, helper-local fluid-placeholder, and total forces. It does not create
production runtime X/V/A state, insert production hooks, advance or commit state,
spread force to Eulerian RHS, call Stage 14 RHS injection, modify IBM/DNS-core /
projection / Poisson / RK3 channel forcing, modify production restart/statistics /
visualization I/O, run MPI, run DNS, or introduce contact/collision/production
multifibre logic.

## Candidate equations

The physical structure equation boundary is
`rho_l * X_tt = d_s(T X_s) - B * X_ssss + F_h`.

Stage 19.4 validates the candidate-only definitions:

* `F_b_candidate = -B * X_ssss`, or nondimensionally `F_b_candidate = -gamma * X_ssss`.
* `F_T_candidate = d_s(T X_s)`.
* `F_h_candidate` is a helper-local placeholder only.
* `F_total_candidate = F_b_candidate + F_T_candidate + F_h_candidate`.

## Required force-candidate fields

* `X_prod`
* `V_prod`
* `A_prod`
* `F_b_candidate`
* `F_T_candidate`
* `F_h_candidate`
* `F_total_candidate`
* `X_candidate`
* `V_candidate`
* `A_candidate`
* `owner_rank`
* `global_point_id`
* `local_point_id`
* `n_fibre`
* `n_point`
* `component_dim`
* `fibre_length`
* `ds`
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
* `state_valid = true`
* `container_initialized = true`
* `commit_allowed = false`
* `rhs_spreading_allowed = false`
* `stage14_rhs_injection_allowed = false`

## Candidate validation rules

The helper uses the Stage 19.3 small-sine initialized state with `X_candidate =
X_prod`, `V_candidate = V_prod`, and `A_candidate = A_prod`. A deterministic
finite-difference fourth-derivative stencil validates `F_b_candidate = -gamma *
X_ssss`. For `tension_mode = constant`, a deterministic second-derivative form
validates `F_T_candidate = d_s(T X_s)`; the safe default `tension_value = 0`
therefore produces zero tension force. The helper-local placeholder validates
`F_h_candidate = 0` by default, or a sine-shaped helper-local placeholder when
`controlled_force_amplitude > 0`. This never activates DNS fluid-force input,
RHS spreading, or Stage 14 RHS injection.

Shape checks require each vector candidate array to have shape `(N, 3)` and the
ownership arrays to have shape `(N,)`. Numeric checks require `n_fibre = 1`,
`n_point >= 8`, `component_dim = 3`, positive length/densities, nonnegative
bending/gamma/tension/placeholder amplitudes, finite arrays, deterministic owner
rank, complete `global_point_id` coverage, no duplicate global IDs, total-force
consistency, unchanged candidate state arrays, no state advance, and no commit.

## Fail-closed policy

If `diagnostic_only = true`, `commit_allowed` must be false. If
`rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed` must be false.
If `force_candidate_only = true`, no structure advance or commit API may be
called. Stage 19.4 remains inactive with respect to production runtime when Stage
19.1 physical-structure configuration is default-disabled. Force-candidate schema
text and local helper arrays are not production runtime force activation.

Stage 19.4 preserves the Stage 19.0 source-only Stage 18 closure acceptance fix
and does not require rerunning Stage 18.0 through Stage 18.11.

## Manual command

```bash
stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh
```

Expected PASS evidence:

```text
STAGE 19.4 BENDING TENSION FORCE CANDIDATE API VERDICT: PASS
STAGE 19.4 FINAL VERDICT: PASS
final_status PASS
```
