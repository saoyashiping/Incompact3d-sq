# Stage 19.3 production physical structure initialization

Stage 19.3 is an initialization-boundary diagnostic only. It validates
helper-local initialization of the Stage 19.2 state-container schema, but it does
not create production runtime X/V/A state, insert production hooks, advance or
commit structure state, call Stage 14 RHS injection, spread force to Eulerian
RHS, modify IBM/DNS-core/projection/Poisson/RK3 channel forcing, modify
production restart/statistics/visualization I/O, run MPI, run DNS, or introduce
contact/collision/production multifibre logic.

## Required initialization modes

* `straight_fibre_zero_velocity`
* `small_sine_fibre_zero_velocity`
* `straight_fibre_controlled_velocity`
* `small_sine_fibre_controlled_velocity`
* `controlled_helper_force_placeholder`

## Required initialized arrays

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

## Required metadata and safe defaults

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
* `velocity_amplitude = 0.0`
* `controlled_force_amplitude = 0.0`
* `diagnostic_only = true`
* `state_valid = true` after helper-local initialization
* `container_initialized = true` after helper-local initialization
* `commit_allowed = false`
* `rhs_spreading_allowed = false`
* `stage14_rhs_injection_allowed = false`

## Initialization definitions

Straight-fibre initialization sets `X_prod(i,:) = [i * ds, 0, 0]`, zeroes
`A_prod`, `F_b_prod`, and `F_T_prod`, zeroes `F_h_prod` unless the helper-local
controlled force placeholder is enabled, sets `F_total_prod = F_b_prod +
F_T_prod + F_h_prod`, and copies `X_candidate = X_prod`, `V_candidate = V_prod`,
and `A_candidate = A_prod`.

Small-sine initialization sets `X_prod(i,0) = i * ds`, `X_prod(i,1) =
sine_amplitude * sin(2*pi*sine_mode*i*ds/fibre_length)`, and `X_prod(i,2) = 0`.
It uses the same zero acceleration, zero force placeholders, candidate copy, and
fail-closed behavior as straight-fibre initialization.

Controlled velocity sets `V_prod(i,1) = velocity_amplitude *
sin(2*pi*sine_mode*i*ds/fibre_length)` when `velocity_amplitude > 0`; otherwise
`V_prod(:,:) = 0`.

The controlled helper-local force placeholder sets `F_h_prod(i,1) =
controlled_force_amplitude * sin(2*pi*sine_mode*i*ds/fibre_length)` when
`controlled_force_amplitude > 0`; otherwise `F_h_prod(:,:) = 0`. This is
helper-local initialization evidence only: it does not activate DNS fluid-force
input, does not spread force to Eulerian RHS, and does not call Stage 14 RHS
injection.

## Shape, numeric, ownership, and consistency rules

For `N = n_point` and `component_dim = 3`, `X_prod`, `V_prod`, `A_prod`,
`F_b_prod`, `F_T_prod`, `F_h_prod`, `F_total_prod`, `X_candidate`,
`V_candidate`, and `A_candidate` have shape `(N, 3)`. The `owner_rank`,
`global_point_id`, and `local_point_id` arrays have shape `(N,)`.

Numeric checks require `n_fibre = 1`, `n_point >= 2`, `component_dim = 3`,
`fibre_length > 0`, `rho_l > 0`, `rho_tilde > 0`, `bending_stiffness >= 0`,
`gamma >= 0`, finite small `sine_amplitude`, positive integer `sine_mode`,
`velocity_amplitude >= 0`, `controlled_force_amplitude >= 0`, finite arrays,
`global_point_id` covering `0..N-1` exactly once, no duplicate global IDs,
deterministic `owner_rank`, `F_total_prod = F_b_prod + F_T_prod + F_h_prod`, and
candidate arrays equal to production arrays at initialization.

## Fail-closed policy

If `diagnostic_only = true`, `commit_allowed` must be false. If
`rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed` must be false.
If `commit_allowed = false`, Stage 19.3 must not write back to existing
production X/V/A state. Controlled helper-local `F_h_prod` does not imply DNS
fluid-force input activation. Stage 19.3 remains inactive with respect to
production runtime when Stage 19.1 physical-structure configuration is
default-disabled.

Stage 19.3 preserves the Stage 19.0 source-only Stage 18 closure acceptance fix
and does not require rerunning Stage 18.0 through Stage 18.11. Helper-local stage
outputs are not production I/O, documentation strings are not runtime activation,
and local helper initialization is not a production runtime state update.

## Manual command

```bash
stage19_checks/run_stage19_3_physical_structure_initialization.sh
```

Expected PASS evidence:

```text
STAGE 19.3 PHYSICAL STRUCTURE INITIALIZATION VERDICT: PASS
STAGE 19.3 FINAL VERDICT: PASS
final_status PASS
```
