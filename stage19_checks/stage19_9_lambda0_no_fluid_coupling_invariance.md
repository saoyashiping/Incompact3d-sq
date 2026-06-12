# Stage 19.9 lambda=0 / no-fluid-coupling invariance

Stage 19.9 defines a diagnostic/helper-local lambda=0 no-fluid-coupling invariance boundary.  It recreates the Stage 19.8 controlled one-fibre response, computes a helper-local would-be coupling force, scales it by `lambda_coupling = 0`, and verifies that the effective coupling is exactly zero while every fluid/DNS/projection/Poisson/RK3/I/O marker remains unchanged.

For each helper-local step:

```text
F_total_candidate = F_b_candidate + F_T_candidate + F_h_candidate
A_advance_candidate = F_total_candidate / rho_l
V_next_candidate = V_current + dt * A_advance_candidate
X_next_candidate = X_current + dt * V_current + 0.5 * dt**2 * A_advance_candidate
F_coupling_effective = lambda_coupling * F_coupling_would_be
lambda_coupling = 0
```

This is lambda=0 / no-fluid-coupling invariance only.  It does not insert production runtime hooks, does not activate production advance/commit paths, does not spread force to Eulerian RHS, does not call Stage 14 RHS injection, does not modify IBM/DNS-core/projection/Poisson/RK3, does not modify production restart/statistics/visualization I/O, and does not run MPI or DNS.

## Required arrays, histories, and markers

The schema includes structure arrays, force and advance candidates, `F_coupling_would_be`, `F_coupling_effective`, displacement/velocity/acceleration/force/coupling norm histories, fluid RHS markers, IBM markers, DNS-core markers, projection/Poisson/RK3 markers, restart/statistics/visualization markers, and ownership metadata.

## Safe defaults and bounds

The default run uses one fibre, `np_semantics = 1`, `n_point = 64`, `n_steps = 5`, `dt = 1.0e-5`, `fibre_length = 1.0`, `lambda_coupling = 0.0`, `rho_l = rho_tilde = 1.0`, `bending_stiffness = gamma = 1.0e-5`, `sine_amplitude = 1.0e-4`, zero tension, and zero helper fluid-force placeholder.

Conservative helper-local bounds are:

* `max_abs_displacement <= 1.0e-3`
* `max_abs_velocity <= 1.0`
* `max_abs_acceleration <= 1.0e3`
* `max_abs_force <= 1.0e3`
* `max_abs_effective_coupling <= STAGE19_9_ZERO_TOL`

The helper-local structure response gates are open only for Stage 19.9 evidence.  Coupling gates remain closed: `hook_enable = false`, `rhs_spreading_allowed = false`, `stage14_rhs_injection_allowed = false`, `fluid_force_input_allowed = false`, `restart_io_allowed = false`, `statistics_io_allowed = false`, and `visualization_io_allowed = false`.

## Validation policy

The helper validates zero lambda, zero effective coupling at every step, zero effective-coupling norm history, bounded finite response, formula consistency, exact helper-local commit count, unchanged fluid/IBM/DNS/projection/Poisson/RK3/I/O markers, preserved Stage 18/19 evidence, no closed-stage file modifications, no production Fortran/CMake changes, no runtime hook/advance/commit/RHS coupling, no Stage 14 RHS call, and no DNS/MPI/contact/collision/multifibre activation.

Documentation strings, helper-local arrays, helper-local `physical_structure_enable=1`, helper-local `commit_allowed=1`, helper-local would-be coupling forces, and lambda-scaled zero effective coupling arrays are not production runtime activation.

## Manual command

From the repository root, run:

```bash
bash stage19_checks/run_stage19_9_lambda0_no_fluid_coupling_invariance.sh
```

The wrapper writes:

```text
stage19_outputs/fibre_stage19_9_lambda0_no_fluid_coupling_invariance.dat
```

Expected PASS evidence:

```text
STAGE 19.9 LAMBDA0 NO-FLUID-COUPLING INVARIANCE VERDICT: PASS
STAGE 19.9 FINAL VERDICT: PASS
```
