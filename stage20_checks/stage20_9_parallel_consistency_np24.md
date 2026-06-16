# Stage 20.9: parallel consistency np=2/4 for two-way coupling

Stage 20.9 is a diagnostic-only helper-level decomposition consistency audit. It simulates `np=1`, `np=2`, and `np=4` ownership semantics for a controlled one-fibre closed-loop candidate path and verifies that global force, Eulerian integral, RHS-delta norm, structure statistics, action-reaction residual, and force-conservation residual remain invariant with respect to helper-rank decomposition.

This stage does not run actual MPI, `mpirun`, `mpiexec`, production DNS, production validation, production RHS updates, Stage 14 RHS injection, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, restart, statistics, or visualization production paths.

## Source-only and no-rerun policy

Stage 20.9 accepts Stage 20.8 PASS evidence when present and preserves source-only acceptance for prior Stage 20 closure behavior. It does not require old closure files, does not require all old stage outputs, and does not rerun Stage 10-19 or Stage 20.0-20.8.

## Helper decomposition

The helper simulates `np_values = 1,2,4` without MPI. The owner rule is documented as `modulo`:

```text
owner_rank(q) = global_point_id(q) mod np
local_point_id(q) = count of prior points with the same owner_rank
```

For each simulated `np`, the helper computes the same deterministic one-fibre candidate physics and then reduces rank-local diagnostics in deterministic Python arrays.

```text
u_relative_candidate = u_interp_candidate - V_current
F_fs_candidate = C_fs * u_relative_candidate
F_total_structure_candidate = F_b_candidate + F_T_candidate - F_fs_candidate
A_candidate = F_total_structure_candidate / rho_l
V_next_candidate = V_current + dt * A_candidate
X_next_candidate = X_current + dt * V_current + 0.5 * dt**2 * A_candidate
F_on_fluid_from_structure_candidate = +F_fs_candidate
f_eulerian_candidate = nearest-grid-point spread of F_on_fluid_from_structure_candidate * ds / dV
RHS_delta_zero = 0.0 * f_eulerian_candidate
RHS_delta_small = lambda_small * f_eulerian_candidate
```

The audit compares `np=2` and `np=4` against the `np=1` reference for global force, Eulerian integral force, RHS-delta-small norm, and structure global statistics.

## Safety boundary

All decomposition, X/V/A commits, reductions, force-density arrays, and RHS arrays are helper-local diagnostics only. Production two-way coupling, production structure-to-fluid coupling, production RHS coupling, production RHS updates, Stage 14 RHS injection, wall contact, fibre-fibre collision, and production multifibre logic remain disabled.

## Next stage

Stage 20.10: restart/statistics/visualization compatibility for active coupling.

## Manual command

```bash
stage20_checks/run_stage20_9_parallel_consistency_np24.sh
```

## Expected PASS evidence

```text
STAGE 20.9 PARALLEL CONSISTENCY NP24 VERDICT: PASS
STAGE 20.9 FINAL VERDICT: PASS
```
