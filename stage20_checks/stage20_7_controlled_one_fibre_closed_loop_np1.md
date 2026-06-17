# Stage 20.7: controlled one-fibre closed-loop response np=1

Stage 20.7 is a diagnostic-only, helper-local closed-loop audit for one fibre under `np=1` semantics. It connects the Stage 20.2 through Stage 20.6 helper-local concepts for several small steps: controlled fibre-point fluid velocity, fluid-to-structure force, hydrodynamic structure advance, helper-local controlled X/V/A commit, structure-to-fluid reaction force, Lagrangian-to-Eulerian force-density candidate, and lambda-gated RHS candidate response.

This stage does not run production DNS, MPI, production validation, production RHS writes, Stage 14 RHS injection, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, restart, statistics, or visualization production paths.

## Source-only and no-rerun policy

Stage 20.7 accepts Stage 20.6 PASS evidence when present and preserves source-only acceptance for prior Stage 20 closure behavior. It does not require old closure files, does not require all old stage outputs, and does not rerun Stage 10-19 or Stage 20.0-20.6.

## Helper-local loop formulas

At each step, the helper computes:

```text
u_relative_candidate = u_interp_candidate - V_current
F_fs_candidate = C_fs * u_relative_candidate
F_total_structure_candidate = F_b_candidate + F_T_candidate - F_fs_candidate
A_candidate = F_total_structure_candidate / rho_l
V_next_candidate = V_current + dt * A_candidate
X_next_candidate = X_current + dt * V_current + 0.5 * dt**2 * A_candidate
```

Then it performs a helper-local controlled commit only:

```text
X_current = X_next_candidate
V_current = V_next_candidate
A_current = A_candidate
```

The reaction and RHS candidate checks are also helper-local:

```text
F_on_structure_from_fluid_candidate = -F_fs_candidate
F_on_fluid_from_structure_candidate = +F_fs_candidate
F_action_reaction_sum_candidate = 0
f_eulerian_candidate = nearest-grid-point spread of F_on_fluid_from_structure_candidate * ds / dV
RHS_delta_zero = 0.0 * f_eulerian_candidate
RHS_delta_small = lambda_small * f_eulerian_candidate
```

`lambda_zero=0` must remain a strict RHS no-op at every step. `lambda_small` may produce only a bounded helper-local RHS response and is never written into a production RHS array.

## Safety boundary

`STAGE20_7_CONTROLLED_COMMIT_ENABLE=1` is permitted only for helper-local X/V/A state inside the audit. Production two-way coupling, production structure-to-fluid coupling, production RHS coupling, production RHS update, Stage 14 RHS injection, wall contact, fibre-fibre collision, and production multifibre logic remain disabled.

## Next stage

Stage 20.8: lambda=0 regression and small-lambda response comparison.

## Manual command

```bash
stage20_checks/run_stage20_7_controlled_one_fibre_closed_loop_np1.sh
```

## Expected PASS evidence

```text
STAGE 20.7 CONTROLLED ONE-FIBRE CLOSED-LOOP NP1 VERDICT: PASS
STAGE 20.7 FINAL VERDICT: PASS
```
