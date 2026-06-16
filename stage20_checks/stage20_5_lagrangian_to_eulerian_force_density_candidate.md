# Stage 20.5 Lagrangian-to-Eulerian force-density candidate

Stage 20.5 is a helper-local, diagnostic-only force-density candidate audit. It
spreads the helper-local Lagrangian reaction-force candidate to a helper-local
Eulerian force-density candidate, but it does not add that candidate to
production RHS, call Stage 14 RHS injection, activate production IBM force
application, or activate production two-way coupling.

## Source-only and previous-stage policy

Stage 20.5 accepts Stage 20.4 PASS evidence when present and preserves Stage
20.0/20.1/20.2/20.3/20.4 source-only acceptance behavior. It does not rerun
Stage 10-19 or Stage 20.0-20.4, does not require old closure files, and does not
require all old stage output directories in source-only archives.

## Kernel choice

The helper uses a deterministic nearest-grid-point kernel. Each Lagrangian point
has one dimensionless kernel weight equal to one at its nearest helper-local grid
cell, so the weights sum to one per point. The force density is assembled as
`F_on_fluid_from_structure_candidate * weight * ds / dV`, making the Eulerian
integral match the Lagrangian total reaction force.

## Helper-local formulas

```text
u_relative_candidate = u_interp_candidate - V_current
F_fs_candidate = C_fs * u_relative_candidate
F_on_structure_from_fluid_candidate = -F_fs_candidate
F_on_fluid_from_structure_candidate = +F_fs_candidate
F_action_reaction_sum_candidate = F_on_structure_from_fluid_candidate + F_on_fluid_from_structure_candidate
lagrangian_total_reaction_force = sum_q F_on_fluid_from_structure_candidate(q) * ds
f_eulerian_candidate(i) = sum_q F_on_fluid_from_structure_candidate(q) * weight(q,i) * ds / dV
eulerian_integral_force_candidate = sum_i f_eulerian_candidate(i) * dV
force_conservation_residual_candidate = eulerian_integral_force_candidate - lagrangian_total_reaction_force
f_eulerian_effective = lambda_coupling * f_eulerian_candidate
```

Because lambda is zero and RHS coupling is disabled, `f_eulerian_effective` is
zero and no production RHS update occurs.

## Safety boundary

Stage 20.5 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 or Stage 20.0-20.4 file. It does not run MPI, DNS, previous
stages, builds, or production validation.

## Next stage

Stage 20.6: production RHS coupling activation with lambda gate.

## Manual command

```bash
stage20_checks/run_stage20_5_lagrangian_to_eulerian_force_density_candidate.sh
```

Expected output includes:

```text
STAGE 20.5 LAGRANGIAN TO EULERIAN FORCE DENSITY CANDIDATE VERDICT: PASS
STAGE 20.5 FINAL VERDICT: PASS
```
