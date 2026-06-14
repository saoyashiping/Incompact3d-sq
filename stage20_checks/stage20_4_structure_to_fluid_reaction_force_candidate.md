# Stage 20.4 structure-to-fluid reaction force candidate

Stage 20.4 is a helper-local, diagnostic-only Lagrangian reaction-force
candidate audit. It derives the fluid-side candidate from the Stage 20 sign
convention and verifies Newton third-law consistency, but it does not spread to
an Eulerian grid, add to fluid RHS, activate production structure-to-fluid
coupling, or activate production two-way coupling.

## Source-only and previous-stage policy

Stage 20.4 accepts Stage 20.3 PASS evidence when present and preserves Stage
20.0/20.1/20.2/20.3 source-only acceptance behavior. It does not rerun Stage
10-19 or Stage 20.0-20.3, does not require old closure files, and does not
require all old stage output directories in source-only archives.

## Helper-local formulas

```text
u_relative_candidate = u_interp_candidate - V_current
F_fs_candidate = C_fs * u_relative_candidate
F_on_structure_from_fluid_candidate = -F_fs_candidate
F_on_fluid_from_structure_candidate = +F_fs_candidate
F_action_reaction_sum_candidate = F_on_structure_from_fluid_candidate + F_on_fluid_from_structure_candidate
F_total_structure_candidate_with_fluid = F_b_candidate + F_T_candidate - F_fs_candidate
reaction_force_norm_candidate = norm(F_on_fluid_from_structure_candidate)
action_reaction_residual_norm_candidate = norm(F_action_reaction_sum_candidate)
```

## Gate distinction

`STAGE20_4_STRUCTURE_TO_FLUID_CANDIDATE_ENABLE=1` is allowed only for
helper-local Lagrangian reaction-force construction. `STAGE20_4_STRUCTURE_TO_FLUID_ENABLE`
remains false, production reaction force remains disabled, RHS spreading remains
disabled, and lambda remains zero. The helper-local reaction candidate is not an
Eulerian RHS force density.

## Safety boundary

Stage 20.4 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 or Stage 20.0-20.3 file. It does not run MPI, DNS, previous
stages, builds, or production validation.

## Next stage

Stage 20.5: Lagrangian-to-Eulerian force-density coupling candidate.

## Manual command

```bash
stage20_checks/run_stage20_4_structure_to_fluid_reaction_force_candidate.sh
```

Expected output includes:

```text
STAGE 20.4 STRUCTURE TO FLUID REACTION FORCE CANDIDATE VERDICT: PASS
STAGE 20.4 FINAL VERDICT: PASS
```
