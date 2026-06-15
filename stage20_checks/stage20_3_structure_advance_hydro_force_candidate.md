# Stage 20.3 structure advance with hydrodynamic force candidate

Stage 20.3 is a helper-local, diagnostic-only hydrodynamic-force-included
structure advance candidate. It computes candidate acceleration, velocity, and
position updates from `F_b_candidate`, `F_T_candidate`, and `F_fs_candidate`, but
it does not production-commit those candidates or activate production structure
advance.

## Source-only and previous-stage policy

Stage 20.3 accepts Stage 20.2 PASS evidence when present and preserves Stage
20.0/20.1/20.2 source-only acceptance behavior. It does not rerun Stage 10-19 or
Stage 20.0-20.2, does not require old closure files, and does not require all
old stage output directories in source-only archives.

## Helper-local formulas

```text
u_relative_candidate = u_interp_candidate - V_current
F_fs_candidate = C_fs * u_relative_candidate
F_total_structure_candidate_without_fluid = F_b_candidate + F_T_candidate
F_total_structure_candidate_with_fluid = F_b_candidate + F_T_candidate - F_fs_candidate
A_hydro_candidate = F_total_structure_candidate_with_fluid / rho_l
V_next_hydro_candidate = V_current + dt * A_hydro_candidate
X_next_hydro_candidate = X_current + dt * V_current + 0.5 * dt**2 * A_hydro_candidate
delta_V_hydro_candidate = V_next_hydro_candidate - V_current
delta_X_hydro_candidate = X_next_hydro_candidate - X_current
```

## Gate distinction

`STAGE20_3_FLUID_TO_STRUCTURE_ENABLE=1` and
`STAGE20_3_STRUCTURE_ADVANCE_HYDRO_CANDIDATE_ENABLE=1` are allowed only for
helper-local candidate force and X/V/A computation. They do not imply production
runtime structure advance, production commit, two-way coupling,
structure-to-fluid reaction, RHS coupling, Stage 14 RHS injection, or production
runtime hook insertion. Production commit and production structure advance gates
remain disabled.

## Safety boundary

Stage 20.3 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 or Stage 20.0-20.2 file. It does not run MPI, DNS, previous
stages, builds, or production validation.

## Next stage

Stage 20.4: structure-to-fluid reaction force candidate.

## Manual command

```bash
stage20_checks/run_stage20_3_structure_advance_hydro_force_candidate.sh
```

Expected output includes:

```text
STAGE 20.3 STRUCTURE ADVANCE HYDRO FORCE CANDIDATE VERDICT: PASS
STAGE 20.3 FINAL VERDICT: PASS
```
