# Stage 20.2 fluid-to-structure force input adapter

Stage 20.2 is a helper-local, diagnostic-only adapter audit. It constructs a
controlled candidate `F_fs_candidate` from helper-local fibre velocity and
controlled fluid-point velocity arrays. It does not activate production DNS
fluid-force input, structure-to-fluid reaction, Eulerian RHS spreading, Stage 14
RHS injection, IBM, DNS-core, projection, Poisson, RK3/channel forcing, or
production restart/statistics/visualization I/O.

## Source-only and previous-stage policy

Stage 20.2 accepts Stage 20.1 PASS evidence when present and preserves Stage
20.0/20.1 source-only acceptance behavior. It does not rerun Stage 10-19 or
Stage 20.0-20.1, does not require old closure files, and does not require all
old stage output directories in source-only archives.

## Helper-local fields

Required helper-local arrays are `X_fibre`, `V_fibre`, `u_interp_candidate`,
`u_relative_candidate`, `F_fs_candidate`, `F_b_candidate`, `F_T_candidate`,
`F_total_structure_candidate_without_fluid`,
`F_total_structure_candidate_with_fluid`, `owner_rank`, `global_point_id`, and
`local_point_id`.

The audited sign convention is:

```text
u_relative_candidate = u_interp_candidate - V_fibre
F_fs_candidate = C_fs * u_relative_candidate
F_total_structure_candidate_without_fluid = F_b_candidate + F_T_candidate
F_total_structure_candidate_with_fluid = F_b_candidate + F_T_candidate - F_fs_candidate
```

## Gate distinction

`STAGE20_2_FLUID_TO_STRUCTURE_ENABLE=1` is allowed only for helper-local
candidate `F_fs_candidate` construction. It does not imply two-way coupling,
structure-to-fluid reaction, RHS coupling, Stage 14 RHS injection, or production
runtime hook insertion. `STAGE20_2_TWOWAY_COUPLING_ENABLE`,
`STAGE20_2_STRUCTURE_TO_FLUID_ENABLE`, and `STAGE20_2_RHS_COUPLING_ENABLE`
remain false, while `STAGE20_2_LAMBDA_COUPLING` remains zero.

## Safety boundary

Stage 20.2 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 or Stage 20.0-20.1 file. It does not run MPI, DNS, previous
stages, builds, or production validation.

## Next stage

Stage 20.3: structure advance with hydrodynamic force candidate.

## Manual command

```bash
stage20_checks/run_stage20_2_fluid_to_structure_force_input_adapter.sh
```

Expected output includes:

```text
STAGE 20.2 FLUID TO STRUCTURE FORCE INPUT ADAPTER VERDICT: PASS
STAGE 20.2 FINAL VERDICT: PASS
```
