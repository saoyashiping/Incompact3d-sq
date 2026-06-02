# Stage 13.9 I/O, restart, stats, and visualization compatibility

## Target

Stage 13.9 verifies restart / statistics / visualization / coarse I/O
compatibility while the Stage 13 Eulerian force-density candidate diagnostic hook
is enabled together with the closed Stage 11 one-way sampling hook and Stage 12
Lagrangian feedback-force candidate hook.

The Stage 13 hook may compute diagnostic force-density data, but the production
fluid fields, restart state, statistics, visualization, coarse I/O, RHS,
projection, Poisson solve, RK3 path, channel forcing, IBM forcing, and fibre
structure must remain uncontaminated.

## Mathematical / physical meaning

The diagnostic force-density candidate path remains

```text
U_f = I_h[u](X_f)
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

No production forcing is applied:

```text
f_fsi = 0
RHS_stage13.9 = RHS_stage13.8 = RHS_stage13.7 = RHS_stage13.6
              = RHS_stage13.5 = RHS_stage13.4 = RHS_stage13.3
              = RHS_stage13.2 = RHS_stage13.1 = RHS_stage13.0
              = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
Delta restart_state = 0
Delta stats_visu_coarse_io = 0
```

The forbidden production update remains

```text
RHS <- RHS + f_i_cand
```

## Reused gates

Stage 13.9 reuses the Stage 9 production I/O gates as black-box runtime paths:

```text
Stage 9.7 stats / visu / coarse I/O smoke
Stage 9.8 restart I/O regression
```

Both gates run with Stage 11, Stage 12, and Stage 13 hook environment variables
enabled and `STAGE9_SKIP_PREREQS=1` so closed-stage prerequisite scripts are not
rerun by default.

## Stage 13 diagnostic evidence

After each Stage 9 gate, Stage 13.9 requires
`stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat` to show:

- force-density candidate computed;
- finite force-density;
- finite force-density norm;
- finite integrated force;
- integrated-force conservation;
- spreading input uses the fibre-on-fluid sign;
- wrong-sign rejection;
- no field modification;
- no RHS modification;
- no RHS injection;
- no production IBM forcing;
- no feedback application;
- no two-way force;
- no fibre structure advance.

## Intentionally not done

- no RHS injection;
- no production IBM forcing application;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance;
- no Stage 10 script reruns;
- no closed-stage edits.

## Pass criteria

The gate passes only if all required targets build, Stage 9.7 reports
`STAGE 9.7 FINAL VERDICT: PASS` under the Stage 11 + Stage 12 + Stage 13 hook
environment, Stage 9.8 reports `STAGE 9.8 FINAL VERDICT: PASS` under the same
hook environment, Stage 13 force-density diagnostics are present and valid after
both runs, and the final Stage 13.9 `.dat` reports no restart, stats / visu /
coarse I/O, RHS, IBM, feedback, two-way, or structure contamination.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE13_9_RUN_STAGE13_8=0 \
bash stage13_checks/run_stage13_9_io_restart_stats_visu_force_density.sh
```
