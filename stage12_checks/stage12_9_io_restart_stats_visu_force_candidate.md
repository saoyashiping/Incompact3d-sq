# Stage 12.9 I/O restart stats/visu compatibility with force candidate

## Target

Stage 12.9 verifies restart, statistics, visualization, and coarse-I/O compatibility while the Stage 11 one-way sampling hook and Stage 12 production feedback-force candidate hook are enabled.

## Mathematical / physical meaning

The production no-fibre DNS remains unchanged:

```text
U_f = I_h[u](X_f)
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
f_fsi = 0
RHS_stage12.9 = RHS_stage12.8 = RHS_stage12.7 = RHS_stage12.6 = RHS_stage12.5 = RHS_stage12.4 = RHS_stage12.3 = RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
Delta restart_state = 0
Delta stats_visu_coarse_io = 0
```

Stage 12.9 allows diagnostic Lagrangian force-candidate and power-diagnostic computations, but performs no spreading, no RHS injection, and no modification to restart, statistics, visualization, or coarse I/O behavior.

## Reused gates

Stage 12.9 reuses these production smoke/regression paths with Stage 11 and Stage 12 hooks enabled:

```text
Stage 9.7 stats / visu / coarse I/O smoke
Stage 9.8 restart I/O regression
```

## Stage 12 diagnostic evidence

After both reused gates, `stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat` must prove:

- force candidate computed;
- finite force;
- finite force norm;
- finite power diagnostics;
- action-reaction consistency;
- pair-power consistency;
- no field modification;
- no RHS modification;
- no Eulerian force density;
- no RHS injection;
- no IBM spreading;
- no feedback application;
- no two-way force;
- no structure advance.

## Intentionally not done

- No Eulerian force density.
- No RHS injection.
- No IBM spreading.
- No feedback force application to fluid.
- No two-way force.
- No fibre structure advance.
- No Stage 10 script reruns.
- No closed-stage edits.

## Pass criteria

The Stage 12.9 gate passes only if the required build targets compile, Stage 9.7 reports `STAGE 9.7 FINAL VERDICT: PASS`, Stage 9.8 reports `STAGE 9.8 FINAL VERDICT: PASS`, Stage 12.6 diagnostics are present and passing after each reused gate, and the final gate dat reports no restart, stats, visualization, coarse-I/O, RHS, IBM, feedback, two-way, or structure contamination.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE12_9_RUN_STAGE12_8=0 \
bash stage12_checks/run_stage12_9_io_restart_stats_visu_force_candidate.sh
```
