# Stage 14.9 — Restart / Stats / Visu / Coarse I/O Compatibility under RHS Injection

## Objective

Stage 14.9 adds the evidence path for I/O and restart compatibility while the full diagnostic-to-RHS chain is enabled:

1. Stage 11 one-way sampling.
2. Stage 12 Lagrangian feedback-force candidate.
3. Stage 13 Eulerian force-density candidate.
4. Stage 14 controlled RHS injection with a very small nonzero `lambda_14`.

The stage verifies that the small-lambda RHS-injection path remains active and bounded while statistics output, visualization output, coarse-mesh/decomp I/O, and restart write/read behavior continue to work.

Stage 14.9 is an I/O and restart compatibility stage only. It does not introduce new physics and does not advance the fibre structure.

## Formula and physical meaning

The controlled update remains:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

where `f_i_cand` is the already-audited Stage 13 Eulerian force-density candidate with the fibre-on-fluid sign convention.

The default Stage 14.9 lambda is intentionally tiny:

```text
STAGE14_9_SMALL_LAMBDA=1.0e-8
```

The wrapper rejects non-positive values and values larger than `1.0e-4`.

## Run layout

The Stage 14.9 wrapper builds `xcompact3d` and then performs three production runs at `STAGE14_9_NP` ranks. The default is:

```text
STAGE14_9_NP=2
```

The runs are:

| Run | Purpose | Stage 9 evidence path |
| --- | --- | --- |
| stats/visu/coarse | Exercise statistics, visualization, and coarse/decomp I/O paths with Stage 14 RHS injection enabled. | Stage 9.7 |
| restart cold | Write restart/checkpoint files with Stage 14 RHS injection enabled. | Stage 9.8 cold phase |
| restart read | Read restart/checkpoint files and continue with Stage 14 RHS injection enabled. | Stage 9.8 restart phase |

`STAGE14_9_RUN_STAGE14_8=1` may optionally run the closed Stage 14.8 prerequisite gate first. The default is `0`, so no closed-stage script is run unless explicitly requested.

## Required artifacts

The wrapper writes the final gate file:

```text
stage14_outputs/stage14_9_io_restart_stats_visu_rhs_injection.dat
```

It also records/copies the evidence artifacts:

```text
stage14_outputs/stage14_9_stats_visu_coarse_io.log
stage14_outputs/stage14_9_stats_visu_coarse_io.dat
stage14_outputs/stage14_9_stats_visu_rhs_hook.dat
stage14_outputs/stage14_9_stats_visu_stage13_force_density.dat
stage14_outputs/stage14_9_restart_cold.log
stage14_outputs/stage14_9_restart_cold.dat
stage14_outputs/stage14_9_restart_cold_rhs_hook.dat
stage14_outputs/stage14_9_restart_cold_stage13_force_density.dat
stage14_outputs/stage14_9_restart_read.log
stage14_outputs/stage14_9_restart_read.dat
stage14_outputs/stage14_9_restart_read_rhs_hook.dat
stage14_outputs/stage14_9_restart_read_stage13_force_density.dat
stage14_outputs/stage14_9_restart_signature_np${STAGE14_9_NP}.dat
```

The wrapper fails closed if any required log, `.dat` file, restart/checkpoint file, Stage 13 diagnostic, Stage 14 diagnostic, or evidence key is missing.

## Pass/fail criteria

The gate passes only when all of the following are true.

### Stage 14 RHS-injection evidence

- Stage 14 hook diagnostics exist for the stats/visu/coarse run, restart cold run, and restart-read continuation run.
- Stage 14 hook is initialized and applied.
- `lambda_14` is recorded, finite, nonzero, and equal to `STAGE14_9_SMALL_LAMBDA` within a strict small-lambda tolerance.
- RHS increment L2 and max-abs diagnostics are finite and nonzero.
- RHS increment max-abs is bounded by `STAGE14_9_MAX_RHS_INCREMENT`.
- Stage 13 diagnostics confirm the fibre-on-fluid sign convention and wrong-sign rejection.

### Statistics compatibility

- Stage 9.7 statistics path is executed.
- Expected statistics output is reported present and nonempty.
- Statistics diagnostics are finite and report no NaN/Inf.
- Statistics I/O does not disable or corrupt Stage 13 or Stage 14 diagnostics.

### Visualization compatibility

- Stage 9.7 visualization path is executed.
- Expected visualization output is reported present and nonempty.
- Visualization field diagnostics are finite and report no NaN/Inf.
- Visualization I/O does not disable or corrupt Stage 13 or Stage 14 diagnostics.

### Coarse/decomp I/O compatibility

- Stage 9.7 coarse-mesh/decomp I/O descriptor, open, write, and close statuses are present and true.
- Coarse/decomp I/O output is nonempty and finite.
- Coarse/decomp I/O does not disable or corrupt Stage 13 or Stage 14 diagnostics.

### Restart compatibility

- Stage 9.8 cold phase writes restart/checkpoint files.
- Restart/checkpoint files exist and are nonempty before the restart-read phase.
- Stage 9.8 restart phase reads restart/checkpoint files and completes.
- Restart signature status is true.
- Restart signature deltas are finite and bounded by `STAGE14_9_MAX_RESTART_DELTA`.
- Restart-derived fluid response deltas are bounded by `STAGE14_9_MAX_FLUID_DELTA`.
- Stage 14 diagnostics remain active after restart continuation.

### No-contamination boundaries

The gate requires diagnostic evidence for:

- no pressure/projection modification;
- no Poisson modification;
- no RK3/channel-forcing modification;
- no production IBM forcing;
- no feedback application beyond diagnostics;
- no two-way force;
- no fibre structure advance.

Stage 14.9 also records no-bending, no-tension, no-fibre-position-update, no-fibre-velocity-structural-update, and no-wall/contact statuses as Stage 14.9 no-structure-advance boundary evidence. The wrapper does not call or modify any structure solver.

## Manual command

Run from the repository root:

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_9_SMALL_LAMBDA=1.0e-8 \
STAGE14_9_MAX_RHS_INCREMENT=1.0e-4 \
STAGE14_9_MAX_FLUID_DELTA=1.0e-4 \
STAGE14_9_MAX_RESTART_DELTA=1.0e-8 \
STAGE14_9_MAX_IO_SIGNATURE_DELTA=1.0e-8 \
STAGE14_9_NP=2 \
STAGE14_9_RUN_STAGE14_8=0 \
bash stage14_checks/run_stage14_9_io_restart_stats_visu_rhs_injection.sh
```

The script may also be run with defaults exactly as:

```sh
bash stage14_checks/run_stage14_9_io_restart_stats_visu_rhs_injection.sh
```

## Expected verdict

A successful run prints:

```text
STAGE 14.9 IO/RESTART/STATS/VISU RHS-INJECTION VERDICT: PASS
STAGE 14.9 FINAL VERDICT: PASS
```
