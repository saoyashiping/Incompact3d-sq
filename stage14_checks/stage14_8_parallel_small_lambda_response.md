# Stage 14.8 — Parallel Small-Lambda RHS-Response Consistency

## Objective

Stage 14.8 adds the evidence path for small-lambda Stage 14 RHS-response consistency across MPI decompositions `np=1`, `np=2`, and `np=4`.

For every decomposition, the gate compares:

- a lambda-zero production run with Stage 11, Stage 12, Stage 13, and Stage 14 enabled;
- a very small nonzero-lambda production run with the same Stage 11/12/13/14 path enabled.

The script records all run logs and diagnostic snapshots under `stage14_outputs/` and fails closed if any required log, diagnostic file, or contamination/status key is missing.

## Formula and physical meaning

The controlled RHS update remains:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

where `f_i_cand` is the already-audited Stage 13 Eulerian force-density candidate with the fibre-on-fluid sign convention.

For each `np` in `{1, 2, 4}`, Stage 14.8 evaluates the evidence for:

```text
Delta_RHS_np = RHS(lambda_14, np) - RHS(0, np)
Delta_RHS_np / lambda_14 ~= f_i_cand_np
```

Stage 14.8 does not reimplement the Stage 13 spreading path. It consumes the existing Stage 13 and Stage 14 production diagnostics.

## Run matrix

The script runs the following cases directly through `xcompact3d`:

| Decomposition | Baseline | Candidate |
| --- | --- | --- |
| `np=1` | `lambda_14=0` | `lambda_14=STAGE14_8_SMALL_LAMBDA` |
| `np=2` | `lambda_14=0` | `lambda_14=STAGE14_8_SMALL_LAMBDA` |
| `np=4` | `lambda_14=0` | `lambda_14=STAGE14_8_SMALL_LAMBDA` |

The default small lambda is intentionally tiny:

```text
STAGE14_8_SMALL_LAMBDA=1.0e-8
```

The script rejects non-positive lambda values and values larger than `1.0e-4`.

## Required artifacts

For each `np`, the gate writes or copies:

```text
stage14_outputs/stage14_8_lambda0_np${np}.dat
stage14_outputs/stage14_8_small_lambda_np${np}.dat
stage14_outputs/stage14_8_lambda0_rhs_hook_np${np}.dat
stage14_outputs/stage14_8_small_lambda_rhs_hook_np${np}.dat
stage14_outputs/stage14_8_lambda0_stage13_force_density_np${np}.dat
stage14_outputs/stage14_8_small_lambda_stage13_force_density_np${np}.dat
stage14_outputs/stage14_8_lambda0_np${np}.log
stage14_outputs/stage14_8_small_lambda_np${np}.log
```

The final Stage 14.8 gate file is:

```text
stage14_outputs/stage14_8_parallel_small_lambda_response.dat
```

## Pass/fail criteria

The gate passes only if all of the following are true.

### Lambda-zero evidence for np=1/2/4

- Stage 14 hook diagnostics exist.
- Stage 14 hook is initialized and applied.
- `stage14_5_lambda_zero_status = 1`.
- RHS increment L2 and max-abs are within zero tolerance.
- Stage 9.9 final fluid signatures are finite.
- No pressure/projection/Poisson/RK3/channel-forcing modification status is missing or false.
- No production IBM forcing, feedback application, two-way force, or structure advance status is missing or false.

### Small-lambda evidence for np=1/2/4

- Stage 14 hook diagnostics exist.
- Stage 14 hook is initialized and applied.
- `stage14_5_lambda_zero_status = 0`.
- The recorded `stage14_5_injection_gain` matches `STAGE14_8_SMALL_LAMBDA`.
- RHS increment L2 and max-abs are nonzero, finite, and bounded.
- Stage 13 sign diagnostics confirm the fibre-on-fluid candidate path and wrong-sign rejection.
- Stage 9.9 final fluid-signature deltas versus the corresponding lambda-zero run are finite and bounded.
- No pressure/projection/Poisson/RK3/channel-forcing modification status is missing or false.
- No production IBM forcing, feedback application, two-way force, or structure advance status is missing or false.

If the small-lambda hook diagnostics are absent or still show zero-lambda/zero-increment behavior, the script fails explicitly with `nonzero_lambda_still_blocked_by_stage14_5_np${np}`.

### Cross-decomposition evidence

- Normalized RHS increment diagnostics `stage14_5_rhs_increment_l2 / lambda_14` and `stage14_5_rhs_increment_max_abs / lambda_14` are consistent across `np=1`, `np=2`, and `np=4`.
- Stage 13 force-density diagnostics are consistent across `np=1`, `np=2`, and `np=4`.
- Final fluid signatures for lambda-zero and small-lambda cases remain within inherited Stage 9/Stage 13/Stage 14 tolerances.

## No-contamination boundaries

Stage 14.8 intentionally does not:

- modify source files;
- modify closed Stage 10, Stage 11, Stage 12, Stage 13, or Stage 14.0–14.7 files;
- reimplement Stage 13 Eulerian force-density spreading;
- alter pressure projection, Poisson, RK3 coefficients, channel forcing, or DNS numerics;
- call production IBM forcing outside the existing diagnostic path;
- advance fibre structure, bending, tension, wall/contact, or two-way structural coupling.

## Manual command

Run the gate manually from the repository root:

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_8_SMALL_LAMBDA=1.0e-8 \
STAGE14_8_MAX_RHS_INCREMENT=1.0e-4 \
STAGE14_8_MAX_FLUID_DELTA=1.0e-4 \
bash stage14_checks/run_stage14_8_parallel_small_lambda_response.sh
```

Expected successful verdict lines:

```text
STAGE 14.8 PARALLEL SMALL LAMBDA RESPONSE VERDICT: PASS
STAGE 14.8 FINAL VERDICT: PASS
```
