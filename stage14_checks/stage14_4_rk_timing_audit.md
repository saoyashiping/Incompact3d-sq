# Stage 14.4 — RK Timing Audit

## Target

Stage 14.4 audits the planned RK-substep timing for future controlled RHS injection. It verifies, in standalone/static arrays only, that a future injection point is planned once per intended RK substep and that the abstract event order is:

```text
Stage 13 candidate ready
Stage 14 RHS increment ready
RHS addition planned
projection planned
```

No production RHS addition is performed in Stage 14.4.

## Mathematical / physical meaning

Stage 13 provides the Eulerian force-density candidate:

```text
f_i_cand = S_h[F_fibre_to_fluid_cand]
```

Stage 14.1 provides the standalone accumulator:

```text
f_rhs_stage14 = lambda_14 * f_i_cand
```

Stage 14.2 verifies the standalone additive formula:

```text
RHS_new = RHS_old + f_rhs_stage14
```

Stage 14.3 verifies the fibre-on-fluid sign and lambda scaling.

Stage 14.4 verifies timing only. The future production operation must be ordered as:

```text
build / update Stage 13 force-density candidate
compute Stage 14 increment lambda_14 * f_i_cand
add increment to RHS once per intended RK substep
pressure / projection correction
```

At Stage 14.4, production FSI remains disabled:

```text
f_fsi = 0
RHS_stage14.4 = RHS_stage14.3 = RHS_stage14.2 = RHS_stage14.1 = RHS_stage14.0 = RHS_stage13 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Controlled checks

The Stage 14.4 standalone check covers:

- valid RK3 event sequence;
- Stage 13 candidate-before-increment order;
- Stage 14 increment-before-RHS-addition order;
- RHS-addition-before-projection order;
- once-per-substep planned injection;
- duplicate injection detection;
- skipped substep detection;
- wrong projection-before-RHS order detection;
- invalid event/substep rejection;
- `lambda=0` zero increment through the Stage 14.1 standalone accumulator path.

## What is intentionally not done

Stage 14.4 intentionally performs no production coupling. In particular, it does not:

- add an `xcompact3d` hook call;
- insert code in the production main loop;
- add any increment to production RHS arrays;
- call production IBM forcing;
- modify pressure, projection, Poisson, RK3, or channel forcing code;
- access or modify production fluid fields;
- advance the fibre structure.

## Pass criteria

The check passes only if all required status keys are `1`, the valid sequence records exactly three planned RHS additions, the valid-sequence duplicate/skipped/wrong-order counts remain zero, and the lambda-zero RHS increment is no larger than `1.0e-12`.

The shell gate writes:

```text
stage14_outputs/stage14_4_rk_timing_audit_gate.dat
```

and the standalone executable writes:

```text
stage14_outputs/fibre_stage14_4_rk_timing_audit.dat
```

Both verdicts must report PASS.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE14_4_RUN_STAGE14_3=0 \
bash stage14_checks/run_stage14_4_rk_timing_audit.sh
```
