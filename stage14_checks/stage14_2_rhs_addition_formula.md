# Stage 14.2 — RHS Addition Formula Checks

## Target

Stage 14.2 verifies the controlled RHS-addition formula in standalone analytic
arrays only.  It checks that

```text
RHS_new = RHS_old + f_rhs_stage14
```

with the Stage 14.1-style relation

```text
f_rhs_stage14 = lambda_14 * f_i_cand
```

while performing no production RHS addition.

## Mathematical / physical meaning

Stage 13 provides the Eulerian force-density candidate `f_i_cand`.  Stage 14.1
provides the diagnostic accumulator `f_rhs_stage14 = lambda_14 * f_i_cand`.
A future Stage 14 production step may use:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

Stage 14.2 verifies this formula only in controlled standalone arrays.  It does
not update production RHS arrays.  At Stage 14.2:

```text
f_fsi = 0
RHS_stage14.2 = RHS_stage14.1 = RHS_stage14.0 = RHS_stage13 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Controlled checks

The standalone Stage 14.2 check covers:

- `lambda_14 = 0` RHS invariance;
- `lambda_14 = 1` exact addition;
- `lambda_14 = 0.1` scaled addition;
- component-wise `x/y/z` addition;
- additive-not-overwrite behavior (`RHS_new - increment = RHS_old`);
- preservation of the old RHS arrays;
- finite old RHS, increment, and new RHS arrays;
- invalid-shape rejection before the valid formula checks.

## Intentionally not done

Stage 14.2 intentionally does not add:

- an `xcompact3d` hook call;
- a production main-loop insertion;
- production RHS addition;
- production IBM forcing application;
- pressure / projection / Poisson / RK3 / channel-forcing modification;
- feedback force application to fluid;
- two-way force application;
- fibre structure advance.

## Pass criteria

The Stage 14.2 gate passes only when:

1. `xcompact3d` and all required closed-stage check targets still build.
2. `fibre_stage14_config_check` still builds.
3. `fibre_stage14_rhs_accumulator_check` still builds.
4. `fibre_stage14_rhs_addition_formula_check` builds and prints PASS.
5. The diagnostic `.dat` file reports all required Stage 14.2 status keys as
   passing.
6. The lambda-zero, lambda-one, lambda-fractional, component, and old-RHS
   preservation scalar errors are all at or below `1.0e-12`.
7. All no-production-RHS-modification, no-hook, no-fluid-field-access,
   no-fluid-field-write, no-pressure/projection/RK3/channel-forcing, no-IBM,
   no-feedback, no-two-way, and no-structure statuses remain `1`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE14_2_RUN_STAGE14_1=0 \
bash stage14_checks/run_stage14_2_rhs_addition_formula.sh
```
