# Stage 14.1 — RHS Accumulator Skeleton

## Target

Stage 14.1 adds the RHS injection candidate accumulator / buffer skeleton for
future controlled Stage 14 forcing.  It can form the diagnostic standalone
quantity

```text
f_rhs_stage14 = lambda_14 * f_i_cand
```

from controlled Stage 13-like force-density candidates, but it does not add the
accumulator to production RHS.

## Mathematical / physical meaning

Stage 13 provides the Eulerian candidate force density `f_i_cand`.  A later
Stage 14 step may use:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

Stage 14.1 does not perform `RHS_new`.  It only forms the standalone diagnostic
accumulator:

```text
f_rhs_stage14 = lambda_14 * f_i_cand
```

At Stage 14.1:

```text
f_fsi = 0
RHS_stage14.1 = RHS_stage14.0 = RHS_stage13 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

No production RHS addition is performed and no structure equation is advanced.

## State variables

The accumulator owns only diagnostic standalone buffers:

```text
rhs_fx_cand(i,j,k)
rhs_fy_cand(i,j,k)
rhs_fz_cand(i,j,k)
rhs_increment_norm(i,j,k)
rhs_increment_valid_flag(i,j,k)
```

These buffers are not production RHS arrays and are not connected to
`xcompact3d`.

## Controlled checks

The Stage 14.1 standalone check verifies:

- zero initialization after allocation;
- clear behavior after a compute path;
- `lambda_14 = 0` gives a zero accumulator;
- `lambda_14 = 1` reproduces the candidate field;
- `lambda_14 = 0.1` produces the expected scaled field;
- component-wise `x/y/z` scaling;
- finite accumulator values and finite norms;
- valid flags remain in the set `{0, 1}`;
- invalid shapes are rejected without allocation.

## Intentionally not done

Stage 14.1 intentionally does not add:

- an `xcompact3d` hook call;
- a production main-loop insertion;
- production RHS addition;
- production IBM forcing application;
- pressure / projection / Poisson / RK3 / channel-forcing modification;
- feedback force application to fluid;
- two-way force application;
- fibre structure advance.

## Pass criteria

The Stage 14.1 gate passes only when:

1. `xcompact3d` and all required closed-stage check targets still build.
2. `fibre_stage14_config_check` still builds.
3. `fibre_stage14_rhs_accumulator_check` builds and prints PASS.
4. The diagnostic `.dat` file reports all required Stage 14.1 status keys as
   passing.
5. The lambda-zero, lambda-one, lambda-fractional, component, and post-clear
   scalar errors are all at or below `1.0e-12`.
6. All no-RHS-addition, no-hook, no-fluid-field-access, no-fluid-field-write,
   no-pressure/projection/RK3/channel-forcing, no-IBM, no-feedback, no-two-way,
   and no-structure statuses remain `1`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE14_1_RUN_STAGE14_0=0 \
bash stage14_checks/run_stage14_1_rhs_accumulator.sh
```
