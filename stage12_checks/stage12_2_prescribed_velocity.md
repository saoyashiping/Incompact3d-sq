# Stage 12.2 Prescribed Velocity Skeleton

## Stage 12.2 target

- prescribed fibre/control-point velocity skeleton.

## Mathematical / physical meaning

- Future relation: `F_fs_cand = alpha * (U_f - V_f)`.
- Stage 12.2 only allocates/initializes `V_f`.
- Stage 12.2 does not compute force.
- `f_fsi = 0`.
- `RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9`.

## State variables

- `V_f(k,1:3)`
- `velocity_valid_flag(k)`
- `velocity_norm(k)`

## What is intentionally not done

- no xcompact3d hook call
- no production main-loop insertion
- no force candidate computation
- no Eulerian force density
- no RHS injection
- no IBM spreading
- no feedback force application
- no two-way force
- no fibre structure advance

## Pass criteria

Pass requires the Stage 12.2 prescribed-velocity check to allocate the Lagrangian placeholder velocity arrays, verify initial zero velocity, set a deterministic constant velocity, verify finite velocity norms and valid flags, clear back to zero, and prove no force/RHS/IBM/two-way/structure/fluid-field activity is enabled.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE12_2_RUN_STAGE12_1=0 \
bash stage12_checks/run_stage12_2_prescribed_velocity.sh
```
