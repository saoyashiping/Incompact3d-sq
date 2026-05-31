# Stage 12.1 Lagrangian Force Candidate Buffer

## Stage 12.1 target

- Lagrangian feedback force candidate buffer skeleton.

## Mathematical / physical meaning

- Future relation: `F_fs_cand = alpha * (U_f - V_f)`.
- Stage 12.1 only allocates/initializes a zero-valued `F_fs_cand` buffer.
- Stage 12.1 does not compute force.
- `f_fsi = 0`.
- `RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9`.

## State variables

- `F_fs_cand(k,1:3)`
- `force_valid_flag(k)`
- `force_norm(k)`

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

Pass requires the Stage 12.1 force-buffer check to allocate the Lagrangian placeholder arrays, keep all force entries and norms zero/finite, keep validity flags safe, clear the buffer, and prove no force/RHS/IBM/two-way/structure/fluid-field activity is enabled.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE12_1_RUN_STAGE12_0=0 \
bash stage12_checks/run_stage12_1_force_buffer.sh
```
