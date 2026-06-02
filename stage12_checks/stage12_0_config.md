# Stage 12.0 Feedback Candidate Configuration

## Stage 12.0 target

- Stage 12 feedback force candidate configuration
- global switches
- readonly-only mode

## Mathematical / physical meaning

- Future relation: `F_fs_cand = alpha * (U_f - V_f)`.
- Stage 12.0 does not compute `F_fs_cand`.
- `f_fsi = 0`.
- `RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9`.

## Environment variables

- `X3D_STAGE12_FEEDBACK_CANDIDATE`
- `X3D_STAGE12_FORCE_READONLY`
- `X3D_STAGE12_FEEDBACK_GAIN`
- `X3D_STAGE12_FORCE_SIGN`
- `X3D_STAGE12_MAX_POINTS`

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

Pass requires the Stage 12.0 config check to build, run with controlled environment variables, write all required diagnostic keys, and prove readonly/no-force/no-coupling statuses remain enabled.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE12_0_RUN_STAGE11_CLOSURE=0 \
bash stage12_checks/run_stage12_0_config.sh
```
