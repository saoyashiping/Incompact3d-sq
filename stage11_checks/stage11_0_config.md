# Stage 11.0 Configuration Gate

## Stage 11.0 target

- Stage 11 one-way hook configuration
- global switches
- readonly-only mode

## Mathematical / physical meaning

- future one-way relation: `U_f = I_h[u](X_f)`
- Stage 11.0 does not compute `U_f`
- `f_fsi = 0`
- `RHS_stage11.0 = RHS_stage10 = RHS_stage9`

## Environment variables

- `X3D_STAGE11_ONEWAY_HOOK`
- `X3D_STAGE11_FORCE_READONLY`
- `X3D_STAGE11_MAX_POINTS`
- `X3D_STAGE11_MAX_STEPS`

## Intentionally not done

- no `xcompact3d` hook call
- no production main-loop insertion
- no Lagrangian point allocation
- no interpolation
- no fluid-field read
- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance

## Pass criteria

- required build targets succeed
- `fibre_stage11_config_check` prints `STAGE 11.0 CONFIG VERDICT: PASS`
- all required dat keys are `1`
- `stage11_outputs/stage11_0_config_gate.dat` final status is `1`

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE11_0_RUN_STAGE10_CLOSURE=0 \
bash stage11_checks/run_stage11_0_config.sh
```
