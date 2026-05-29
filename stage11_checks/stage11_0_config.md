# Stage 11.0 Configuration Gate

## Target

Stage 11.0 introduces only the one-way hook configuration layer and global switches for a future fluid-to-fibre read-only path. This stage is readonly-only.

## Mathematical / physical meaning

Future one-way sampling target:

\[
U_f = I_h[u](X_f)
\]

In Stage 11.0, this is **not** computed. Also:

- \(f_{fsi}=0\)
- \(RHS_{stage11.0} = RHS_{stage10} = RHS_{stage9}\)

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
- no fluid-field modification
- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance

## Pass criteria

The stage passes when:

1. `fibre_stage11_config_check` builds and runs in serial;
2. the run log contains `STAGE 11.0 CONFIG VERDICT: PASS`;
3. `stage11_outputs/fibre_stage11_0_config.dat` contains required `stage11_0_*` keys with pass values;
4. `stage11_outputs/stage11_0_config_gate.dat` reports `stage11_0_gate_status 1`.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE11_0_RUN_STAGE10_CLOSURE=0 \
bash stage11_checks/run_stage11_0_config.sh
```
