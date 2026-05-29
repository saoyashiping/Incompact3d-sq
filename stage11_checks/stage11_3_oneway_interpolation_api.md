# Stage 11.3 Oneway Interpolation API Gate

## Stage 11.3 target

One-way interpolation API compile-only / no-call skeleton.

## Mathematical / physical meaning

Future target relation:

\[
U_f = I_h[u](X_f)
\]

Stage 11.3 defines the interpolation API only; it does not compute `U_f`.

- \(f_{fsi}=0\)
- \(RHS_{stage11.3} = RHS_{stage10} = RHS_{stage9}\)

## API roles

- init
- prepare
- guarded sample interface
- status getter
- diagnostics writer
- finalize

## What is intentionally not done

- no `xcompact3d` hook call
- no production main-loop insertion
- no interpolation
- no fluid-field read
- no fluid-field modification
- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance

## Pass criteria

1. `fibre_stage11_oneway_interpolation_check` builds and runs.
2. Log contains `STAGE 11.3 ONEWAY INTERPOLATION API VERDICT: PASS`.
3. `.dat` file contains required keys set to `1`.
4. Gate `.dat` file reports `stage11_3_gate_status 1`.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE11_3_RUN_STAGE11_2=0 \
bash stage11_checks/run_stage11_3_oneway_interpolation_api.sh
```
