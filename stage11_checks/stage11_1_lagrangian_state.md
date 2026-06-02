# Stage 11.1 Lagrangian State Gate

## Stage 11.1 target

Stage 11.1 introduces a minimal Lagrangian/control-point state skeleton for later one-way fluid-to-fibre sampling.

## Mathematical / physical meaning

Future target relation:

\[
U_f = I_h[u](X_f)
\]

Stage 11.1 defines only the `X_f` and `U_f` containers and does not compute interpolation. Also:

- \(f_{fsi}=0\)
- \(RHS_{stage11.1} = RHS_{stage10} = RHS_{stage9}\)

## State variables

- `X_f(k,1:3)`
- `U_f(k,1:3)`
- `valid_point_flag(k)`

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

1. `fibre_stage11_lagrangian_state_check` builds and runs.
2. Log contains `STAGE 11.1 LAGRANGIAN STATE VERDICT: PASS`.
3. `stage11_outputs/fibre_stage11_1_lagrangian_state.dat` contains required keys with value `1`.
4. `stage11_outputs/stage11_1_lagrangian_state_gate.dat` reports `stage11_1_gate_status 1`.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE11_1_RUN_STAGE11_0=0 \
bash stage11_checks/run_stage11_1_lagrangian_state.sh
```
