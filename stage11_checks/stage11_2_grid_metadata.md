# Stage 11.2 Grid Metadata Gate

## Stage 11.2 target

Stage 11.2 introduces a grid metadata bridge skeleton for future one-way interpolation.

## Mathematical / physical meaning

Future target relation:

\[
U_f = I_h[u](X_f)
\]

Stage 11.2 defines grid metadata only and does not compute `U_f`.

- \(f_{fsi}=0\)
- \(RHS_{stage11.2} = RHS_{stage10} = RHS_{stage9}\)

## Metadata variables

- global size
- local bounds
- physical extents
- `dx/dy/dz`
- periodic flags
- staggered layout policy
- nonuniform-y policy

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

1. `fibre_stage11_grid_metadata_check` builds and runs.
2. Log contains `STAGE 11.2 GRID METADATA VERDICT: PASS`.
3. `stage11_outputs/fibre_stage11_2_grid_metadata.dat` contains required keys with value `1`.
4. `stage11_outputs/stage11_2_grid_metadata_gate.dat` reports `stage11_2_gate_status 1`.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE11_2_RUN_STAGE11_1=0 \
bash stage11_checks/run_stage11_2_grid_metadata.sh
```
