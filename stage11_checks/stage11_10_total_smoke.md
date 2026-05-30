# Stage 11.10 Total Smoke / Closure

## Stage 11.10 target

- Stage 11 total smoke
- Stage 11 closure
- production one-way fluid-to-fibre sampling path closure

## Mathematical / physical meaning

- `U_f = I_h[u](X_f)`
- `f_fsi = 0`
- `RHS_stage11 = RHS_stage10 = RHS_stage9`
- no feedback force
- no two-way coupling

## Closure evidence

- Stage 11.5 production one-way hook smoke
- Stage 11.6 np=1 invariance
- Stage 11.7 np=1/2/4 parallel consistency
- Stage 11.8 restart / stats / visu / coarse I/O compatibility
- Stage 11.9 RHS / coupling contamination audit

## Explicit closed-stage warning

- do not modify Stage 10
- do not modify Stage 11.0–11.9
- do not modify production solver sources during Stage 11.10

## What is intentionally not done

- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance
- no new production physics

## Pass criteria

Pass requires all Stage 11.5–11.9 closure gates to report final PASS, hook and contamination diagnostics to match expected read-only statuses, and final Stage 11.10 status to pass.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
bash stage11_checks/run_stage11_10_total_smoke.sh
```
