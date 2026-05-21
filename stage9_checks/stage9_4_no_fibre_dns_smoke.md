# Stage 9.4 no-fibre DNS smoke gate

## Goal
Run real production channel DNS beyond init for a few complete outer time steps (default 3) with no fibre/coupling injection.

## Why still no-fibre/no-coupling
This stage is diagnostics-only smoke on production loop continuity. It does not connect Stage 8 coupling routines and does not inject fibre forces into production RHS.

## Difference from Stage 9.3
- Stage 9.3: initialization dry-run, stops before time loop.
- Stage 9.4: executes real production outer steps, then exits after requested completed outer steps.

## Production path exercised
Initialization + real production time loop body (including RK substeps), pressure solve/correction, postprocessing, restart/stats flow for completed outer steps.

## Not tested yet
Long-time stability/accuracy, coupling physics validation, production fibre two-way coupling.

## Pass criteria
- np=1/2/4 logs contain `STAGE 9.4 NO-FIBRE DNS SMOKE VERDICT: PASS`
- np-specific dat files contain `stage9_4_no_fibre_dns_smoke_status 1`
- script prints `STAGE 9.4 FINAL VERDICT: PASS`

## Manual run
```bash
cd /home/sq/Incompact3d-sq-Xcompact3D
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE9_4_MAX_STEPS=3 \
bash stage9_checks/run_stage9_4_no_fibre_dns_smoke.sh
```

## Expected PASS lines
- `STAGE 9.4 NO-FIBRE DNS SMOKE VERDICT: PASS`
- `STAGE 9.4 FINAL VERDICT: PASS`
