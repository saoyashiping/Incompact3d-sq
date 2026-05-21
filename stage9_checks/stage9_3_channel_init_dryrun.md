# Stage 9.3 channel initialization dry-run gate

## Goal

Stage 9.3 adds a production **initialization dry-run** gate for channel case that exercises real xcompact3d startup and stops before any real timestep/RK3/RHS/projection advance.

## Why this remains no-fibre/no-coupling

This stage does not connect Stage 8 coupling to production DNS loop. It only audits initialization state. No fibre advance, no RHS force injection, no production IBM coupling path is introduced.

## Initialization path exercised

With `X3D_STAGE9_3_CHANNEL_INIT_DRYRUN=1`, `xcompact3d` runs normal `init_xcompact3d()` path including MPI init, parameter read, decomp init, decomp I/O init, coarse stat/visu init, `ph1/ph2/ph3/ph4/phG` setup, arrays/schemes/poisson/channel init, then performs Stage 9.3 audit and exits before time loop.

## Not tested yet

- Real timestep advance loop
- RK3 / full RHS/projection iteration behavior
- Production fibre coupling injection
- Stage 7/8 physics validation (unchanged and out of this gate)

## Pass criteria

- each np run log contains `STAGE 9.3 CHANNEL INIT DRY-RUN VERDICT: PASS`
- `stage9_outputs/fibre_stage9_3_channel_init_dryrun.dat` exists
- dat file contains `stage9_3_channel_init_dryrun_status 1`
- gate script final line is `STAGE 9.3 FINAL VERDICT: PASS`

## Manual run

```bash
cd /home/sq/Incompact3d-sq-Xcompact3D
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
bash stage9_checks/run_stage9_3_channel_init_dryrun.sh
```
