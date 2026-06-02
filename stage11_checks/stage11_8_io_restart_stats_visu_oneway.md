# Stage 11.8 I/O + Restart + Stats/Visu/Coarse Compatibility

## Stage 11.8 target

Validate restart / statistics / visualization / coarse I/O compatibility under the Stage 11 one-way hook.

## Mathematical / physical meaning

- `U_f = I_h[u](X_f)` is read-only sampling.
- `f_fsi = 0`.
- `RHS_stage11.8 = RHS_stage10 = RHS_stage9`.
- `Delta restart_state = 0` and `Delta stats_visu_coarse_io = 0`.

## Reused gates

- Stage 9.7 stats / visu / coarse I/O smoke
- Stage 9.8 restart I/O regression

## Hook evidence

- sampled velocity finite
- sample performed
- no field modification
- no RHS modification
- no IBM/feedback/two-way/structure advance

## What is intentionally not done

- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance
- no Stage 10 script reruns
- no closed-stage edits

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE11_8_RUN_STAGE11_7=0 \
bash stage11_checks/run_stage11_8_io_restart_stats_visu_oneway.sh
```
