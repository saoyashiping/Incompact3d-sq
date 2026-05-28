# Stage 11.6 Oneway Sampling Invariance (np=1)

## Target

Compare np=1 deterministic no-fibre baseline against Stage 11 one-way hook-enabled run.

## Mathematical / physical meaning

Read-only sampling path: `U_f = I_h[u](X_f)` with:

- `f_fsi = 0`
- `RHS_stage11.6 = RHS_stage10 = RHS_stage9`
- invariance requirement: `Delta u = 0`

## Compared signatures

- `sum_ux`, `sum_uy`, `sum_uz`
- `max_ux`, `max_uy`, `max_uz`
- `l2_ux`, `l2_uy`, `l2_uz`

## Tolerance formula

`abs(delta) <= max(abs_tol, rel_tol * max(1, abs(reference)))`

## Hook evidence

- sampled velocity finite
- sample performed
- no field modification
- no RHS modification
- no IBM/feedback/two-way/structure advance

## Intentionally not done

- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no structure advance
- no Stage 10 script reruns
- no closed-stage edits

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE11_6_RUN_STAGE11_5=0 \
STAGE11_6_INVARIANCE_ABS_TOL=1.0e-12 \
STAGE11_6_INVARIANCE_REL_TOL=1.0e-14 \
bash stage11_checks/run_stage11_6_oneway_sampling_invariance_np1.sh
```
