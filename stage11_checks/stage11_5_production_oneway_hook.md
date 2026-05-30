# Stage 11.5 Production Oneway Hook

Stage 11.5 target:
- guarded production one-way sampling hook
- read production velocity
- write sampled velocity diagnostics
- no force / no RHS modification

Math/physics:
- `U_f = I_h[u](X_f)` read-only data path
- `f_fsi = 0`
- `RHS_stage11.5 = RHS_stage10 = RHS_stage9`

Implementation policy:
- conservative production smoke sampling
- finite sampled velocity signatures
- no physical feedback yet
- no production grid-accurate interpolation claim yet

Intentionally not done:
- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance
- no pressure/projection/RK3/channel forcing modifications

Manual command:
```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE11_5_RUN_STAGE11_4=0 \
bash stage11_checks/run_stage11_5_production_oneway_hook.sh
```
