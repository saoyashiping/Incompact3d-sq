# Stage 9.11 total smoke

Stage 9.11 is the Stage 9 total-smoke closure gate for no-fibre production DNS regression.

## Mathematical/physical model
No new physics is added:

- `div(u)=0`
- `du/dt + u·grad(u) = -grad(p) + nu laplacian(u) + f_drive`
- `F_coupling = 0`

## Stage 9.11 intentionally does not do
- no IBM force;
- no fibre structure advance;
- no feedback force;
- no two-way coupling;
- no Stage 10 production coupling hook.

## Pass criteria
- `xcompact3d` build PASS;
- Stage 9.1–9.10 all PASS;
- no-coupling status PASS;
- `stage9_checks/STAGE9_CLOSED.md` generated only on full PASS.

## Current known blocker
If Stage 9.8 is in known regression state, Stage 9.11 must fail honestly (blocked by Stage 9.8) until Stage 9.8 is repaired in a separate step.

## Manual command
```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE9_11_TIMEOUT_SEC=240 \
bash stage9_checks/run_stage9_11_total_smoke.sh
```
