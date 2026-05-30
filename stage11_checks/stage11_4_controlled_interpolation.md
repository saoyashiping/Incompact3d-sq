# Stage 11.4 Controlled Interpolation

Stage 11.4 target: controlled analytic interpolation checks.

- Verifies `U_f = I_h[u](X_f)` on synthetic fields only.
- `f_fsi=0`, so `RHS_stage11.4 = RHS_stage10 = RHS_stage9`.

Analytic tests: constant, linear, shear, Poiseuille-like, periodic x/z safety, near-wall y safety, out-of-domain rejection, weight-sum check.

Tolerances:
- constant <= 1e-12
- linear <= 1e-12
- shear <= 1e-12
- weight sum <= 1e-12
- Poiseuille-like <= 5e-3

Not done: no xcompact3d hook call, no production fluid read/modify, no RHS injection, no IBM/feedback/twoway, no structure advance.

Manual command:
```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE11_4_RUN_STAGE11_3=0 \
bash stage11_checks/run_stage11_4_controlled_interpolation.sh
```
