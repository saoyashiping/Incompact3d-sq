# Stage 9.5 pressure projection regression gate

## Goal
Audit real no-fibre production channel pressure-projection behavior over a few complete outer steps.

## Why still no-fibre/no-coupling
This is diagnostics-only instrumentation around existing production projection path. No fibre advance/coupling/RHS injection is connected.

## Difference from Stage 9.4
- Stage 9.4: run continuity smoke for few steps.
- Stage 9.5: specifically records divergence before/after projection and checks tolerance/finite regression.

## Projection path exercised
Real production path in `xcompact3d` around:
- `calc_divu_constraint(divu3,...)` before projection,
- `solve_poisson(...)`,
- `cor_vel(...)`,
- `calc_divu_constraint(divu3,...)` after correction.

## Recorded divergence quantities
Per completed step:
- div_before_max / div_before_mean
- div_after_max / div_after_mean
- div_reduction_ratio_max / div_reduction_ratio_mean

## Finite-field criteria
Checks finite status for pressure/projection fields (`pp3`,`px1`,`py1`,`pz1`), corrected velocity (`ux1`,`uy1`,`uz1`), and divergence snapshots.

## Tolerance criteria
From env:
- `X3D_STAGE9_5_DIV_MAX_TOL` (default `1.0e-8`)
- `X3D_STAGE9_5_DIV_MEAN_TOL` (default `1.0e-9`)
Pass requires post-projection divergence meeting both tolerances for sampled steps.

## Not tested yet
Long horizon accuracy/stability and any fibre-coupled physics remain out of scope.

## Pass criteria
- np=1/2/4 logs contain `STAGE 9.5 PRESSURE PROJECTION REGRESSION VERDICT: PASS`
- np dat files contain `stage9_5_projection_regression_status 1`
- script prints `STAGE 9.5 FINAL VERDICT: PASS`

## Manual run
```bash
cd /home/sq/Incompact3d-sq-Xcompact3D
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE9_5_MAX_STEPS=3 \
STAGE9_5_DIV_MAX_TOL=1.0e-8 \
STAGE9_5_DIV_MEAN_TOL=1.0e-9 \
bash stage9_checks/run_stage9_5_projection_regression.sh
```

## Expected PASS lines
- `STAGE 9.5 PRESSURE PROJECTION REGRESSION VERDICT: PASS`
- `STAGE 9.5 FINAL VERDICT: PASS`


## Machine-readable data requirement
Stage 9.5 closure requires `.dat` per-step divergence values to match the real production projection diagnostics (`DIV U*` / `DIV U`) within normal floating-point formatting differences.


## Data source and compile-safety notes
- Stage 9.5 records real `DIV U*` and `DIV U` diagnostics directly from `navier::divergence` (same `tmax` and global mean values as production log prints).
- Pipe bulk/volume-average routines are not part of Stage 9.5 divergence capture and must not import Stage 9.5 diagnostic module symbols.
- Machine-readable `.dat` per-step divergence values must be non-placeholder and match production log diagnostics within floating-point formatting differences.

- `pipe_volume_avg` must not use a dummy argument named `var` when importing module `var`; it uses `field` to avoid ambiguity.
