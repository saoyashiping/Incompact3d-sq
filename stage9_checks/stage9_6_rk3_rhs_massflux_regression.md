# Stage 9.6 RK3/RHS/Mass-flux regression

Stage 9.6 adds a **diagnostic-only** gate for production channel DNS with no fibre/coupling.

- Goal: verify RK3 outer-step progression, RHS finiteness, velocity finiteness, CFL finiteness/limit, and mass-flux proxy drift.
- Not Stage 9.5: Stage 9.5 checks pressure-projection divergence pairs; Stage 9.6 checks RK3/RHS/mass-flux behavior.
- Real path exercised: production `xcompact3d` time loop (`calculate_transeq_rhs`, RK loop, projection/correction, post-step checks).
- No-fibre/no-coupling: reports `stage9_6_no_fibre_coupling_status 1` and does not connect Stage 8/IBM/fibre feedback paths.

## Recorded quantities
- RHS finite status from real `drho1/dux1/duy1/duz1` arrays.
- Velocity finite status from real `ux1/uy1/uz1` arrays.
- Per-step CFL proxy: `max(abs(u))*dt/dx`, `max(abs(v))*dt/dy`, `max(abs(w))*dt/dz` (MPI-global max).
- Mass-flux proxy: MPI-global domain mean of `ux1` at completed outer steps.
- Mass-flux tolerance: step-to-step absolute drift threshold `X3D_STAGE9_6_MASS_FLUX_TOL`.

## Criteria
PASS requires requested mode active, channel case, no coupling, >=1 completed steps, RK3 execution observed, finite RHS/velocity/CFL/mass-flux, CFL <= `X3D_STAGE9_6_CFL_MAX_LIMIT`, mass-flux drift <= tolerance, and finalise reached.

## Not tested yet
- No fibre coupling behavior, IBM feedback, or structural solver integration.
- Not a numerical-accuracy research benchmark.

## Manual run
`bash stage9_checks/run_stage9_6_rk3_rhs_massflux_regression.sh`

Expected lines include:
- `STAGE 9.6 RK3 RHS MASS-FLUX REGRESSION VERDICT: PASS`
- `STAGE 9.6 FINAL VERDICT: PASS`
