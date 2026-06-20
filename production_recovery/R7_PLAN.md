# Production Recovery R7 Plan

## R7 scope

R7 adds a standalone two-way FSI closed-loop micro-case that connects the R3 grid adapter, R4 IBM interpolation, R5 dry structure solver, R6 IBM spreading, a fluid surrogate, and FSI diagnostics for short-time force-exchange validation.

## R7 non-goals

R7 does not connect to `xcompact3d.f90`, does not modify Xcompact3D RK3/pressure/projection/channel forcing, does not write the real Navier-Stokes RHS, does not implement wall contact, does not implement fibre-fibre collision, does not validate restart/statistics/visualization, does not validate `np=1/2/4` MPI production consistency, and does not claim production DNS-FSI closure.

## Files created/modified

Created:

* `src/fibre_prod_fsi_config.f90`
* `src/fibre_prod_fluid_surrogate.f90`
* `src/fibre_prod_fsi_coupling.f90`
* `src/fibre_prod_fsi_diagnostics.f90`
* `src/fibre_prod_fsi_closed_loop_check.f90`
* `production_recovery/R7_PLAN.md`
* `production_recovery/R7_BUILD_LOG.txt`
* `production_recovery/R7_RUN_LOG.txt`
* `production_recovery/R7_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R7_PASS_FAIL.md`
* `production_recovery/R7_evidence/README.md`

Modified:

* `src/CMakeLists.txt`

## Build strategy

Build only the standalone `fibre_prod_fsi_closed_loop_check` target from the R2-R6 standalone modules plus R7 config, surrogate, coupling, diagnostics, and closed-loop check sources.

## Run strategy

Run only the standalone `fibre_prod_fsi_closed_loop_check` executable. Success requires real output containing `R7_FIBRE_PROD_FSI_CLOSED_LOOP_CHECK PASS`.

## Pass/fail criteria

PASS requires the standalone target to compile, run, and print `R7_FIBRE_PROD_FSI_CLOSED_LOOP_CHECK PASS` while checking lambda=0 no-contamination, small-lambda bounded response, action-reaction residual, finite short-time stability, and fail-closed invalid configuration/out-of-domain/NaN cases.

FAIL applies if the standalone target builds and runs but any required FSI closed-loop check fails or R7 violates forbidden-scope boundaries.

BLOCKED applies if external environment prerequisites prevent building or running the standalone target.

## Evidence boundary

R7 PASS only means the standalone two-way FSI closed-loop micro-case independent validation passed. It does not mean xcompact3d main-loop coupling, real DNS RHS coupling, wall contact, fibre-fibre collision, restart/statistics/visualization, `np=1/2/4` MPI production consistency, or production DNS-FSI final closure passed.
