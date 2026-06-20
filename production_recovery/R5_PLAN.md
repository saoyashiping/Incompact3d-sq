# Production Recovery R5 Plan

## R5 scope

R5 adds a standalone dry flexible-fibre structure solver foundation for production fibre state. It covers free-free endpoint handling, bending force diagnostics, tension/length diagnostics, explicit dry-structure stepping, finite-state checks, and energy diagnostics.

## R5 non-goals

R5 does not connect to `xcompact3d.f90`, does not call fluid velocity, does not call IBM interpolation, does not implement IBM spreading, does not write RHS, does not implement two-way FSI, does not implement wall contact, and does not implement fibre-fibre collision.

## Files created/modified

Created:

* `src/fibre_prod_boundary_conditions.f90`
* `src/fibre_prod_bending_solver.f90`
* `src/fibre_prod_tension_solver.f90`
* `src/fibre_prod_structure_solver.f90`
* `src/fibre_prod_structure_solver_check.f90`
* `production_recovery/R5_PLAN.md`
* `production_recovery/R5_BUILD_LOG.txt`
* `production_recovery/R5_RUN_LOG.txt`
* `production_recovery/R5_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R5_PASS_FAIL.md`
* `production_recovery/R5_evidence/README.md`

Modified:

* `src/CMakeLists.txt`

## Build strategy

Build only the standalone `fibre_prod_structure_solver_check` target from the R5 dry-structure modules plus the existing production fibre state modules.

## Run strategy

Run only the standalone `fibre_prod_structure_solver_check` executable. Success requires real output containing `R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK PASS`.

## Pass/fail criteria

PASS requires the standalone target to compile, run, and print `R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK PASS` while checking allocation, straight initialization, near-zero straight bending force, no-force invariance, finite state after step, external-force response direction, finite segment/stretch diagnostics, non-negative energy diagnostics, fail-closed invalid `dt`/`ds`, and destroy/deallocation.

FAIL applies if the standalone target builds and runs but any required structure-solver check fails or R5 violates forbidden-scope boundaries.

BLOCKED applies if external environment prerequisites prevent building or running the standalone target.

## Evidence boundary

R5 PASS only means the production dry structure solver passed independent validation. It does not mean IBM interpolation, IBM spreading, RHS coupling, two-way FSI, wall contact, fibre-fibre collision, MPI, DNS, or production DNS-FSI closure passed.
