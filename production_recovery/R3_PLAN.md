# Production Recovery R3 Plan

## R3 scope

R3 adds a standalone production fibre grid adapter that describes Xcompact3D-style grid coordinates, local pencil ranges, periodic flags, grid spacing, and local cell volumes for later fibre-system use.

## R3 non-goals

R3 does not connect to `xcompact3d.f90`, does not enter the production time loop, does not perform IBM interpolation, does not perform IBM spreading, does not inject RHS coupling, does not advance structure, does not implement wall contact, and does not implement fibre-fibre collision.

## Files created/modified

Created:

* `src/fibre_prod_grid_adapter.f90`
* `src/fibre_prod_grid_adapter_check.f90`
* `production_recovery/R3_PLAN.md`
* `production_recovery/R3_BUILD_LOG.txt`
* `production_recovery/R3_RUN_LOG.txt`
* `production_recovery/R3_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R3_PASS_FAIL.md`
* `production_recovery/R3_evidence/README.md`

Modified:

* `src/CMakeLists.txt`
* R2 evidence documents listed in `R3_SOURCE_DIFF_SUMMARY.md`

## Build strategy

Add and build only the standalone `fibre_prod_grid_adapter_check` target from `fibre_prod_grid_adapter.f90` and `fibre_prod_grid_adapter_check.f90`.

## Run strategy

Run only the standalone `fibre_prod_grid_adapter_check` executable. Success requires real output containing `R3_FIBRE_PROD_GRID_ADAPTER_CHECK PASS`.

## Pass/fail criteria

PASS requires the standalone target to compile, run, and print `R3_FIBRE_PROD_GRID_ADAPTER_CHECK PASS` while checking coordinate validity, local ranges, positive spacing, positive cell volume, point lookup, out-of-domain failure, and destroy/deallocation.

FAIL applies if the standalone target builds and runs but any required adapter check fails or R3 violates forbidden-scope boundaries.

BLOCKED applies if external environment prerequisites prevent building or running the standalone target.

## Evidence boundary

R3 PASS only means the grid adapter independent validation passed. It does not mean IBM interpolation, IBM spreading, RHS coupling, structure advancement, wall contact, fibre-fibre collision, MPI, DNS, FSI, or production DNS-FSI closure passed.
