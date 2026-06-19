# R3 Source Diff Summary

## Added files

* `src/fibre_prod_grid_adapter.f90`
* `src/fibre_prod_grid_adapter_check.f90`
* `production_recovery/R3_PLAN.md`
* `production_recovery/R3_BUILD_LOG.txt`
* `production_recovery/R3_RUN_LOG.txt`
* `production_recovery/R3_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R3_PASS_FAIL.md`
* `production_recovery/R3_evidence/README.md`
* `production_recovery/PRODUCTION_RECOVERY_R2_CLOSED.md`

## Modified files

* `src/CMakeLists.txt`
* `production_recovery/R2_RUN_LOG.txt`
* `production_recovery/R2_PASS_FAIL.md`
* `production_recovery/R2_SOURCE_DIFF_SUMMARY.md`
* `PRODUCTION_RECOVERY_STATUS.md`

## `src/xcompact3d.f90` modified?

No.

## `src/CMakeLists.txt` modified?

Yes. Only the standalone `fibre_prod_grid_adapter_check` target was added for R3.

## Connected to `xcompact3d` executable?

No. `fibre_prod_grid_adapter.f90` was not added to the `xcompact3d` executable source list and is not connected to the main time loop.

## IBM/RHS/FSI/structure/contact/collision implementation?

No. R3 does not implement IBM interpolation, IBM spreading, RHS coupling, FSI coupling, structure advancement, wall contact, or fibre-fibre collision.
