# R4 Source Diff Summary

## Added files

* `src/fibre_prod_ibm_delta.f90`
* `src/fibre_prod_ibm_interpolation.f90`
* `src/fibre_prod_ibm_interpolation_check.f90`
* `production_recovery/R4_PLAN.md`
* `production_recovery/R4_BUILD_LOG.txt`
* `production_recovery/R4_RUN_LOG.txt`
* `production_recovery/R4_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R4_PASS_FAIL.md`
* `production_recovery/R4_evidence/README.md`
* `production_recovery/PRODUCTION_RECOVERY_R3_CLOSED.md`

## Modified files

* `src/CMakeLists.txt`
* `production_recovery/R3_BUILD_LOG.txt`
* `production_recovery/R3_RUN_LOG.txt`
* `production_recovery/R3_PASS_FAIL.md`
* `production_recovery/R3_SOURCE_DIFF_SUMMARY.md`

## `src/xcompact3d.f90` modified?

No.

## `src/CMakeLists.txt` modified?

Yes. Only the standalone `fibre_prod_ibm_interpolation_check` target was added for R4.

## Connected to `xcompact3d` executable?

No. `fibre_prod_ibm_delta.f90` and `fibre_prod_ibm_interpolation.f90` were not added to the `xcompact3d` executable source list and are not connected to the main time loop.

## Spreading/RHS/structure/contact/collision implementation?

No. R4 does not implement IBM spreading, RHS coupling, structure advancement, two-way FSI, wall contact, or fibre-fibre collision.
