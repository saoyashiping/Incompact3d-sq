# R6 Source Diff Summary

## Added files

* `src/fibre_prod_ibm_force_buffer.f90`
* `src/fibre_prod_ibm_spreading.f90`
* `src/fibre_prod_ibm_spreading_check.f90`
* `production_recovery/R6_PLAN.md`
* `production_recovery/R6_BUILD_LOG.txt`
* `production_recovery/R6_RUN_LOG.txt`
* `production_recovery/R6_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R6_PASS_FAIL.md`
* `production_recovery/R6_evidence/README.md`

## Modified files

* `src/CMakeLists.txt`

## `src/xcompact3d.f90` modified?

No.

## `src/CMakeLists.txt` modified?

Yes. Only the standalone `fibre_prod_ibm_spreading_check` target was added for R6.

## Connected to `xcompact3d` executable?

No. The R6 modules were not added to the `xcompact3d` executable source list and are not connected to the main time loop.

## Real RHS coupling implemented?

No. R6 writes only to a standalone force-density buffer and does not write the real Navier-Stokes RHS.

## Two-way FSI/structure/contact/collision implementation?

No. R6 does not implement two-way FSI, structure advancement, wall contact, or fibre-fibre collision.
