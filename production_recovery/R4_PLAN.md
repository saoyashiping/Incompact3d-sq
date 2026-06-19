# Production Recovery R4 Plan

## R4 scope

R4 implements a standalone production IBM interpolation layer using the R3 grid adapter. It provides compact-support regularized delta weights plus scalar and velocity interpolation for production fibre-system use.

## R4 non-goals

R4 does not connect to `xcompact3d.f90`, does not enter the production time loop, does not implement IBM spreading, does not inject RHS forcing, does not advance structure, does not implement two-way FSI, does not implement wall contact, and does not implement fibre-fibre collision.

## Files created/modified

Created:

* `src/fibre_prod_ibm_delta.f90`
* `src/fibre_prod_ibm_interpolation.f90`
* `src/fibre_prod_ibm_interpolation_check.f90`
* `production_recovery/R4_PLAN.md`
* `production_recovery/R4_BUILD_LOG.txt`
* `production_recovery/R4_RUN_LOG.txt`
* `production_recovery/R4_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R4_PASS_FAIL.md`
* `production_recovery/R4_evidence/README.md`

Modified:

* `src/CMakeLists.txt`
* R3 evidence documents listed in `R4_SOURCE_DIFF_SUMMARY.md`

## Build strategy

Build only the standalone `fibre_prod_ibm_interpolation_check` target from the grid adapter, IBM delta, IBM interpolation, and independent check driver sources.

## Run strategy

Run only the standalone `fibre_prod_ibm_interpolation_check` executable. Success requires real output containing `R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK PASS`.

## Pass/fail criteria

PASS requires the standalone target to compile, run, and print `R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK PASS` while checking compact-support delta behavior, constant scalar/velocity preservation, affine scalar sanity, periodic x/z wrapping, y fail-closed behavior, positive weight sum, and destroy/deallocation.

FAIL applies if the standalone target builds and runs but any required interpolation check fails or R4 violates forbidden-scope boundaries.

BLOCKED applies if external environment prerequisites prevent building or running the standalone target.

## Evidence boundary

R4 PASS only means standalone IBM interpolation validation passed. It does not mean IBM spreading, RHS coupling, structure advancement, two-way FSI, wall contact, fibre-fibre collision, MPI, DNS, or production DNS-FSI closure passed.
