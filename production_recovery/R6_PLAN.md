# Production Recovery R6 Plan

## R6 scope

R6 adds standalone production IBM spreading from Lagrangian fibre point forces to an Eulerian force-density buffer using the R3 grid adapter and R4 compact-support regularized delta kernel.

## R6 non-goals

R6 does not connect to `xcompact3d.f90`, does not write the real Navier-Stokes RHS, does not modify time advancement, does not implement a two-way FSI closed loop, does not advance structure, does not implement wall contact, does not implement fibre-fibre collision, and does not perform production DNS-FSI validation.

## Files created/modified

Created:

* `src/fibre_prod_ibm_force_buffer.f90`
* `src/fibre_prod_ibm_spreading.f90`
* `src/fibre_prod_ibm_spreading_check.f90`
* `production_recovery/R6_PLAN.md`
* `production_recovery/R6_BUILD_LOG.txt`
* `production_recovery/R6_RUN_LOG.txt`
* `production_recovery/R6_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R6_PASS_FAIL.md`
* `production_recovery/R6_evidence/README.md`

Modified:

* `src/CMakeLists.txt`

## Build strategy

Build only the standalone `fibre_prod_ibm_spreading_check` target from the grid adapter, delta kernel, interpolation module, force buffer, spreading module, and independent check driver.

## Run strategy

Run only the standalone `fibre_prod_ibm_spreading_check` executable. Success requires real output containing `R6_FIBRE_PROD_IBM_SPREADING_CHECK PASS`.

## Pass/fail criteria

PASS requires the standalone target to compile, run, and print `R6_FIBRE_PROD_IBM_SPREADING_CHECK PASS` while checking buffer allocation/reset/destroy, single-point and multi-point force conservation, x/z periodic wrapping, y fail-closed behavior, NaN/Inf fail-closed behavior, positive weight sum, finite buffer, zero integral after reset, and constant-velocity power consistency.

FAIL applies if the standalone target builds and runs but any required spreading check fails or R6 violates forbidden-scope boundaries.

BLOCKED applies if external environment prerequisites prevent building or running the standalone target.

## Evidence boundary

R6 PASS only means production IBM spreading independent validation passed. It does not mean RHS coupling, structure advancement, two-way FSI, wall contact, fibre-fibre collision, MPI, DNS, or production DNS-FSI closure passed.
