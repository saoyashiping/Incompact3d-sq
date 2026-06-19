# R4 Pass/Fail Record

Result: BLOCKED

## R4 PASS condition

R4 PASS requires the standalone `fibre_prod_ibm_interpolation_check` target to build and run, with `R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK PASS` present in `R4_RUN_LOG.txt`.

## Current result rationale

R4 is BLOCKED in this environment because the requested configure command could not find `mpif90`, the build command had no generated Makefile, and the standalone executable did not exist for the run command.

## R4 PASS boundary

R4 PASS only means standalone IBM interpolation validation passed.

R4 PASS does not indicate IBM spreading PASS.

R4 PASS does not indicate RHS coupling PASS.

R4 PASS does not indicate structure advancement PASS.

R4 PASS does not indicate two-way FSI PASS.

R4 PASS does not indicate wall contact PASS.

R4 PASS does not indicate fibre-fibre collision PASS.

R4 PASS does not indicate production DNS-FSI closure.
