# R3 Pass/Fail Record

Result: BLOCKED

## R3 PASS boundary

R3 PASS only means the grid adapter independent validation passed.

R3 PASS does not indicate IBM interpolation PASS.

R3 PASS does not indicate IBM spreading PASS.

R3 PASS does not indicate RHS coupling PASS.

R3 PASS does not indicate structure advancement PASS.

R3 PASS does not indicate wall contact PASS.

R3 PASS does not indicate fibre-fibre collision PASS.

R3 PASS does not indicate production DNS-FSI closure.

## Current result rationale

R3 is BLOCKED in this environment because CMake cannot configure the project without an available Fortran compiler, so the standalone `fibre_prod_grid_adapter_check` target could not be built or run here.
