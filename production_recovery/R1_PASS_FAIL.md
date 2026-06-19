# R1 Pass/Fail Record

## R1 PASS conditions

R1 PASS requires successful CMake configure, successful build of the baseline production executable, and real baseline no-fibre execution for `np=1`, `np=2`, and `np=4` with no NaN/Inf or abnormal termination.

## R1 FAIL conditions

R1 FAIL applies when required tooling and dependencies are present but the baseline production path cannot configure, build, or run due to project errors that are not corrected within R1's allowed minimal compile-fix scope.

## R1 BLOCKED conditions

R1 BLOCKED applies when missing compiler, MPI tools, 2decomp-fft/decomp2d, input prerequisites, or other external environment limitations prevent real build/run evidence from being generated.

## Current conclusion

**BLOCKED.**

The build did not succeed; therefore R1 cannot be PASS. The `np=1`, `np=2`, and `np=4` runs were not executed; therefore R1 cannot be a complete PASS.

## Boundary statement

R1 PASS does not equal FSI solver PASS.

R1 PASS only indicates that the baseline build/run path has the foundation needed to enter R2.
