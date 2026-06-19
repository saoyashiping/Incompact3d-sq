# R1 Blocked Record

## Current conclusion

R1 is BLOCKED.

## Blocking reasons

* `cmake` is present, but CMake configure failed because no Fortran compiler was found.
* `gfortran` is not available in `PATH`.
* `mpif90` is not available in `PATH`.
* `mpirun` is not available in `PATH`.
* No usable 2decomp-fft/decomp2d installation was discovered by the environment probe.

## Evidence

The raw configure error is preserved in `R1_BUILD_LOG.txt`. Environment probes are preserved in `R1_ENVIRONMENT.md`. The `np=1`, `np=2`, and `np=4` run logs each contain an explicit BLOCKED reason instead of fabricated runtime output.

## Required unblockers

Install or expose a Fortran compiler, MPI Fortran wrapper, MPI launcher, and the required 2decomp-fft/decomp2d dependency. Then rerun R1 from a clean build directory and capture real configure, build, and baseline no-fibre run logs.
