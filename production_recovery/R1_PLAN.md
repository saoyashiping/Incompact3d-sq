# Production Recovery R1 Plan

## R1 scope

R1 attempts to establish real baseline build/run evidence for the current Xcompact3D production path with fibre/FSI disabled. The goal is to prove, or honestly attempt to prove, that the no-fibre baseline can configure, compile, and run under real tooling.

## R1 non-goals

R1 does not validate FSI, does not implement new FSI physics, does not implement IBM interpolation/spreading, does not implement a structure solver, does not implement wall contact, and does not implement fibre-fibre collision.

## Allowed source changes

Only minimal compile fixes are allowed if they directly block baseline compilation. Any such change must be listed in `R1_SOURCE_DIFF_SUMMARY.md` with its reason and compile-fix classification.

## Forbidden changes

R1 must not add `fibre_prod_*` production modules, modify Stage 20-22 synthetic checks to create PASS evidence, remove R0 warnings, bypass core solver calls, bypass Poisson/projection, bypass momentum/RHS, or claim production DNS-FSI closure.

## Build strategy

Use an isolated build directory:

```sh
mkdir -p build_r1_baseline
cmake -S . -B build_r1_baseline
cmake --build build_r1_baseline -j 2
```

All configure and build output is saved in `R1_BUILD_LOG.txt`.

## Run strategy

If the build succeeds, select a small existing baseline test input, copy it into `production_recovery/R1_evidence/`, confirm fibre/FSI is disabled or absent, and run the built executable with `mpirun -np 1`, `mpirun -np 2`, and `mpirun -np 4`. Each run log must contain real command output or an explicit BLOCKED reason.

## Pass/fail criteria

R1 PASS requires successful configure, successful build, and real no-fibre baseline execution for `np=1`, `np=2`, and `np=4` without NaN/Inf or abnormal termination.

R1 FAIL applies when required tooling and inputs are available but baseline configure/build/run fails due to project errors not fixed within R1's allowed minimal compile-fix scope.

R1 BLOCKED applies when missing external tooling, compiler, MPI, 2decomp-fft/decomp2d, or runnable baseline prerequisites prevent real build/run evidence from being generated.

## Evidence files

* `production_recovery/R1_ENVIRONMENT.md`
* `production_recovery/R1_BUILD_LOG.txt`
* `production_recovery/R1_SOURCE_DIFF_SUMMARY.md`
* `production_recovery/R1_BASELINE_RUN_PLAN.md`
* `production_recovery/R1_RUN_LOG_np1.txt`
* `production_recovery/R1_RUN_LOG_np2.txt`
* `production_recovery/R1_RUN_LOG_np4.txt`
* `production_recovery/R1_BASELINE_NO_FIBRE_STATUS.md`
* `production_recovery/R1_PASS_FAIL.md`
* `production_recovery/R1_BLOCKED.md` when blocked
* `production_recovery/R1_evidence/`
