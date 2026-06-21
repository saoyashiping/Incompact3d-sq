# Production Recovery R11 Plan — MPI Consistency Under Controlled R10 Hook

## R11 scope

R11 validates the controlled R10 main-loop hook under real MPI np=1/2/4 execution.  It checks the standalone hook target, `xcompact3d` build, lambda=0 no-contamination runs, small-lambda response runs, hook diagnostics, finite status, `last_status=0`, `failed_calls=0`, and normal program termination.

## R11 non-goals

R11 does not add or change physical models, does not change R10 hook mathematics, does not modify RK3, does not modify pressure/projection, does not modify channel forcing, does not modify restart/statistics/visualization semantics, does not enter R12, and does not claim production DNS-FSI final closure.

## Validation strategy

R11 uses `production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh` from the repository root.  The script configures a clean build, builds and runs `fibre_prod_main_hook_check`, builds `xcompact3d`, prepares isolated run directories with copied input files, and runs lambda=0 and small-lambda cases for np=1/2/4 using `mpirun --oversubscribe`.

## np=1/2/4 commands

Lambda=0 command form:

```bash
FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=0 FIBRE_PROD_DIAGNOSTICS=1 \
  FIBRE_PROD_DIAGNOSTICS_DIR=<absolute R11_evidence path> \
  mpirun --oversubscribe -np <1|2|4> <xcompact3d executable>
```

Small-lambda command form:

```bash
FIBRE_PROD_ENABLE=1 FIBRE_PROD_LAMBDA=1.0e-8 FIBRE_PROD_DIAGNOSTICS=1 \
  FIBRE_PROD_DIAGNOSTICS_DIR=<absolute R11_evidence path> \
  mpirun --oversubscribe -np <1|2|4> <xcompact3d executable>
```

## Pass/fail criteria

R11 PASS requires the hook check PASS, xcompact3d build success, six main-program logs ending successfully, three lambda=0 audit PASS files, three small-lambda audit PASS files, diagnostics with finite signatures, lambda=0 `modified_cells=0`, small-lambda `modified_cells>0`, `last_status=0`, `failed_calls=0`, and `R11_MPI_CONSISTENCY_AUDIT.md` PASS.

R11 FAIL applies if build/run completes but any required consistency, finite, no-contamination, response, or status condition fails.

R11 BLOCKED applies if configure/build/run cannot complete because the required compiler, MPI wrapper, runtime launcher, dependency path, input case, or executable is unavailable.

## Evidence boundary

R11 PASS would mean the R10 controlled hook passed np=1/2/4 main-program consistency smoke validation.  R11 PASS would not mean R12 paper-level validation matrix PASS and would not mean production DNS-FSI final closure.
