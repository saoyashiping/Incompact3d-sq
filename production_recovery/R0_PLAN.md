# Production Recovery R0 Plan

## R0 scope

R0 corrects project status language and resets the evidence boundary. It records that Stage 22 source-only closure exists while production DNS-FSI closure is not complete.

## R0 non-goals

R0 does not implement solver functionality, production coupling, physics validation, restart integration, statistics output, visualization output, MPI execution, or build/run validation.

## Files to be created/modified

Created files:

* `PRODUCTION_RECOVERY_STATUS.md`
* `production_recovery/R0_PLAN.md`
* `production_recovery/R0_SOURCE_ONLY_AUDIT.md`
* `production_recovery/R0_PASS_FAIL.md`
* `production_recovery/README.md`

Modified historical files:

* `STAGE22_CLOSED.md`
* `PROJECT_FINAL_CLOSED.md`

## Pass/fail criteria

R0 passes only if the authoritative status is corrected, required R0 files exist, historical closure files retain their original content with a top warning block, and no forbidden production/build/DNS/MPI actions are performed.

R0 fails if it modifies solver source, changes production build targets, changes production inputs, deletes Stage 10-22 evidence, runs build/DNS/MPI validation, creates synthetic PASS evidence as a substitute for real execution, or leaves unconditional production-ready status as authoritative.

## Forbidden operations

* Modify `src/*.f90`.
* Modify `src/CMakeLists.txt`.
* Modify `input/*.i3d` or production run input files.
* Run `cmake`, `make`, `ninja`, `mpirun`, `xcompact3d`, or DNS/MPI/build validation.
* Generate new synthetic PASS evidence to replace real execution.
* Delete existing Stage 10-22 evidence files.
* Rewrite Stage 22 history as if it never occurred.

## Evidence boundary

R0 establishes that Stage 22 closure files are historical source-only audit artifacts. They may support source-audit traceability, but they do not certify production DNS-FSI closure, real build/run validation, real MPI execution, or production contact/collision coupling.
