# R12 Fix Note

## Problem diagnosed

The previous R12 validation script failed for script-level reasons rather than Fortran source-level reasons:

1. it used `rg` even though `ripgrep` was not installed on the target Ubuntu environment;
2. it wrote fixed FAIL/BLOCKED-style audit text even after helper targets and xcompact3d were successfully built and run;
3. `R12_STANDALONE_HELPER_AUDIT.md`, `R12_NO_FIBRE_BASELINE_AUDIT.md`, and the final matrix were not computed from the actual run logs;
4. no-fibre np=1/2/4 baselines were never actually launched by the script; placeholder FAIL logs were written instead.

## Fix applied

`production_recovery/R12_evidence/R12_VALIDATION_COMMAND_FIXED.sh` was rewritten to:

1. use POSIX `grep` instead of `rg`;
2. build and run R2-R10 helper targets and compute PASS/FAIL from real PASS tokens;
3. rerun R11 and compute PASS/FAIL from exact `Result: PASS` lines;
4. build `xcompact3d` and run no-fibre np=1/2/4 with `FIBRE_PROD_ENABLE=0`;
5. check restart, snapshot, and statistics smoke evidence from actual no-fibre logs;
6. compute `R12_FINAL_VALIDATION_MATRIX.md` from real audit files;
7. write `R12_PASS_FAIL.md` from the final matrix and closure report.

## Boundary

No `src/*.f90` file was modified.
No `src/xcompact3d.f90` change was made.
No `src/CMakeLists.txt` change was made.
No RK3, pressure/projection, channel forcing, restart, statistics, or visualization code was modified.
R13 was not entered.
