# Production Recovery R12 Source Diff Summary

## Modified files

- `production_recovery/R12_evidence/R12_VALIDATION_COMMAND_FIXED.sh`
- `production_recovery/R12_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R12_PASS_FAIL.md`
- `production_recovery/R12_evidence/README.md`
- `production_recovery/R12_evidence/R12_FIX_NOTE.md`

## Source boundary

- `src/*.f90` modified: no
- `src/xcompact3d.f90` modified: no
- `src/CMakeLists.txt` modified: no
- RK3 modified: no
- pressure/projection modified: no
- channel forcing modified: no
- restart/statistics/visualization modified: no
- closed stage directories modified: no
- R13 entered: no

## Fix summary

The R12 failure was caused by the validation script, not by the Fortran production source. The script used `rg`, which was unavailable in the user's environment, and wrote fixed FAIL audits instead of computing results from actual helper, R11, no-fibre, restart/statistics/visualization evidence.

The revised script uses `grep`, computes all audit booleans from actual logs, runs no-fibre np=1/2/4 baselines, and writes the final PASS/FAIL from the validation matrix.
