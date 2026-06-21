# Production Recovery R12 Plan — Final Validation Matrix and Closure

## R12 scope

R12 performs the final Production Recovery validation matrix.  It reruns/reaudits the helper chain, reruns R11 consistency, checks no-fibre baseline np=1/2/4, checks restart/statistics/visualization smoke evidence, audits source boundaries, and writes final matrix and closure reports.

## R12 non-goals

R12 does not add physical models, does not modify source code, does not modify `src/xcompact3d.f90`, does not modify CMake, does not modify RK3, does not modify pressure/projection, does not modify channel forcing, does not modify restart/statistics/visualization semantics, does not enter R13, and does not claim long-duration production DNS or paper-grade physical validation.

## Validation matrix

The final matrix records: source boundary, standalone helper chain, R11 rerun, no-fibre np=1/2/4, lambda=0 np=1/2/4, small-lambda np=1/2/4, restart/statistics/visualization smoke, no source modification, no RK3/projection/channel forcing modification, and no R13 entry.

## Build strategy

Use `build_r12_final_validation`, configure with `mpif90` and the configured 2decomp root, then build the helper targets and `xcompact3d`.

## Run strategy

Run the R2-R10 standalone helper checks, rerun R11 validation, run no-fibre baseline np=1/2/4 with copied input files, and inspect restart/statistics/visualization smoke evidence.

## Pass/fail criteria

R12 PASS requires all matrix booleans to be `1` and both `R12_FINAL_VALIDATION_MATRIX.md` and `R12_FINAL_CLOSURE_REPORT.md` to contain `Result: PASS`.

R12 FAIL applies if any required validation matrix item fails or remains unavailable.

## Evidence boundary

R12 PASS would close the R0-R12 Production Recovery validation matrix only.  It would not claim long-time DNS completion, paper-grade physical statistics, all research cases, grid/time-step independence, or experimental comparison completion.
