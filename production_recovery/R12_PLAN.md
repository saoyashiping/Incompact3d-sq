# Production Recovery R12 Plan

## Scope

R12 is the final production-recovery validation matrix for the controlled Xcompact3D fibre-FSI recovery branch.

R12 validates, in one reproducible workflow:

1. R11 controlled main-loop hook np=1/2/4 MPI consistency is re-run and must PASS;
2. R2-R10 standalone helper/check targets build and run;
3. no-fibre main-program baseline runs for np=1/2/4;
4. lambda=0 no-contamination evidence from R11 is PASS for np=1/2/4;
5. small-lambda response evidence from R11 is PASS for np=1/2/4;
6. restart/statistics/visualization smoke evidence exists in all main-program runs;
7. final closure boundaries are recorded.

## Non-goals

R12 does not introduce new physics.
R12 does not modify `src/xcompact3d.f90`.
R12 does not modify RK3, pressure/projection, channel forcing, restart, statistics, or visualization implementation.
R12 does not perform mesh refinement.
R12 does not run a long production DNS.
R12 does not claim high-Reynolds-number physical validation against experiments.

## Files created

- `production_recovery/R12_PLAN.md`
- `production_recovery/R12_PASS_FAIL.md`
- `production_recovery/R12_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R12_evidence/README.md`
- `production_recovery/R12_evidence/R12_VALIDATION_COMMAND_FIXED.sh`
- `production_recovery/R12_evidence/R12_FINAL_VALIDATION_MATRIX.md`
- `production_recovery/R12_evidence/R12_STANDALONE_HELPER_AUDIT.md`
- `production_recovery/R12_evidence/R12_NO_FIBRE_BASELINE_AUDIT.md`
- `production_recovery/R12_evidence/R12_RESTART_STATS_VISU_AUDIT.md`
- `production_recovery/R12_evidence/R12_FINAL_CLOSURE_REPORT.md`

## Pass/fail criteria

R12 PASS requires all of the following:

1. R11 rerun writes `production_recovery/R11_PASS_FAIL.md` with `Result: PASS`;
2. R11 MPI audit writes `Result: PASS`;
3. all six R11 lambda=0/small-lambda audits write `Result: PASS`;
4. R2-R10 standalone check targets all print their exact PASS strings;
5. no-fibre np=1/2/4 main-program runs all finish with `Good job! Xcompact3d finished successfully!`;
6. restart, snapshot/visualization, and statistics messages are found in main-program logs;
7. no R12 source modifications are required.

R12 FAIL applies if the workflow runs but any required validation criterion fails.
R12 BLOCKED applies if build, MPI, input case, or executable prerequisites are unavailable.
