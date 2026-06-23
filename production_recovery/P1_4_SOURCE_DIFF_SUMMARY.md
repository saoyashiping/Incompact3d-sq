# P1_4 source/script diff summary

Modified files:
- src/fibre_prod_p1_np_consistency_closure_case.f90
- production_recovery/P1_4_evidence/P1_4_VALIDATION_COMMAND.sh
- production_recovery/run_P1_4.sh
- production_recovery/P1_4_evidence/P1_4_SELF_CONTAINED_SCRIPT_SOURCE_FIX_NOTE.md
- production_recovery/P1_4_SOURCE_DIFF_SUMMARY.md

Purpose:
- Fix P1_4 validation false failure from harmless diagnostic text containing `uniform RHS contribution`.
- Harden self-contained P1_4 validation so it reports explicit failure context and performs actual np=1/2/4 signature comparisons.
- Keep P1_4 independent from P1_0-P1_3 evidence files.
- Preserve production-run status as STILL BLOCKED FOR PAPER-SCALE DNS until later P stages.
