# P1_3 source diff summary — self-contained validation

Modified files:

- `production_recovery/P1_3_evidence/P1_3_VALIDATION_COMMAND.sh`
- `production_recovery/run_P1_3.sh`
- `production_recovery/P1_3_PASS_FAIL.md`
- `production_recovery/P1_3_evidence/P1_3_VALIDATION_RESULT.txt`
- `production_recovery/P1_3_evidence/P1_3_P1_2_REGRESSION_AUDIT.txt`
- `production_recovery/P1_3_evidence/P1_3_SELF_CONTAINED_FIX_NOTE.md`

Changes:

1. Removed the P1_2 PASS-file and P1_2 log dependency from P1_3 pass/fail logic.
2. P1_3 now validates itself only through its own 128x129x96 real DNS-FSI segment1/restart-segment2 evidence.
3. Added robust `ERR` trapping so unexpected failures write `Result: FAIL` instead of leaving `Result: PENDING`.
4. Kept `P1_3_P1_2_REGRESSION_AUDIT.txt` as a non-gating `SKIPPED` note only.
5. Did not modify Fortran physics source code.

This patch does not claim P1_3 PASS. The user must rerun the P1_3 validation script on Ubuntu.

## P1_3 RHS formula grep validation fix

- File changed: `production_recovery/P1_3_evidence/P1_3_VALIDATION_COMMAND.sh`.
- Added note: `production_recovery/P1_3_evidence/P1_3_RHS_FORMULA_GREP_FIX_NOTE.md`.
- Reason: P1_3 runtime logs emit `P1_3 RHS increment diagnostic: nonzero finite bounded PASS formula=lambda*penalty_beta*force_buffer values=...`, but the previous script searched for the nonexistent token `RHS increment formula: PASS`.
- Effect: P1_3 no longer fails a valid segment1 run with `segment 1 missing RHS formula PASS`.
- Source scope: validation script only; no DNS-FSI physics source is modified.
