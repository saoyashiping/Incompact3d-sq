# P1_1 source diff summary

Modified files:

- `production_recovery/P1_1_evidence/P1_1_VALIDATION_COMMAND.sh`
- `production_recovery/P1_1_evidence/P1_1_FIX_NOTE.md`
- `production_recovery/run_P1_1.sh`

No Fortran production source was modified.

The validation script now:

- auto-refreshes P1_0 if P1_0 PASS evidence was overwritten by a fresh unzip;
- avoids false failure from whole-source RHS API name matches;
- keeps P1_1 as a real xcompact3d channel DNS run on `96x97x96`, `dt=5e-5`, `ilast=100`, `1 fibre`, `49 nodes`, `lambda=0`;
- writes explicit PASS/FAIL reasons.
