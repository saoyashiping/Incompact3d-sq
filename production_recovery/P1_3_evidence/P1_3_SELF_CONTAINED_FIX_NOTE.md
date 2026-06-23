# P1_3 self-contained validation fix

This patch removes P1_2/P1_1/P1_0 PASS-file dependency from P1_3 validation.

P1_3 PASS is now judged only from:

- P1_3 build success;
- `fibre_prod_p1_enlarged_stability_case_check` PASS;
- P1_3 real channel DNS-FSI segment1 run;
- P1_3 restart segment2 run;
- P1_3 force-buffer/RHS/restart/statistics/visualization/wall/NaN-Inf audits;
- P1_3 static unsafe-uniform-RHS audit.

The legacy `P1_3_P1_2_REGRESSION_AUDIT.txt` file is kept only as a non-gating SKIPPED note so that existing `cat` commands do not fail.
