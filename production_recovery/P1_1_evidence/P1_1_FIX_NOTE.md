# P1_1 validation script fix

This fix makes P1_1 validation self-contained and removes two script-level false-fail mechanisms.

## Fixed issues

1. P1_1 previously required pre-existing P1_0 PASS evidence. When a fresh source zip overwrote evidence files back to `PENDING`, P1_1 failed before the real 96x97x96 run. The script now automatically refreshes P1_0 if the required P1_0 PASS evidence is missing or overwritten.

2. The static RHS feedback audit searched the whole `src/` tree for names such as `fibre_prod_force_buffer_rhs_gate_apply`. That incorrectly matched legitimate definitions and non-P1 modules. The audit is now narrowed to the P1_1 one-way module and exact old unsafe uniform RHS formula.

3. The NaN/Inf audit now avoids false positives from safety text such as `No NaN/Inf detected`.

## Scope

No Fortran production physics source is changed. The fix changes validation/runner scripts only.
