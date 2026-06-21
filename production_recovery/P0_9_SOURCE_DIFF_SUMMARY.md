# P0.9 source diff summary

## Files changed

- `production_recovery/P0_9_evidence/P0_9_VALIDATION_COMMAND.sh`
- `production_recovery/P0_9_evidence/P0_9_FIX_NOTE.md`
- `production_recovery/P0_9_SOURCE_DIFF_SUMMARY.md`

## Summary

This repair is validation-script-only. It keeps the previous executable path discovery fix and adds target-specific runtime environment handling for the P0.2 regression executable `fibre_prod_force_buffer_rhs_path_check`.

The executable itself requires `FIBRE_PROD_LAMBDA` and `FIBRE_PROD_PENALTY_BETA`; running it without these variables causes `ERROR STOP 9`. The P0.9 script now invokes this check with deterministic small-lambda settings while leaving all other checks unchanged.

## Production impact

No Fortran production source was changed. No RHS/structure/IBM behavior was modified. P0.9 remains a guarded validation stage and does not imply production DNS-FSI readiness.
