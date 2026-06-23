# P1_3 source/script diff summary

Modified files:
- production_recovery/P1_3_evidence/P1_3_VALIDATION_COMMAND.sh
- production_recovery/run_P1_3.sh
- production_recovery/P1_2_evidence/P1_2_VALIDATION_COMMAND.sh
- production_recovery/P1_3_evidence/P1_3_PENDING_FIX_NOTE.md
- production_recovery/P1_3_SOURCE_DIFF_SUMMARY.md

Scope:
- Validation-script robustness only.
- No change to xcompact3d pressure/projection/RK3/channel-forcing logic.
- No change to fibre force, structure, IBM spreading, or RHS-coupling Fortran source.
- P1_3 remains a real 128x129x96 channel DNS-FSI guarded stability/restart/statistics/visualization validation.
