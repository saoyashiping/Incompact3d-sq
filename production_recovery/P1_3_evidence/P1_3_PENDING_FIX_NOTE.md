# P1_3 PENDING-result validation-script fix

Root cause:
- P1_3_VALIDATION_COMMAND.sh used `set -euo pipefail` but only explicit `fail()` calls rewrote P1_3_PASS_FAIL.md.
- CMake/build/dependency/P1_2 rerun/xcompact3d/grep/python failures could terminate the shell before `fail()` was called, leaving the original `Result: PENDING` files untouched.
- The old unsafe-uniform-RHS audit searched the plain text phrase `uniform RHS contribution` across all src files. This can match harmless diagnostic strings such as `no uniform RHS contribution audit`, producing a false failure.
- The old NaN/Inf audit could also match safe strings such as `no NaN/Inf detected`.

Fix:
- Add ERR trap and RUNNING state so every early exit becomes explicit Result: FAIL with line/command context.
- Add P1_3_FAILURE_CONTEXT.txt for debugging unhandled exits.
- Replace broad uniform-RHS text search with targeted unsafe formula-pattern search.
- Make NaN/Inf audit ignore safe no-NaN diagnostic lines.
- Make P1_2 dependency handling explicit: use existing PASS evidence, accept valid lambda logs, or rerun P1_2; never leave P1_3 as PENDING.
- Update run_P1_3.sh to call the fixed validation script through bash.

No Fortran physics source is changed by this patch.
