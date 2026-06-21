# P0.9 validation environment fix

This patch fixes the latest P0.9 validation failure after executable-path discovery was repaired.

Observed failure:

```text
Running fibre_prod_force_buffer_rhs_path_check (.../build_p0_9/bin/fibre_prod_force_buffer_rhs_path_check)
ERROR STOP 9
... fibre_prod_force_buffer_rhs_path_check.f90:73
FAIL: run fibre_prod_force_buffer_rhs_path_check
```

Root cause:

`fibre_prod_force_buffer_rhs_path_check` is the P0.2 regression executable and intentionally requires the coupling environment variables `FIBRE_PROD_LAMBDA` and `FIBRE_PROD_PENALTY_BETA`. P0.9 ran this executable without those variables, so its local environment read failed at line 73.

Fix:

`P0_9_VALIDATION_COMMAND.sh` now runs only `fibre_prod_force_buffer_rhs_path_check` with the same deterministic P0.2 micro-coupling environment used by earlier passing P0.2/P0.3 validations:

```bash
FIBRE_PROD_ENABLE=1
FIBRE_PROD_LAMBDA=1.0e-3
FIBRE_PROD_PENALTY_BETA=2.0
FIBRE_PROD_DIAGNOSTICS=0
```

No Fortran production source was changed.
