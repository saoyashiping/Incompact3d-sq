# R10 Fix Note

## Diagnosed problem

The R10 Fortran hook-check and `xcompact3d` np=1 runs can complete, but the R10 audit files may remain as stale `BLOCKED` files generated before the real run. The previous shell-side validation also used broad `grep -qi "PASS"`, which can falsely pass a `BLOCKED` audit because the text contains phrases such as `Required PASS evidence`.

## Fix applied

1. Added `FIBRE_PROD_DIAGNOSTICS_DIR` support in `fibre_prod_runtime_config.f90`.
2. Made `fibre_prod_main_hook_finalize` write runtime diagnostics and the correct R10 audit file to the configured diagnostics directory.
3. Added exact audit writers for:
   - `R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt`
   - `R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt`
4. Audit files now contain explicit `Result: PASS` or `Result: FAIL` lines based on runtime diagnostics.
5. Direct standalone gfortran validation confirms the R10 hook check still prints `R10_FIBRE_PROD_MAIN_HOOK_CHECK PASS`.
6. Direct runtime-hook validation confirms the lambda=0 and small-lambda audit writers produce `Result: PASS` under controlled artificial RHS buffers.

## Boundary

This fix does not enter R11.
This fix does not perform np=1/2/4 MPI consistency.
This fix does not claim final production DNS-FSI closure.
This fix does not change pressure projection, RK3 coefficients, channel forcing, restart, statistics, or visualization logic.
