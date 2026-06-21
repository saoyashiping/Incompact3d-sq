# R11 MPI Consistency Audit

Result: BLOCKED

## Reason

R11 build/run prerequisites were unavailable or incomplete; configure_status=1, hook_build_status=2, hook_run_status=127, xcompact_build_status=2.

## np=1/2/4 lambda=0 result

BLOCKED for np=1, np=2, and np=4.

## np=1/2/4 small-lambda result

BLOCKED for np=1, np=2, and np=4.

## Hook diagnostics

No real hook diagnostics were produced because the build/run workflow did not complete.

## NaN/Inf status

No runtime NaN/Inf evidence exists because no xcompact3d MPI run completed.

## R12 status

R12 was not entered.

## Evidence boundary

R11 evidence is limited to np=1/2/4 consistency for the controlled R10 hook.  BLOCKED does not prove PASS or FAIL and does not certify production DNS-FSI final closure.
