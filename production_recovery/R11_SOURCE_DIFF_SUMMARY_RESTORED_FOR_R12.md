# Production Recovery R11 Source Diff Summary

## Files modified by the R11 fix

- `src/fibre_prod_runtime_config.f90`
- `src/fibre_prod_main_diagnostics.f90`
- `src/fibre_prod_main_hook.f90`
- `src/xcompact3d.f90`
- `production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh`
- `production_recovery/R11_evidence/R11_FIX_NOTE.md`
- `production_recovery/R11_SOURCE_DIFF_SUMMARY.md`

## Why source files were modified

R11 originally should have been script-only. However, the latest uploaded source had regressed to an older R10 diagnostic implementation that could not produce the fields required by the R11 criteria:

- `last_status`
- `failed_calls`
- exact runtime audit files for lambda=0 and small lambda
- mode-specific diagnostics files

Therefore, the minimum R10 diagnostic/hook source restoration was necessary before R11 could be evaluated honestly.

## `src/xcompact3d.f90`

Modified only to restore the safer R10 finalize behavior:

- keep `fibre_prod_r10_hook_initialized` separate from `fibre_prod_r10_hook_ready`;
- allow `fibre_prod_main_hook_apply` to be disabled after a nonzero status;
- still call `fibre_prod_main_hook_finalize` if the hook was initialized, so audit files are written even after an apply-side failure.

No RK3 coefficient was changed.
No pressure/projection logic was changed.
No channel forcing logic was changed.
No restart/statistics/visualization semantics were changed.

## R11 validation script

`R11_VALIDATION_COMMAND_FIXED.sh` was replaced with a complete validation workflow:

1. clean stale evidence;
2. configure/build;
3. run hook check;
4. build xcompact3d;
5. run lambda=0 np=1/2/4;
6. run small-lambda np=1/2/4;
7. copy per-run R10 audit/diagnostic files into R11 evidence names;
8. check exact pass/fail criteria;
9. write `R11_MPI_CONSISTENCY_AUDIT.md`;
10. write `R11_PASS_FAIL.md`.

## Closed-stage boundary

The following directories were not modified:

- `stage20_checks/`
- `stage21_checks/`
- `stage22_checks/`

## R12 boundary

R12 was not entered.

## Meaning of a future R11 PASS

R11 PASS means only that the R10 controlled main-loop hook passes np=1/2/4 MPI smoke consistency checks for lambda=0 and small lambda.

R11 PASS does not mean R12 paper-level validation matrix PASS.
R11 PASS does not mean production DNS-FSI final closure.
