# R11 Fix Note

## Diagnosis

The latest R11 terminal output shows that the R11 build and MPI runs reached real execution, including successful Xcompact3D completions under MPI. However, the R11 evidence files remained at the old `Result: BLOCKED` state.

The failure was caused by two coupled issues:

1. `production_recovery/R11_evidence/R11_VALIDATION_COMMAND_FIXED.sh` executed the six MPI runs but did not perform post-processing after the run loop. It did not copy the R10 runtime audit/diagnostic files into R11-named evidence files, did not evaluate the six np=1/2/4 criteria, did not rewrite `R11_MPI_CONSISTENCY_AUDIT.md`, and did not rewrite `R11_PASS_FAIL.md`.
2. The latest source package had regressed to an older R10 hook diagnostic implementation. It did not expose `last_status` and `failed_calls`, and did not generate exact `R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt` / `R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt` files from `fibre_prod_main_hook_finalize`. R11 therefore could not satisfy its own diagnostics criteria.

## Fix

This patch restores the stricter R10 runtime diagnostics and replaces the R11 validation script with a complete build-run-postprocess-audit workflow.

The fixed R11 script now:

1. deletes stale R10/R11 audit files before each validation attempt;
2. builds `fibre_prod_main_hook_check` and `xcompact3d`;
3. runs lambda=0 and small-lambda cases for np=1, np=2, and np=4;
4. uses a per-run `FIBRE_PROD_DIAGNOSTICS_DIR` so each MPI run writes its own R10 audit/diagnostic files;
5. copies those generated R10 files into R11-named evidence files;
6. checks exact `Result: PASS` audit lines;
7. checks diagnostics fields including `modified_cells`, `last_status`, `failed_calls`, finite flags, `no_contamination`, and `small_lambda_response`;
8. rewrites `R11_MPI_CONSISTENCY_AUDIT.md` and `R11_PASS_FAIL.md` from real evidence.

## Boundary

This patch does not enter R12.
It does not change the intended R10 hook location in `xcompact3d.f90`; it only restores the safer finalize behavior so diagnostics are always written after an initialized hook.
It does not modify RK3 coefficients, pressure/projection, channel forcing, restart/statistics/visualization semantics, or any closed Stage 20/21/22 files.
