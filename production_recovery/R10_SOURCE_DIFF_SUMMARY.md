# Production Recovery R10 Source Diff Summary

## New source files

- `src/fibre_prod_runtime_config.f90` — environment-gated runtime config with default disabled state and fail-closed validation.
- `src/fibre_prod_rhs_adapter.f90` — minimal gated RHS adapter with before/after signatures and no-write behavior for disabled or lambda=0 states.
- `src/fibre_prod_main_diagnostics.f90` — hook call, signature, finite, modified-cell, no-contamination, and small-lambda response diagnostics.
- `src/fibre_prod_main_hook.f90` — init/apply/finalize API used by the controlled main-loop hook.
- `src/fibre_prod_main_hook_check.f90` — standalone hook check for default disabled, lambda=0 no-op, small-lambda response, negative-lambda failure, NaN failure, and diagnostics finite behavior.

## New evidence/documentation files

- `production_recovery/R10_PLAN.md`
- `production_recovery/R10_BUILD_LOG.txt`
- `production_recovery/R10_RUN_LOG_hook_check.txt`
- `production_recovery/R10_RUN_LOG_lambda0_np1.txt`
- `production_recovery/R10_RUN_LOG_smalllambda_np1.txt`
- `production_recovery/R10_MAIN_HOOK_AUDIT.md`
- `production_recovery/R10_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R10_PASS_FAIL.md`
- `production_recovery/R10_evidence/README.md`
- `production_recovery/R10_evidence/R10_HOOK_SITE_AUDIT.txt`
- `production_recovery/R10_evidence/R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt`
- `production_recovery/R10_evidence/R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt`

## Modified files

- `src/xcompact3d.f90` — adds minimal R10 use/init/apply/finalize calls.
- `src/CMakeLists.txt` — adds R10 modules to `xcompact3d` and adds standalone `fibre_prod_main_hook_check` target.

## `src/xcompact3d.f90` modified?

Yes.  R10 is the first controlled main-loop hook stage and modifies `src/xcompact3d.f90` only for the R10 hook use/init/apply/finalize calls.

## `src/CMakeLists.txt` modified?

Yes.  R10 modules are added to the `xcompact3d` executable compile chain, and `fibre_prod_main_hook_check` is added as a standalone target.

## Connected to `xcompact3d` executable?

Yes.  R10 runtime config, diagnostics, RHS adapter, and main hook modules are added to the `xcompact3d` executable compile chain.

## Real RHS coupling implemented?

R10 implements only a strictly gated, diagnostic RHS adapter contribution at the audited hook site.  It is default-off, lambda=0 no-op, and small-lambda bounded.  It is not final production DNS-FSI closure.

## RK3 modified?

No.  R10 does not modify RK3 coefficients, `iadvance_time`, or `int_time(...)`.

## Pressure/projection modified?

No.  R10 does not modify `pre_correc(...)`, `calc_divu_constraint(...)`, `solve_poisson(...)`, or `cor_vel(...)`.

## Restart/statistics/visualization modified?

No.  R10 does not modify restart, statistics, or visualization semantics.

## np=1/2/4 final consistency executed?

No.  R10 does not execute R11 np=1/2/4 final consistency.

## production DNS-FSI final closure executed?

No.  R10 does not execute or claim production DNS-FSI final closure.

## R10 build/run final result

R10 is BLOCKED in this environment because `mpif90` is unavailable, so configure/build did not produce `fibre_prod_main_hook_check` or `xcompact3d`, and lambda=0/small-lambda np=1 runs could not be executed.
