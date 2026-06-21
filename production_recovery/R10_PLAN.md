# Production Recovery R10 Plan — Controlled Main-Loop Hook

## R10 scope

R10 adds the first controlled `xcompact3d` main-loop production fibre hook.  The hook is default-off, environment-gated, lambda=0 no-op, and small-lambda diagnostic-only.  It applies a bounded RHS contribution only through the R10 adapter and records signatures for no-contamination and response audits.

## R10 non-goals

R10 does not complete production DNS-FSI closure, does not modify RK3 coefficients, does not modify pressure/projection mathematics, does not modify channel forcing logic, does not modify restart/statistics/visualization semantics, does not run R11 np=1/2/4 final consistency, and does not enter R12.

## Files created/modified

Created source files:

- `src/fibre_prod_runtime_config.f90`
- `src/fibre_prod_rhs_adapter.f90`
- `src/fibre_prod_main_diagnostics.f90`
- `src/fibre_prod_main_hook.f90`
- `src/fibre_prod_main_hook_check.f90`

Created evidence files:

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

Modified source/build files:

- `src/xcompact3d.f90`
- `src/CMakeLists.txt`

## xcompact3d hook strategy

The hook is inserted after RHS construction and existing Stage 14 RHS injection, before `int_time(...)`.  This is a controlled RHS hook site and does not bypass pressure/projection or RK3.

## Build strategy

Configure a dedicated `build_r10_main_hook` directory, build the standalone `fibre_prod_main_hook_check` target, then build `xcompact3d`.

## Run strategy

Run `fibre_prod_main_hook_check`, then run `xcompact3d` np=1 twice: once with `FIBRE_PROD_LAMBDA=0` and once with `FIBRE_PROD_LAMBDA=1.0e-8`.  Do not modify user input files; copy any required input into R10 evidence run directories.

## Pass/fail criteria

R10 PASS requires the standalone hook check PASS, `xcompact3d` build success, lambda=0 np=1 run success with no-contamination audit PASS, small-lambda np=1 run success with finite response audit PASS, no NaN/Inf, and no R11/R12 claims.

R10 FAIL applies if build/run completes but required hook safety checks fail.

R10 BLOCKED applies if configure/build/run cannot complete due to missing compiler, MPI wrapper, dependency, input case, or executable.

## Evidence boundary

R10 PASS would mean controlled main-loop hook initial integration only.  It would not mean np=1/2/4 MPI production consistency, restart/statistics/visualization final validation, paper-level validation matrix, or production DNS-FSI final closure.
