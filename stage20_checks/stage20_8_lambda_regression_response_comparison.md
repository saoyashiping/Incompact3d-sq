# Stage 20.8: lambda=0 regression and small-lambda response comparison

Stage 20.8 is a diagnostic-only helper-local regression/comparison audit. It reconstructs the controlled one-fibre `np=1` closed-loop response for lambda-zero, lambda-small, and lambda-double cases to verify strict no-op behavior at `lambda=0`, bounded helper-local RHS response at small lambda, and deterministic scaling for a frozen force-density subtest.

This stage is not a production DNS run. It does not run MPI, production validation, production RHS updates, Stage 14 RHS injection, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, restart, statistics, or visualization production paths.

## Source-only and no-rerun policy

Stage 20.8 accepts Stage 20.7 PASS evidence when present and preserves source-only acceptance for prior Stage 20 closure behavior. It does not require old closure files, does not require all old stage outputs, and does not rerun Stage 10-19 or Stage 20.0-20.7.

## Comparison cases

* `lambda_zero = 0.0`
* `lambda_small = 1.0e-6` by default
* `lambda_double = 2.0e-6` by default when the double-lambda scaling check is enabled

At each helper-local step the audit computes the controlled fluid velocity, `F_fs_candidate`, hydrodynamic structure advance, helper-local X/V/A commit, reaction force, nearest-grid-point force-density candidate, and RHS candidate response.

```text
RHS_delta_zero = 0.0 * f_eulerian_candidate_zero
RHS_after_zero = RHS_before + RHS_delta_zero
RHS_zero_residual = RHS_after_zero - RHS_before

RHS_delta_small = lambda_small * f_eulerian_candidate_small
RHS_small_formula_residual = RHS_delta_small - lambda_small * f_eulerian_candidate_small

RHS_delta_double = lambda_double * f_eulerian_candidate_double
RHS_double_formula_residual = RHS_delta_double - lambda_double * f_eulerian_candidate_double
```

The frozen force-density scaling subtest uses identical helper-local force-density data for `lambda_small` and `lambda_double`, so `norm(RHS_delta_double) / norm(RHS_delta_small)` matches `lambda_double / lambda_small` within tolerance.

## Safety boundary

All X/V/A commits are helper-local diagnostic state only. Production two-way coupling, production structure-to-fluid coupling, production RHS coupling, production RHS updates, Stage 14 RHS injection, wall contact, fibre-fibre collision, and production multifibre logic remain disabled.

## Next stage

Stage 20.9: parallel consistency np=2/4 for two-way coupling.

## Manual command

```bash
stage20_checks/run_stage20_8_lambda_regression_response_comparison.sh
```

## Expected PASS evidence

```text
STAGE 20.8 LAMBDA REGRESSION RESPONSE COMPARISON VERDICT: PASS
STAGE 20.8 FINAL VERDICT: PASS
```
