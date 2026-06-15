# Stage 20.6 RHS coupling lambda gate

Stage 20.6 is a diagnostic-only RHS-coupling lambda-gate boundary. It constructs
helper-local RHS candidate arrays from the Stage 20.5 helper-local Eulerian
force-density candidate and verifies lambda-gated behavior. It does not write to
production RHS, call Stage 14 RHS injection, run production DNS, or activate
production two-way coupling.

## Source-only and previous-stage policy

Stage 20.6 accepts Stage 20.5 PASS evidence when present and preserves Stage
20.0/20.1/20.2/20.3/20.4/20.5 source-only acceptance behavior. It does not
rerun Stage 10-19 or Stage 20.0-20.5, does not require old closure files, and
does not require all old stage output directories in source-only archives.

## Helper-local lambda cases

```text
f_eulerian_effective_zero = 0.0 * f_eulerian_candidate
RHS_delta_zero = f_eulerian_effective_zero
RHS_after_zero = RHS_before + RHS_delta_zero
RHS_zero_residual = RHS_after_zero - RHS_before

f_eulerian_effective_small = lambda_small * f_eulerian_candidate
RHS_delta_small = f_eulerian_effective_small
RHS_after_small = RHS_before + RHS_delta_small
RHS_small_formula_residual = RHS_delta_small - lambda_small * f_eulerian_candidate
```

The zero-lambda case is a strict RHS no-op. The small-lambda case produces only
a bounded helper-local RHS delta equal to `lambda_small * f_eulerian_candidate`.

## Safety boundary

Stage 20.6 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 or Stage 20.0-20.5 file. It does not run MPI, DNS, previous
stages, builds, or production validation.

## Next stage

Stage 20.7: controlled one-fibre closed-loop response np=1.

## Manual command

```bash
stage20_checks/run_stage20_6_rhs_coupling_lambda_gate.sh
```

Expected output includes:

```text
STAGE 20.6 RHS COUPLING LAMBDA GATE VERDICT: PASS
STAGE 20.6 FINAL VERDICT: PASS
```
