# Stage 14.6 — Lambda-Zero No-Contamination Invariance

## Target

Stage 14.6 verifies lambda-zero no-contamination invariance for the guarded Stage 14 production RHS hook.

The gate compares:

- a Stage 13 baseline production smoke with Stage 11 one-way sampling, Stage 12 feedback-force candidate diagnostics, and Stage 13 Eulerian force-density candidate diagnostics enabled, but Stage 14 disabled;
- a Stage 14 candidate production smoke with the same Stage 11/12/13 diagnostics plus the Stage 14 RHS hook enabled with `lambda_14=0`.

The Stage 14 hook must be active in the candidate run, but the RHS increment must be zero.

## Mathematical / physical meaning

Stage 13 provides the Eulerian force-density candidate:

```text
f_i_cand = S_h[F_fibre_to_fluid_cand]
```

The intended Stage 14 production RHS form is:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

At Stage 14.6:

```text
lambda_14 = 0
RHS_new = RHS_old
Delta_RHS = 0
Delta_u = 0
```

Therefore the candidate run must match the Stage 13 baseline final fluid signatures. Stage 14.6 does not advance the structure equation.

## Compared signatures

The gate compares these np=1 final signatures from `stage9_outputs/fibre_stage9_9_parallel_consistency_np1.dat`:

```text
stage9_9_signature_sum_ux
stage9_9_signature_sum_uy
stage9_9_signature_sum_uz
stage9_9_signature_max_ux
stage9_9_signature_max_uy
stage9_9_signature_max_uz
stage9_9_signature_l2_ux
stage9_9_signature_l2_uy
stage9_9_signature_l2_uz
```

The baseline copy is written to:

```text
stage14_outputs/stage14_6_baseline_stage13_np1.dat
```

The candidate copy is written to:

```text
stage14_outputs/stage14_6_candidate_lambda0_np1.dat
```

## Tolerance formula

For each metric, the gate uses the mixed tolerance:

```text
abs(delta) <= max(abs_tol, rel_tol * max(1, abs(reference)))
```

The defaults are:

```text
STAGE14_6_INVARIANCE_ABS_TOL=1.0e-12
STAGE14_6_INVARIANCE_REL_TOL=1.0e-14
```

Every metric comparison prints the metric name, baseline value, candidate value, delta, absolute tolerance, relative tolerance, effective tolerance, and PASS/FAIL verdict.

## Stage 14 diagnostic evidence

The candidate run must produce:

```text
stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
```

The Stage 14 diagnostic file must prove:

- hook active;
- lambda zero;
- RHS increment zero;
- RHS unchanged;
- no pressure modification;
- no projection modification;
- no Poisson modification;
- no RK3 modification;
- no channel-forcing modification;
- no production IBM forcing;
- no feedback application;
- no two-way force;
- no structure advance.

If the diagnostic file is missing, the gate fails explicitly with:

```text
missing_stage14_5_production_rhs_hook_diagnostics
```

## What is intentionally not done

Stage 14.6 intentionally does not:

- test small-lambda production response;
- modify source code;
- modify pressure, projection, Poisson, RK3, or channel forcing source;
- call production IBM forcing application;
- activate two-way force;
- advance the fibre structure;
- edit any closed-stage file.

## Pass criteria

The Stage 14.6 gate passes only if:

- all required build targets compile;
- the Stage 13 baseline Stage 9.9 smoke reports `STAGE 9.9 FINAL VERDICT: PASS`;
- the Stage 14 lambda-zero candidate Stage 9.9 smoke reports `STAGE 9.9 FINAL VERDICT: PASS`;
- all np=1 final signatures are invariant within the mixed tolerance;
- the Stage 14 hook diagnostic file exists;
- Stage 14 diagnostics confirm hook activity, lambda zero, zero RHS increment, unchanged RHS signatures, and all no-modification/no-coupling statuses;
- the final gate file reports `stage14_6_lambda_zero_invariance_status 1`.

The final gate writes:

```text
stage14_outputs/stage14_6_lambda_zero_invariance.dat
```

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_6_RUN_STAGE14_5=0 \
STAGE14_6_INVARIANCE_ABS_TOL=1.0e-12 \
STAGE14_6_INVARIANCE_REL_TOL=1.0e-14 \
bash stage14_checks/run_stage14_6_lambda_zero_invariance.sh
```
