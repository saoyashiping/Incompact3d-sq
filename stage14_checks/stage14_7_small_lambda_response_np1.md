# Stage 14.7 — Small-Lambda Response at np=1

## Target

Stage 14.7 is the first controlled production response gate for a small nonzero Stage 14 RHS-injection gain.

The gate compares:

- a lambda-zero baseline with Stage 11 one-way sampling, Stage 12 feedback-force candidate diagnostics, Stage 13 Eulerian force-density candidate diagnostics, and the Stage 14 RHS hook enabled with `lambda_14=0`;
- a small-lambda candidate with the same Stage 11/12/13/14 path enabled and `lambda_14=epsilon`.

The default epsilon is deliberately tiny:

```text
STAGE14_7_SMALL_LAMBDA=1.0e-8
```

The candidate must show an active hook, nonzero finite bounded RHS increment, and finite bounded final-fluid response.

## Mathematical / physical meaning

Stage 13 provides the Eulerian force-density candidate:

```text
f_i_cand = S_h[F_fibre_to_fluid_cand]
```

The intended Stage 14 production RHS form is:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

Stage 14.6 validated the lambda-zero case:

```text
lambda_14 = 0 -> RHS_new = RHS_old
```

Stage 14.7 validates the controlled small-lambda case:

```text
lambda_14 = epsilon, epsilon > 0, epsilon << 1
Delta_RHS = epsilon * f_i_cand
```

Expected behavior:

- `Delta_RHS` is nonzero;
- `Delta_RHS` is finite;
- `Delta_RHS` uses the fibre-on-fluid Stage 13 candidate sign through the existing Stage 14 production hook path;
- `Delta_RHS` is bounded for the selected epsilon;
- final fluid signatures remain finite and bounded.

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

The lambda-zero copy is written to:

```text
stage14_outputs/stage14_7_lambda0_np1.dat
```

The small-lambda copy is written to:

```text
stage14_outputs/stage14_7_small_lambda_np1.dat
```

For Stage 14.7, exact invariance is not expected. The gate instead requires every fluid-signature delta to be finite and bounded by:

```text
STAGE14_7_MAX_FLUID_DELTA=1.0e-4
```

## Diagnostic evidence

The lambda-zero run must produce:

```text
stage14_outputs/stage14_7_lambda0_rhs_hook.dat
```

The small-lambda candidate must produce:

```text
stage14_outputs/stage14_7_small_lambda_rhs_hook.dat
```

The small-lambda diagnostic evidence must show:

- hook active;
- lambda nonzero (`stage14_5_lambda_zero_status = 0`);
- RHS increment nonzero;
- RHS increment finite;
- reported injection gain matches `STAGE14_7_SMALL_LAMBDA` within `STAGE14_7_RESPONSE_ABS_TOL` / `STAGE14_7_RESPONSE_REL_TOL`;
- RHS increment bounded by `STAGE14_7_MAX_RHS_INCREMENT`, default `1.0e-4`;
- final fluid delta finite and bounded;
- no pressure/projection/Poisson/RK3/channel-forcing modification;
- no production IBM forcing;
- no feedback application beyond Stage 13 diagnostics;
- no two-way force;
- no structure advance.

If the Stage 14 diagnostic file is missing, the gate fails explicitly with:

```text
missing_stage14_5_production_rhs_hook_diagnostics
```

If the closed Stage 14.5 hook still blocks all nonzero lambda, the gate fails explicitly with:

```text
nonzero_lambda_still_blocked_by_stage14_5
```

That failure means a future allowed Stage 14.7 source patch is required to permit strictly gated small-lambda production response. This Stage 14.7 script does not modify closed Stage 14.5 files.

## What is intentionally not done

Stage 14.7 intentionally does not:

- use a large lambda;
- perform the full pressure-projection/divergence regression reserved for Stage 14.8;
- perform np=1/2/4 parallel RHS-injection consistency;
- perform restart/stats/visu/coarse-I/O checks;
- advance the fibre structure;
- edit source files;
- edit closed-stage files.

## Pass criteria

The Stage 14.7 gate passes only if:

- all required build targets compile;
- the lambda-zero Stage 9.9 smoke reports `STAGE 9.9 FINAL VERDICT: PASS`;
- the small-lambda Stage 9.9 smoke reports `STAGE 9.9 FINAL VERDICT: PASS`;
- the lambda-zero Stage 14 hook diagnostics confirm zero increment and unchanged RHS;
- the small-lambda diagnostics confirm active hook, nonzero finite bounded RHS increment, and no forbidden modification/coupling statuses;
- all np=1 final fluid-signature deltas are finite and bounded;
- the final gate file reports `stage14_7_small_lambda_response_np1_status 1`.

The final gate writes:

```text
stage14_outputs/stage14_7_small_lambda_response_np1.dat
```

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_7_RUN_STAGE14_6=0 \
STAGE14_7_SMALL_LAMBDA=1.0e-8 \
STAGE14_7_MAX_RHS_INCREMENT=1.0e-4 \
STAGE14_7_MAX_FLUID_DELTA=1.0e-4 \
bash stage14_checks/run_stage14_7_small_lambda_response_np1.sh
```
