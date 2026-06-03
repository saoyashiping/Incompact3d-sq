# Stage 14.11 Total Smoke and Closure

## Objective

Stage 14.11 is the Stage 14 total-smoke and closure evidence path. It verifies that the controlled RHS-injection chain remains valid and generates `stage14_checks/STAGE14_CLOSED.md` only after every required Stage 14.11 check passes.

The chain under closure is:

1. Stage 11 one-way sampling.
2. Stage 12 Lagrangian feedback-force candidate.
3. Stage 13 Eulerian force-density candidate.
4. Stage 14 controlled RHS injection.

Stage 14.11 introduces no new physics and does not advance the flexible-fibre structure.

## Controlled RHS formula

The Stage 14 controlled update remains:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

Here `f_i_cand` is the Stage 13 Eulerian force-density candidate using the already-audited fibre-on-fluid sign convention.

## Required regression prohibitions

Stage 14.11 fails if any of these regressions are detected:

- The forbidden Stage 14 hook-registration gate `stage14_get_injection_gain() == 0.0` is present.
- Small-lambda Stage 14 hook registration is blocked.
- Stage 11.5, Stage 13, or Stage 14 production diagnostics are missing or not rank0-safe/race-free.
- Stage 13 force-density diagnostics revert to local subdomain center sampling.
- Stage 14 touches pressure/projection/Poisson, RK3, or channel forcing code.
- Production IBM forcing is activated outside the approved Stage 13/14 diagnostic-to-RHS path.
- Any structure advance, bending solve, tension solve, fibre position update, structural velocity update, or wall/contact logic is active in the Stage 14 production path.
- Required PASS evidence is silently skipped.

## Checks performed

### Static regression audit

The wrapper checks that:

- Stage 14 hook registration is not gated by `lambda_14 == 0`.
- The production Stage 14 RHS hook is connected to the Stage 13 force-density candidate and Stage 14 request/enable controls.
- Stage 11.5, Stage 12.6, Stage 13.6, and Stage 14.5 diagnostic write paths exist.
- Stage 11/13/14 production diagnostic writes are rank0-safe or otherwise race-free.
- Stage 13 force-density diagnostics do not use local subdomain center sampling.
- Stage 14 code has no pressure/projection/Poisson, RK3, channel-forcing, legacy IBM, or structure-advance contamination.

### Runtime total smoke

The wrapper runs one small-lambda Stage 14 production smoke case with `STAGE14_11_SMALL_LAMBDA` and `STAGE14_11_NP` MPI ranks. It requires:

- Stage 11 one-way sampling active.
- Stage 12 feedback-force candidate active.
- Stage 13 Eulerian force-density candidate active.
- Stage 14 RHS hook active.
- `lambda_14` recorded and nonzero.
- RHS increment nonzero, finite, sign-correct, and bounded by `STAGE14_11_MAX_RHS_INCREMENT`.
- Fluid response bounded by `STAGE14_11_MAX_FLUID_DELTA`.
- No NaN/Inf in runtime logs or diagnostics.

### Prior Stage 14 evidence preservation

Stage 14.11 verifies existing PASS datfiles by default and can rerun closed prerequisites only when explicitly requested:

- `STAGE14_11_RUN_STAGE14_8=1` reruns the Stage 14.8 np=1/2/4 parallel small-lambda response gate.
- `STAGE14_11_RUN_STAGE14_9=1` reruns the Stage 14.9 restart/statistics/visualization/coarse-I/O compatibility gate.
- `STAGE14_11_RUN_STAGE14_10=1` reruns the Stage 14.10 RHS/IBM/structure contamination audit.

By default, all three flags are `0`, so Stage 14.11 fails closed if the required prior evidence datfiles are missing or failed.

## Required artifacts

A passing run produces:

- `stage14_outputs/stage14_11_total_smoke_closure.dat`
- `stage14_outputs/stage14_11_static_audit_report.txt`
- `stage14_outputs/stage14_11_total_smoke_runtime.log`
- `stage14_outputs/stage14_11_stage11_oneway.dat`
- `stage14_outputs/stage14_11_stage12_feedback_candidate.dat`
- `stage14_outputs/stage14_11_stage13_force_density.dat`
- `stage14_outputs/stage14_11_stage14_rhs_hook.dat`
- `stage14_checks/STAGE14_CLOSED.md`

The closure file is generated only after the total-smoke verdict is fully PASS.

## Pass criteria

The final status requires all of the following:

- Static regression audit PASS.
- Runtime total-smoke PASS.
- Required Stage 11/12/13/14 diagnostics present and internally PASS.
- Small-lambda Stage 14 RHS hook active.
- Small-lambda RHS increment nonzero, finite, sign-correct, and bounded.
- np=1/2/4 Stage 14.8 parallel consistency evidence present and PASS.
- Stage 14.9 I/O/restart evidence present and PASS.
- Stage 14.10 contamination-audit evidence present and PASS.
- No pressure/projection/Poisson/RK3/channel-forcing contamination.
- No production IBM forcing outside the approved Stage 13/14 path.
- No fibre structure advance or structural solver path.

The wrapper prints:

```text
STAGE 14.11 TOTAL SMOKE VERDICT: PASS
STAGE 14.11 FINAL VERDICT: PASS
STAGE14_CLOSED=YES
```

on success, or FAIL with explicit reasons.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_11_RUN_STAGE14_8=0 \
STAGE14_11_RUN_STAGE14_9=0 \
STAGE14_11_RUN_STAGE14_10=0 \
STAGE14_11_SMALL_LAMBDA=1.0e-8 \
STAGE14_11_MAX_RHS_INCREMENT=1.0e-4 \
STAGE14_11_MAX_FLUID_DELTA=1.0e-4 \
STAGE14_11_MAX_PARALLEL_RHS_DIFF=1.0e-12 \
STAGE14_11_MAX_PARALLEL_FORCE_DENSITY_DIFF=1.0e-10 \
STAGE14_11_MAX_RESTART_DELTA=1.0e-8 \
STAGE14_11_MAX_IO_SIGNATURE_DELTA=1.0e-8 \
STAGE14_11_NP=2 \
STAGE14_11_NP_LIST="1 2 4" \
bash stage14_checks/run_stage14_11_total_smoke_closure.sh
```

If prior Stage 14.8/14.9/14.10 datfiles are not already present, rerun them explicitly by setting the corresponding `STAGE14_11_RUN_STAGE14_*` flag(s) to `1`.
