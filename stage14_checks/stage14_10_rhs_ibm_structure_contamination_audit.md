# Stage 14.10 — RHS / IBM / Structure Contamination Audit

## Objective

Stage 14.10 adds a static plus runtime audit path for controlled small-lambda Stage 14 RHS injection. The audit verifies that the only intended Stage 14 production-side update remains:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

where `f_i_cand` is the already-audited Stage 13 Eulerian force-density candidate using the fibre-on-fluid sign convention.

Stage 14.10 is a contamination audit stage. It does not introduce new physics, does not reimplement Stage 13 force-density construction, does not reimplement Stage 14 RHS injection, and does not advance the fibre structure.

## Preservation checks from Stage 14.8 / Stage 14.9

The wrapper explicitly fails if any of these regressions are detected:

- the old small-lambda blocker `stage14_get_injection_gain() == 0.0` appears in the Stage 14 hook-registration path;
- the Stage 14 production RHS hook is disconnected from `xcompact3d` or no longer depends on Stage 13 and Stage 14 request/enable controls;
- Stage 13 or Stage 14 production diagnostics are not rank0-safe/race-free;
- Stage 13 force-density diagnostics revert to local subdomain-center sampling that can cause `np=1/2/4` mismatches;
- small-lambda Stage 14 runtime diagnostics are missing, zero, non-finite, or unbounded.

## Static audit coverage

The static audit checks:

1. `src/xcompact3d.f90` does not gate Stage 14 hook registration on `stage14_get_injection_gain() == 0.0`.
2. The Stage 14 RHS hook remains connected at the production RHS injection site.
3. Stage 14 references do not appear in pressure/projection/Poisson/RK3/channel-forcing source files beyond the configured strict match allowance.
4. Stage 14 production code does not actively `use` or `call` production IBM forcing, pressure/projection/Poisson, RK3/channel-forcing, or structure-advance/bending/tension/wall-contact routines.
5. Stage 13 and Stage 14 production diagnostic writing remains rank0-safe or otherwise race-free.
6. Stage 13 force-density sampling does not use the local subdomain center pattern that previously made diagnostics decomposition-dependent.

The static report is written to:

```text
stage14_outputs/stage14_10_static_audit_report.txt
```

## Runtime audit coverage

The runtime audit performs two production smoke cases at `STAGE14_10_NP` ranks:

| Case | Lambda | Purpose |
| --- | --- | --- |
| `lambda0` | `0.0` | Preserve lambda-zero no-contamination evidence. |
| `small_lambda` | `STAGE14_10_SMALL_LAMBDA` | Verify active small-lambda RHS injection and no IBM/structure contamination. |

Both cases enable the full diagnostic chain:

```text
Stage 11 one-way sampling
Stage 12 Lagrangian feedback-force candidate
Stage 13 Eulerian force-density candidate
Stage 14 controlled RHS injection
```

The runtime audit uses the Stage 9.9 deterministic no-fibre production path as the bounded fluid-signature smoke. It fails closed if Stage 11, Stage 12, Stage 13, or Stage 14 diagnostics are missing.

## Required artifacts

The wrapper writes the final Stage 14.10 gate file:

```text
stage14_outputs/stage14_10_rhs_ibm_structure_contamination_audit.dat
```

It also records/copies:

```text
stage14_outputs/stage14_10_lambda0_runtime.log
stage14_outputs/stage14_10_small_lambda_runtime.log
stage14_outputs/stage14_10_lambda0_stage9_9.dat
stage14_outputs/stage14_10_small_lambda_stage9_9.dat
stage14_outputs/stage14_10_small_lambda_stage11_oneway.dat
stage14_outputs/stage14_10_small_lambda_stage12_feedback_candidate.dat
stage14_outputs/stage14_10_small_lambda_stage13_force_density.dat
stage14_outputs/stage14_10_lambda0_rhs_hook.dat
stage14_outputs/stage14_10_small_lambda_rhs_hook.dat
stage14_outputs/stage14_10_static_audit_report.txt
```

## Pass criteria

The final verdict passes only if all of the following are true:

- static audit status is `1`;
- forbidden lambda-zero registration gate is absent;
- Stage 14 hook is connected to the RHS path;
- pressure/projection/Poisson/RK3/channel-forcing static contamination is absent;
- production IBM forcing static contamination is absent;
- structure/bending/tension/wall-contact static contamination is absent;
- Stage 13 and Stage 14 diagnostics remain rank0-safe/race-free;
- Stage 13 force-density np-sampling repair is preserved;
- lambda-zero runtime case passes and reports zero RHS increment;
- small-lambda runtime case passes and reports nonzero finite bounded RHS increment;
- Stage 11, Stage 12, Stage 13, and Stage 14 diagnostics are active;
- final fluid-signature deltas are finite and bounded by `STAGE14_10_MAX_FLUID_DELTA`;
- no pressure/projection/Poisson/RK3/channel-forcing contamination diagnostics are false;
- no production IBM forcing, feedback application, two-way force, or structure advance diagnostics are false;
- no bending solve, tension solve, fibre position update, structural velocity update, or wall/contact evidence is present.

## Manual command

Run from the repository root:

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_10_RUN_STAGE14_9=0 \
STAGE14_10_SMALL_LAMBDA=1.0e-8 \
STAGE14_10_MAX_RHS_INCREMENT=1.0e-4 \
STAGE14_10_MAX_FLUID_DELTA=1.0e-4 \
STAGE14_10_MAX_STATIC_MATCHES=0 \
STAGE14_10_NP=2 \
bash stage14_checks/run_stage14_10_rhs_ibm_structure_contamination_audit.sh
```

The script may also be run with defaults exactly as:

```sh
bash stage14_checks/run_stage14_10_rhs_ibm_structure_contamination_audit.sh
```

## Expected verdict

A successful run prints:

```text
STAGE 14.10 RHS/IBM/STRUCTURE CONTAMINATION AUDIT VERDICT: PASS
STAGE 14.10 FINAL VERDICT: PASS
```
