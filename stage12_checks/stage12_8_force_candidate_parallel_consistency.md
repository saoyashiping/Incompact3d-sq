# Stage 12.8 force-candidate parallel consistency

## Target

Stage 12.8 verifies np=1/2/4 no-contamination parallel consistency by comparing deterministic no-fibre production runs with Stage 11 one-way sampling enabled against matching runs with both Stage 11 one-way sampling and the Stage 12 production feedback-force candidate hook enabled.

## Mathematical / physical meaning

The no-fibre production DNS remains unchanged:

```text
U_f = I_h[u](X_f)
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
f_fsi = 0
RHS_stage12.8 = RHS_stage12.7 = RHS_stage12.6 = RHS_stage12.5 = RHS_stage12.4 = RHS_stage12.3 = RHS_stage12.2 = RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

Stage 12.8 allows diagnostic Lagrangian force-candidate and power-diagnostic computations, but performs no spreading and no RHS injection.

## Compared signatures

For each of np=1, np=2, and np=4, the gate compares the final Stage 9.9 signatures from the baseline and candidate-enabled runs:

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

## Tolerance formula

Each metric must satisfy:

```text
abs(delta) <= max(abs_tol, rel_tol * max(1, abs(reference)))
```

Defaults:

```sh
STAGE12_8_INVARIANCE_ABS_TOL=${STAGE12_8_INVARIANCE_ABS_TOL:-1.0e-12}
STAGE12_8_INVARIANCE_REL_TOL=${STAGE12_8_INVARIANCE_REL_TOL:-1.0e-14}
```

## Stage 12 diagnostic evidence

The candidate-enabled run must write `stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat` proving:

- the Stage 12 hook was requested, initialized, and sampled;
- the force candidate was computed;
- force, force norm, and power diagnostics are finite;
- production fields were not modified;
- RHS was not modified;
- no Eulerian force density was created;
- no RHS injection occurred;
- no IBM spreading occurred;
- no feedback application, two-way force, or structure advance occurred.

## Intentionally not done

- No Eulerian force density.
- No RHS injection.
- No IBM spreading.
- No feedback force application to fluid.
- No two-way force.
- No fibre structure advance.
- No Stage 10 script reruns.
- No closed-stage edits.

## Pass criteria

The Stage 12.8 gate passes only if the required build targets compile, baseline and candidate Stage 9.9 runs produce np=1/2/4 dat files with local and initial-decomposition statuses passing, the candidate Stage 9.9 log reports `STAGE 9.9 FINAL VERDICT: PASS`, the Stage 12.6 diagnostic statuses prove active finite no-coupling behavior, and all np=1/2/4 final fluid signatures are invariant within tolerance.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE12_8_RUN_STAGE12_7=0 \
STAGE12_8_INVARIANCE_ABS_TOL=1.0e-12 \
STAGE12_8_INVARIANCE_REL_TOL=1.0e-14 \
bash stage12_checks/run_stage12_8_force_candidate_parallel_consistency.sh
```
