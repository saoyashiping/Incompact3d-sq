# Stage 13.8 force-density parallel consistency

## Target

Stage 13.8 verifies np=1/2/4 no-contamination parallel consistency between:

- a deterministic no-fibre production baseline with Stage 11 one-way sampling and
  Stage 12 Lagrangian feedback-force candidate diagnostics enabled, but Stage 13
  disabled; and
- the same deterministic no-fibre production path with the Stage 13 Eulerian
  force-density candidate diagnostic hook enabled.

The Stage 13 hook must remain active for the candidate-enabled run and produce
finite force-density diagnostics without changing final fluid signatures for
np=1, np=2, or np=4.

## Mathematical / physical meaning

The diagnostic force-density candidate path remains

```text
U_f = I_h[u](X_f)
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

No production forcing is applied:

```text
f_fsi = 0
RHS_stage13.8 = RHS_stage13.7 = RHS_stage13.6 = RHS_stage13.5
              = RHS_stage13.4 = RHS_stage13.3 = RHS_stage13.2
              = RHS_stage13.1 = RHS_stage13.0 = RHS_stage12
              = RHS_stage11 = RHS_stage10 = RHS_stage9
```

The forbidden update remains

```text
RHS <- RHS + f_i_cand
```

## Compared signatures

For each np value (1, 2, and 4), Stage 13.8 compares:

```text
sum_ux, sum_uy, sum_uz,
max_ux, max_uy, max_uz,
l2_ux,  l2_uy,  l2_uz
```

## Tolerance formula

Each baseline-vs-candidate metric must satisfy

```text
abs(delta) <= max(abs_tol, rel_tol * max(1, abs(reference)))
```

with defaults

```text
STAGE13_8_INVARIANCE_ABS_TOL=1.0e-12
STAGE13_8_INVARIANCE_REL_TOL=1.0e-14
```

## Stage 13 diagnostic evidence

The candidate-enabled run must produce
`stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat` and
show:

- force-density candidate computed;
- finite force density;
- finite force-density norm;
- finite integrated force;
- integrated-force conservation;
- spreading input uses the fibre-on-fluid sign;
- wrong-sign rejection;
- no field modification;
- no RHS modification;
- no RHS injection;
- no production IBM forcing;
- no feedback application;
- no two-way force;
- no fibre structure advance.

## Intentionally not done

- no RHS injection;
- no production IBM forcing application;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance;
- no Stage 10 script reruns;
- no closed-stage edits.

## Pass criteria

The gate passes only if all required targets build, the Stage 12 baseline and
Stage 13 candidate-enabled Stage 9.9 deterministic no-fibre runs complete, the
candidate-enabled Stage 9.9 np=1/2/4 consistency verdict remains PASS, Stage
13.6 diagnostics are present with all required values, and every np=1/2/4 final
signature matches within the Stage 13.8 tolerance formula.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE13_8_RUN_STAGE13_7=0 \
STAGE13_8_INVARIANCE_ABS_TOL=1.0e-12 \
STAGE13_8_INVARIANCE_REL_TOL=1.0e-14 \
bash stage13_checks/run_stage13_8_force_density_parallel_consistency.sh
```
