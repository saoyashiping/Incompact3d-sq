# Stage 12.10 RHS / Spreading / Structure Contamination Audit

## Target

Stage 12.10 performs a strict contamination audit for the Stage 12 production feedback-force candidate hook. The hook may compute diagnostic Lagrangian force-candidate quantities, but it must remain read-only with respect to production fluid fields, RHS construction, IBM spreading, feedback application, two-way force density, and fibre structure advancement.

## Mathematical and physical meaning

The production no-fibre DNS remains

```text
div(u) = 0

du/dt + u · grad(u) = -grad(p) + nu laplacian(u) + f_drive + f_fsi
```

For Stage 12.10,

```text
f_fsi = 0
RHS_stage12.10 = RHS_stage12.9 = RHS_stage12.8 = RHS_stage12.7 = RHS_stage12.6 =
                 RHS_stage12.5 = RHS_stage12.4 = RHS_stage12.3 = RHS_stage12.2 =
                 RHS_stage12.1 = RHS_stage12.0 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

The diagnostic force-candidate quantities are still permitted:

```text
U_f = I_h[u](X_f)
slip = U_f - V_f
F_fluid_to_fibre_cand = alpha * slip
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
```

However, the stage must prove that no Eulerian force density is created, no force is spread, no RHS injection occurs, and no structure advance is called:

```text
f_IBM = 0
RHS <- RHS + f_IBM is forbidden
rho_tilde * X_ddot = d_s(T d_s X) - gamma d_ssss X - F_fs is not advanced
```

## Static audit scope

The gate inspects these files read-only:

```text
src/xcompact3d.f90
src/fibre_stage12_production_feedback_candidate.f90
```

The audit strips Fortran comments, lowercases active statements, and focuses on executable/import evidence only: `use`, `call`, production-field assignment left-hand sides, production velocity dummy-argument intent, and guarded top-level hook placement. Diagnostic key names such as `no_rhs_injection_status`, `no_ibm_spreading_status`, and `field_modified_status` are not treated as contamination evidence.

## Runtime evidence

The runtime evidence is the existing Stage 9.9 deterministic no-fibre smoke path with both hooks enabled:

- Stage 11 one-way sampling hook enabled in read-only mode.
- Stage 12 production feedback-force candidate hook enabled in read-only mode.
- No Stage 10.2 or Stage 10.3 audit logic is run.
- No closed Stage 11 or Stage 12.0--12.9 scripts are run by default.

The Stage 9.9 gate must report `STAGE 9.9 FINAL VERDICT: PASS`, and Stage 12 diagnostics must be written to `stage12_outputs/fibre_stage12_6_production_feedback_candidate.dat`.

## Diagnostic evidence

Stage 12.10 requires the Stage 12.6 diagnostic file to show:

- hook requested, read-only, initialized, and sampled;
- sampled velocity available;
- force candidate computed;
- finite force candidate;
- finite force norm;
- finite power diagnostics;
- no field modification;
- no RHS modification;
- no Eulerian force density;
- no RHS injection;
- no IBM spreading;
- no feedback force application;
- no two-way force;
- no fibre structure advance.

## What is intentionally not done

Stage 12.10 intentionally does not add, enable, or modify:

- Eulerian force density;
- RHS injection;
- IBM force spreading;
- feedback force application to the fluid;
- two-way force density;
- fibre structure advance;
- pressure, projection, Poisson, RK3, or channel forcing logic;
- Stage 10 reruns;
- closed-stage files or scripts.

## Pass criteria

The gate passes only if all of the following are true:

1. Required build targets compile.
2. Static audit finds no forbidden active `use` or `call` evidence in the Stage 12 production hook path.
3. Production velocity dummy arguments in the Stage 12 hook are read-only `intent(in)`.
4. Static assignment checks find no writes to production velocity, pressure, RHS, or Eulerian force-density fields.
5. Stage 9.9 deterministic no-fibre smoke passes with Stage 11 and Stage 12 hooks enabled.
6. Runtime Stage 12 diagnostics prove finite force/power diagnostics and no field, RHS, Eulerian-force, IBM-spreading, feedback-application, two-way-force, or structure-advance contamination.

The gate writes `stage12_outputs/stage12_10_rhs_spreading_structure_contamination_audit.dat` and prints:

```text
STAGE 12.10 RHS SPREADING STRUCTURE CONTAMINATION AUDIT VERDICT: PASS
STAGE 12.10 FINAL VERDICT: PASS
```

Any failure prints an explicit non-empty `Reasons:` line. Missing Stage 12 feedback-candidate diagnostics fail with `missing_stage12_6_feedback_candidate_diagnostics`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE12_10_RUN_STAGE12_9=0 \
bash stage12_checks/run_stage12_10_rhs_spreading_structure_contamination_audit.sh
```
