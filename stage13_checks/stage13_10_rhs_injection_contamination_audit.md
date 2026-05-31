# Stage 13.10 RHS injection contamination audit

## Target

Stage 13.10 performs a strict RHS injection / production IBM forcing / feedback /
two-way / structure contamination audit for the Stage 13 Eulerian force-density
candidate diagnostic hook.

The hook may compute diagnostic Eulerian force-density candidates, but it must
not activate production RHS injection, production IBM forcing application,
feedback force application to the fluid, two-way coupling, fibre structure
advance, or pressure / projection / Poisson / RK3 / channel-forcing changes.

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
RHS_stage13.10 = RHS_stage13.9 = RHS_stage13.8 = RHS_stage13.7
               = RHS_stage13.6 = RHS_stage13.5 = RHS_stage13.4
               = RHS_stage13.3 = RHS_stage13.2 = RHS_stage13.1
               = RHS_stage13.0 = RHS_stage12 = RHS_stage11
               = RHS_stage10 = RHS_stage9
```

The forbidden update remains

```text
RHS <- RHS + f_i_cand
```

Stage 13.10 also requires no production IBM forcing, no feedback application, no
two-way coupling, and no fibre structure advance.

## Static audit scope

The required read-only static audit scope is:

```text
src/xcompact3d.f90
src/fibre_stage13_production_force_density_candidate.f90
```

The audit strips full-line comments, strips trailing inline comments, lowercases
active statements, and then inspects active `use` / `call` statements, writable
`intent(out/inout)` declarations, production-field assignment left-hand sides,
and top-level hook placement. Diagnostic `no_*` keys, comments, output-key
strings, and diagnostic write lines are not treated as contamination.

## Runtime evidence

Runtime evidence is the Stage 9.9 deterministic no-fibre smoke with Stage 11,
Stage 12, and Stage 13 hooks enabled. The runtime gate must report:

```text
STAGE 9.9 FINAL VERDICT: PASS
```

Stage 13.10 then verifies
`stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat`.

## Diagnostic evidence

The Stage 13 diagnostic file must show:

- force-density candidate computed;
- finite force-density;
- finite force-density norm;
- finite integrated force;
- integrated-force conservation;
- fibre-on-fluid sign;
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

The gate passes only if all required targets build, the static audit finds no
forbidden active calls/imports/writes/intents and confirms guarded hook placement,
Stage 9.9 passes with Stage 11 + Stage 12 + Stage 13 hooks enabled, the Stage
13.6 diagnostic file is present with every required value, and the final Stage
13.10 `.dat` reports no RHS injection, no production IBM forcing, no feedback,
no two-way coupling, and no structure advance.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE13_10_RUN_STAGE13_9=0 \
bash stage13_checks/run_stage13_10_rhs_injection_contamination_audit.sh
```
