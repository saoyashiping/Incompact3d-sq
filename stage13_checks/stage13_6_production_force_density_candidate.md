# Stage 13.6 production force-density candidate diagnostic hook

## Target

Stage 13.6 adds a guarded production Eulerian force-density candidate diagnostic
hook.  The hook reads production velocity arrays in read-only mode, constructs a
finite diagnostic Eulerian force-density candidate, writes diagnostics, and keeps
all production fields unchanged.

Stage 13.6 does not inject RHS, does not apply force to the fluid, and does not
modify production velocity, pressure, projection, Poisson, RK3, channel forcing,
IBM forcing, or fibre-structure state.

## Mathematical / physical meaning

The production hook uses conservative smoke sampling rather than physical-grid
accurate coupling.  It samples a local production velocity

```text
U_f = (u, v, w)
```

uses a placeholder prescribed velocity

```text
V_f = 0
```

and applies the closed Stage 12 / Stage 13.5 sign convention:

```text
F_fluid_to_fibre_cand = alpha * (U_f - V_f)
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
```

The diagnostic Eulerian force-density candidate is

```text
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

with the conservative volume-weighted target

```text
sum_i f_i_cand * DeltaV_i ~= sum_k F_fibre_to_fluid_cand(k) * Delta_s_k
```

Stage 13.6 remains diagnostic-only:

```text
f_fsi = 0
RHS_stage13.6 = RHS_stage13.5 = RHS_stage13.4 = RHS_stage13.3
              = RHS_stage13.2 = RHS_stage13.1 = RHS_stage13.0
              = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Implementation policy

- conservative production smoke sampling;
- diagnostic force-density candidate only;
- finite force-density signatures;
- finite integrated force diagnostics;
- spreading input sign uses `F_fibre_to_fluid_cand`;
- no physical feedback yet;
- no production RHS injection.

## Intentionally not done

- no RHS injection;
- no production IBM forcing application;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance;
- no pressure/projection/RK3/channel forcing modification.

## Pass criteria

The gate passes only if:

1. all required build targets compile;
2. the standalone Stage 13.6 check prints

   ```text
   STAGE 13.6 PRODUCTION FORCE DENSITY CANDIDATE CHECK VERDICT: PASS
   ```

3. the Stage 9.9 deterministic no-fibre smoke path passes with Stage 11, Stage
   12, and Stage 13 hooks enabled;
4. `stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat`
   reports all required active/finite/conservation/no-coupling statuses as `1`,
   while `stage13_6_field_modified_status` and `stage13_6_rhs_modified_status`
   remain `0`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE13_6_RUN_STAGE13_5=0 \
bash stage13_checks/run_stage13_6_production_force_density_candidate.sh
```
