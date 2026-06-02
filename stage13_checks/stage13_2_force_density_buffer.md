# Stage 13.2 force-density candidate buffer skeleton

## Target

Stage 13.2 adds the Eulerian force-density candidate buffer skeleton for the future
Stage 13 spreading path.  It provides zero-valued containers for the candidate
Eulerian force-density components and diagnostics only.

## Mathematical / physical meaning

Later Stage 13 work will construct the Eulerian force-density candidate from the
Stage 12 Lagrangian fibre-on-fluid candidate force using a discrete spreading
operator, for example:

```text
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

Stage 13.2 does not compute this quantity.  It only allocates, initializes,
clears, and validates zero-valued force-density candidate buffers.

The production no-fibre DNS remains unchanged:

```text
f_fsi = 0
RHS_stage13.2 = RHS_stage13.1 = RHS_stage13.0 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

No spreading is performed, no RHS contribution is formed, and no structure
equation is advanced.

## State variables

The Stage 13.2 buffer module owns only diagnostic candidate storage:

```text
fx_cand(i,j,k)
fy_cand(i,j,k)
fz_cand(i,j,k)
force_density_norm(i,j,k)
force_density_valid_flag(i,j,k)
```

The arrays are initialized and cleared to zero, while the validity flags are kept
at one for allocated, valid buffer entries.

## What is intentionally not done

Stage 13.2 intentionally does not add any production physics or coupling:

- no `xcompact3d` hook call;
- no production main-loop insertion;
- no Lagrangian-to-Eulerian spreading;
- no force-density computation from Lagrangian force;
- no RHS injection;
- no IBM spreading;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance.

## Pass criteria

The Stage 13.2 gate passes only if all of the following hold:

1. `xcompact3d` and the closed Stage 11 / Stage 12 / Stage 13.1 check targets build.
2. `fibre_stage13_force_density_buffer_check` builds and runs in controlled
   readonly diagnostic mode.
3. The check reports `STAGE 13.2 FORCE DENSITY BUFFER VERDICT: PASS`.
4. `stage13_outputs/fibre_stage13_2_force_density_buffer.dat` reports allocated
   valid buffers with accepted shape, zero force-density components, finite norm
   values, valid flags, and clear status.
5. All no-spreading, no-RHS-injection, no-IBM-spreading, no-feedback-application,
   no-two-way-force, no-structure-advance, no-fluid-field-access, and
   no-fluid-field-modification statuses are one.
6. `stage13_outputs/stage13_2_force_density_buffer_gate.dat` reports
   `stage13_2_gate_status 1`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE13_2_RUN_STAGE13_1=0 \
bash stage13_checks/run_stage13_2_force_density_buffer.sh
```
