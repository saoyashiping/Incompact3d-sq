# Stage 13.3 controlled spreading kernel checks

## Target

Stage 13.3 adds controlled standalone Lagrangian-to-Eulerian spreading kernel
checks for the future Stage 13 Eulerian force-density candidate path.  The stage
is analytic and diagnostic only: it does not connect the spreading kernel to the
production solver or the main time-stepping loop.

## Mathematical / physical meaning

The controlled standalone formula checked in Stage 13.3 is:

```text
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

The conservation target is:

```text
sum_i f_i_cand * DeltaV_i ~= sum_k F_k * Delta_s_k
```

Stage 13.3 computes this only in standalone controlled checks with compact
trilinear weights on a small Cartesian grid.  The production no-fibre DNS remains
unchanged:

```text
f_fsi = 0
RHS_stage13.3 = RHS_stage13.2 = RHS_stage13.1 = RHS_stage13.0 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

No production RHS injection is performed, no production IBM spreading is called,
and no fibre structure equation is advanced.

## Controlled tests

The Stage 13.3 check covers:

- zero-force spreading;
- single-point finite spreading;
- compact support with at most eight trilinear nodes per Lagrangian point;
- nonnegative weights;
- weight normalization;
- component-wise spreading;
- boundary-safe in-domain spreading;
- integrated-force conservation;
- clear behavior returning diagnostic buffers to zero.

The required tolerances are fixed in the check:

```text
weight_sum_abs_tol         = 1.0e-12
force_conservation_abs_tol = 1.0e-12
zero_force_abs_tol         = 1.0e-12
```

## What is intentionally not done

Stage 13.3 intentionally does not add production coupling or physics:

- no `xcompact3d` hook call;
- no production main-loop insertion;
- no production RHS injection;
- no production IBM spreading;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance.

## Pass criteria

The Stage 13.3 gate passes only if all of the following hold:

1. `xcompact3d` and the closed Stage 11 / Stage 12 / Stage 13.1 / Stage 13.2
   check targets build, along with `fibre_stage13_spreading_kernel_check`.
2. The standalone check reports `STAGE 13.3 SPREADING KERNEL VERDICT: PASS`.
3. `stage13_outputs/fibre_stage13_3_spreading_kernel.dat` reports all required
   status keys as one.
4. Scalar errors satisfy:
   - `stage13_3_zero_force_max_abs <= 1.0e-12`;
   - `stage13_3_max_weight_sum_error <= 1.0e-12`;
   - `stage13_3_integrated_force_error_l2 <= 1.0e-12`;
   - `stage13_3_max_abs_force_density_norm_after_clear <= 1.0e-12`.
5. `stage13_3_single_point_support_count <= 8` for the controlled trilinear
   spreading kernel.
6. `stage13_outputs/stage13_3_spreading_kernel_gate.dat` reports
   `stage13_3_gate_status 1`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE13_3_RUN_STAGE13_2=0 \
bash stage13_checks/run_stage13_3_spreading_kernel.sh
```
