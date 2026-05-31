# Stage 13.4 volume normalization audit

## Target

Stage 13.4 adds a standalone grid / volume normalization audit for the future
Eulerian force-density candidate.  It checks that candidate force density is
interpreted conservatively through the volume-weighted integral and that invalid
cell volumes are rejected before any spreading is attempted with them.

## Mathematical / physical meaning

The audited candidate formula is:

```text
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

The conservation target is volume weighted:

```text
sum_i f_i_cand * DeltaV_i ~= sum_k F_k * Delta_s_k
```

It is not the unweighted density sum.  Local force-density values may change when
positive cell volumes are nonuniform, but the volume-weighted integral must remain
conservative.

The production no-fibre DNS remains unchanged:

```text
f_fsi = 0
RHS_stage13.4 = RHS_stage13.3 = RHS_stage13.2 = RHS_stage13.1 = RHS_stage13.0 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

No production RHS injection is performed, no production IBM spreading is called,
and no fibre structure equation is advanced.

## Controlled tests

The Stage 13.4 check covers:

- uniform positive volume;
- nonuniform positive volume;
- invalid zero-volume rejection;
- invalid negative-volume rejection;
- density-times-volume integral consistency;
- component-wise normalization with mixed-sign force components and unequal
  `Delta_s` values;
- near-boundary volume-normalized spreading;
- finite force-density values;
- clear behavior returning diagnostic candidate arrays and norms to zero.

The required tolerances are fixed in the check:

```text
volume_finite_abs_tol          = 0.0
force_conservation_abs_tol     = 1.0e-12
component_conservation_abs_tol = 1.0e-12
boundary_conservation_abs_tol  = 1.0e-12
zero_after_clear_abs_tol       = 1.0e-12
```

## What is intentionally not done

Stage 13.4 intentionally does not add production coupling or physics:

- no `xcompact3d` hook call;
- no production main-loop insertion;
- no production RHS injection;
- no production IBM spreading;
- no feedback force application to fluid;
- no two-way force;
- no fibre structure advance.

## Pass criteria

The Stage 13.4 gate passes only if all of the following hold:

1. `xcompact3d` and the closed Stage 11 / Stage 12 / Stage 13.1 / Stage 13.2 /
   Stage 13.3 check targets build, along with
   `fibre_stage13_volume_normalization_audit_check`.
2. The standalone check reports
   `STAGE 13.4 VOLUME NORMALIZATION AUDIT VERDICT: PASS`.
3. `stage13_outputs/fibre_stage13_4_volume_normalization_audit.dat` reports all
   required status keys as one.
4. Scalar errors satisfy:
   - `stage13_4_uniform_integrated_force_error_l2 <= 1.0e-12`;
   - `stage13_4_nonuniform_integrated_force_error_l2 <= 1.0e-12`;
   - `stage13_4_boundary_integrated_force_error_l2 <= 1.0e-12`;
   - `stage13_4_max_abs_force_density_norm_after_clear <= 1.0e-12`.
5. `stage13_outputs/stage13_4_volume_normalization_audit_gate.dat` reports
   `stage13_4_gate_status 1`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE13_4_RUN_STAGE13_3=0 \
bash stage13_checks/run_stage13_4_volume_normalization_audit.sh
```
