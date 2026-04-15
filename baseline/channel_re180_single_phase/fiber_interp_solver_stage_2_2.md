# Stage 2.2 supplement: readonly interpolation on real solver grid/field

## What this supplement adds
This supplement adds a **solver-based readonly interpolation verification mode**.

- Uses real solver grid coordinates (`xp`, `yp`, `zp`) for interpolation.
- Uses real solver velocity fields (`ux1`, `uy1`, `uz1`) as interpolation source.
- Interpolates onto static single-fiber points and writes `fiber_interp_solver.dat`.
- Keeps interpolation strictly readonly: no feedback to flow equations.

## Why this is still Stage 2.2
This is still operator verification:
- no spreading,
- no IBM forcing,
- no no-slip closure,
- no rigid/flexible dynamics,
- no multi-fiber logic.

## Mode separation
- `interp_test_active`: synthetic operator tests (constant/linear fields).
- `interp_solver_test_active`: real-solver readonly interpolation verification.

These two switches are separate and cannot be enabled together.

## How to run
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_interp_solver_test.i3d
```

## Outputs to check
- `fiber_points.dat`
- `fiber_interp_solver.dat` (columns: `itime index x y z u_interp v_interp w_interp sumw`)
- console message with output `itime` and `sumw_min/max`

## Current boundaries
- x/z periodic directions use minimum-image distances.
- y direction remains non-periodic (uses solver stretched `yp`).
- Current solver-based readonly verification is implemented for single-process validation path.

## stretched y-grid fix
- Previous solver-readonly path used a single global `hy = yp(2)-yp(1)`; this is not valid for stretched `yp`.
- Now solver-readonly path uses local y-scale `hy_loc(j)`:
  - boundary points: one-sided spacing
  - interior points: `0.5*(yp(j+1)-yp(j-1))`
- The same `hy_loc(j)` is used in:
  - y support check (`|ry| <= 2*hy_loc`)
  - kernel scale in y
  - y part of volume weight (`hx*hy_loc*hz`)
- x/z still use minimum-image periodic distance; y remains non-periodic.
