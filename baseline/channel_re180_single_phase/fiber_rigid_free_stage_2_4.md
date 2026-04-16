# Stage 2.4: single rigid-fiber free dynamics validation

## Scope

Stage 2.4 adds **single-fiber rigid-body free response** on top of the existing Euler↔Lagrange coupling chain.

- still single fiber
- still rigid straight centerline
- still no-slip penalty coupling `Ffs = beta (U - Xdot)`
- still uses existing interpolation and spreading operators

This stage is **not** the flexible-fiber stage.

## Model used in 2.4

State variables:

- `fiber_xc(3)`: center position
- `fiber_uc(3)`: center velocity
- `fiber_p(3)`: unit direction vector
- `fiber_omega(3)`: angular velocity used to evolve `p`
- `fiber_force_total(3)`: total coupled force
- `fiber_torque_total(3)`: total coupled torque

Geometry is rebuilt every step with:

`X_l = X_c + s_ref(l) * p`

where `s_ref` is the reference arc-length coordinate from initialization.

## Force/torque definitions

Using quadrature weights `w_l`:

- `F = sum_l F_l w_l`
- `M = sum_l (X_l - X_c) x F_l w_l`

`F` and `M` are reduced with `MPI_ALLREDUCE` even for `np=1`.

## Dynamics update (minimal explicit)

- `Uc^{n+1} = Uc^n + dt * F / m`
- `Xc^{n+1} = Xc^n + dt * Uc^{n+1}`
- `omega^{n+1} = omega^n + dt * M_perp / I_perp`
- `p^{n+1} = normalize(p^n + dt * (omega^{n+1} x p^n))`

`M_perp` is the torque component perpendicular to `p`.
This is a minimal axisymmetric-rod update for centerline motion; section spin about the rod axis is not modeled in 2.4.

## Input parameters (FiberParam)

- `rigid_free_test_active` (default `F`)
- `rigid_free_case` (default `1`)
- `fiber_mass` (default `1.0`)
- `fiber_inertia_perp` (default `1.0`)
- `free_output_interval` (default `1`)
- optional initial conditions: `fiber_uc`, `fiber_omega`

## Outputs

- `fiber_rigid_free_points.dat`
- `fiber_rigid_free_summary.dat`

Summary includes `xc/uc/p/omega`, total force/torque, slip metrics, spacing error, `p_norm_error`, and failure status.

## Run examples

- free case 1: `baseline/channel_re180_single_phase/input_fiber_rigid_free_case1_test.i3d`
- free case 2: `baseline/channel_re180_single_phase/input_fiber_rigid_free_case2_test.i3d`

## Boundaries of stage 2.4

Not implemented here:

- flexible beam dynamics
- tension constraints
- free-end beam boundary conditions
- collisions or wall contact
- multi-fiber interactions
