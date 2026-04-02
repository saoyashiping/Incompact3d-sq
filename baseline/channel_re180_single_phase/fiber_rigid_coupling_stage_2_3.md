# Stage 2.3: Prescribed rigid-fiber no-slip coupling test

## Scope of this stage
This stage implements only prescribed-motion rigid-fiber no-slip coupling validation:

- Single rigid fiber only.
- Motion is prescribed analytically (not solved dynamically).
- Coupling path: interpolate Eulerian velocity to fiber points, compute slip, build penalty force, spread back to Eulerian grid, and add to momentum RHS only in this test mode.

This is **not** free rigid-body dynamics, not flexible fiber dynamics, and not multi-fiber coupling.

## Prescribed motion cases
- `rigid_motion_case = 1`: stationary rigid fiber (`Xdot = 0`).
- `rigid_motion_case = 2`: constant translation (`Xdot = rigid_translation_velocity`).

Fiber geometry is initialized once by the existing static-fiber initializer. Rigid motion updates all points using a single prescribed mapping, preserving point spacing.

## No-slip coupling model
At each time step:

1. Interpolate fluid velocity to fiber points: `U`.
2. Evaluate prescribed fiber-point velocity: `Xdot`.
3. Compute slip: `slip = U - Xdot`.
4. Compute coupling force at fiber points:

   `Ffs = ibm_beta * (U - Xdot)`

5. Spread `Ffs` to Eulerian body-force fields.
6. Add Eulerian body force to momentum RHS only when `rigid_coupling_test_active = T`.

## How to run
### Stationary test
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_rigid_stationary_test.i3d
```

### Translation test
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_rigid_translation_test.i3d
```

## Outputs to inspect
- `fiber_rigid_coupling_points.dat`
  - per point: position, `Xdot`, interpolated `U`, slip, and coupling force.
- `fiber_rigid_coupling_summary.dat`
  - per output step: `slip_max`, `slip_rms`, Lagrangian total force, Eulerian total force, force-balance error, `spacing_error_max`.

## Current limits
- No free rigid-body translation/rotation solve.
- No flexible beam/tension/bending model.
- No contact/collision/wall-contact model.
- No multi-fiber support.
- Solver-coupled rigid test path currently follows the existing single-process interpolation/spreading assumptions.

## Solver-grid spreading fix (2.3 consistency)
The original 2.3 coupling path mixed:
- interpolation on real solver grid (`xp/yp/zp`), and
- spreading on a synthetic uniform grid (`xg/yg/zg`).

This is now fixed for the rigid-coupling path:
- rigid-coupling spreading uses real solver coordinates `xp`, `yp`, `zp`;
- y-direction uses local stretched-grid spacing `hy_loc` (same local-dy style as solver-readonly interpolation);
- x/z keep minimum-image periodic distance handling.

The conservation test path still uses the synthetic uniform test grid by default, and remains separated from rigid solver-coupling spreading.

## Stabilization and diagnostics update
In short tests, large `|ibm_beta|` can produce very large coupling force on the first coupled step.
This may drive immediate instability before diagnostics are written.

To improve diagnosability without changing the 2.3 model:
- first coupled step now always writes rigid diagnostics (`itime == ifirst` forced output);
- then regular output cadence follows `rigid_output_interval`;
- summary/log now includes max-norm diagnostics:
  - `u_interp_max_norm`
  - `xdot_max_norm`
  - `coupling_force_max_norm`
  - `euler_force_max_norm`

Recommended beta sensitivity sequence (independent runs):
- `ibm_beta = -1`
- `ibm_beta = -10`
- `ibm_beta = -20`
- `ibm_beta = -30`
- `ibm_beta = -40`
- `ibm_beta = -50`
- `ibm_beta = -100`

If the no-ramp run fails early, retry with:
- `coupling_ramp_steps = 5`
- `coupling_ramp_steps = 10`

This beta scan is a test-stability study tool for stage 2.3, not a final physical parameter conclusion.

## Fail-fast and non-finite checks
The rigid-coupling path now performs finite-value checks for:
- `fiber_slip`
- `fiber_coupling_force`
- `fiber_euler_force_x/y/z`
- summary-level diagnostics (`slip_max`, `slip_rms`, `lag/euler totals`, force-balance error, spacing error, norm diagnostics)

If any non-finite value is detected:
- diagnostics are forced to output on that step;
- log prints `RIGID COUPLING TEST FAILED` with step/time/motion case and first failed quantity;
- run exits through `MPI_ABORT` (fail-fast), avoiding false-success completion.

## Effective beta logging
- `beta_input` is always the user value (`ibm_beta`).
- `beta_eff` and `ramp_factor` are logged explicitly.
- `coupling_ramp_steps` is a test-only stabilization knob (default `0`, disabled).
- For `coupling_ramp_steps > 0`, the code uses a linear startup ramp with coupling-step index
  `n = itime - ifirst + 1` (clamped at least to 1):
  `ramp_factor = min(1, n / coupling_ramp_steps)`, so the first coupling step is `1/coupling_ramp_steps`.
  Then `beta_eff = ramp_factor * beta_input`.
- For `coupling_ramp_steps = 0`, `ramp_factor = 1` and `beta_eff = beta_input`.

## Failure code map
- `0`: no failure
- `1`: nonfinite slip
- `2`: nonfinite coupling_force
- `3`: nonfinite euler_force
- `4`: nonfinite summary quantity
- `5`: other numerical blowup (reserved)

## Interpreting outcomes
- **stable and usable**: no failure flag, finite diagnostics throughout target window.
- **diagnosable but unstable**: finite first-step diagnostics + later fail-fast with failure code/logs.
- **immediate failure**: fail-fast at first coupled step; use smaller `|beta|` and/or ramp.
