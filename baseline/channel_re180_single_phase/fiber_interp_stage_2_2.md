# Stage 2.2: Minimal Euler-to-Lagrangian interpolation operator

## What is implemented
1. A single regularized delta kernel (4-point compact support) in `src/fiber_delta.f90`.
2. Euler-to-Lagrangian velocity interpolation in `src/fiber_interp.f90`.
3. Minimal interpolation-result storage in fiber state (`fiber_uinterp`, `fiber_uexact`, `fiber_uerror`, `fiber_interp_max_error`).
4. One-shot interpolation output file `fiber_interp.dat` with interpolated values, exact values, errors and `sumw`.
5. Controlled operator test mode with two cases, and explicit `sumw_min/sumw_max` logging.
   - `interp_test_case = 1`: constant velocity field
   - `interp_test_case = 2`: linear velocity field


## Periodic-direction handling (2.2 fix)
- Interpolation distance in **x** and **z** uses minimum-image mapping with periods `xlx` and `zlz`.
- **y** direction remains non-periodic (no wrapping).
- This fix is local to interpolation distance evaluation only; fiber geometry initialization is unchanged.

## What is intentionally not implemented
- No Lagrangian-to-Euler spreading.
- No IBM forcing / no-slip penalty / coupling terms.
- No rigid-body or flexible-body dynamics.
- No collision or wall-contact model.
- No multi-fiber support.

Fiber remains a static geometry object and does not feed back to the flow solver.

When `interp_test_active = T`, the code now exits immediately after writing interpolation-test outputs and does not enter solver initialization/time integration.

## Files added
- `src/fiber_delta.f90`
- `src/fiber_interp.f90`
- `baseline/channel_re180_single_phase/input_fiber_interp_constant_test.i3d`
- `baseline/channel_re180_single_phase/input_fiber_interp_linear_test.i3d`
- `baseline/channel_re180_single_phase/fiber_interp_stage_2_2.md`

## Existing files minimally changed
- `src/fiber_types.f90`
- `src/fiber_io.f90`
- `src/parameters.f90`
- `src/xcompact3d.f90`
- `src/CMakeLists.txt`

## How to run controlled interpolation tests
Constant-field test:
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_interp_constant_test.i3d
```

Linear-field test:
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_interp_linear_test.i3d
```

Constant-field near x-periodic boundary:
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_interp_constant_near_x0_test.i3d
```

Constant-field near z-periodic boundary:
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_interp_constant_near_z0_test.i3d
```

## How to check outputs
- `fiber_points.dat`: initialized Lagrangian points.
- `fiber_interp.dat`: columns
  - index, x, y, z
  - u_interp, v_interp, w_interp
  - u_exact, v_exact, w_exact
  - err_u, err_v, err_w
  - sumw
- Also check console `max abs error`, `sumw_min`, and `sumw_max` printed for each test case.

## Test meaning
- Constant-field test verifies basic partition-of-unity behavior of interpolation weights.
- Linear-field test verifies operator reasonableness against an analytic target field at fiber points.
