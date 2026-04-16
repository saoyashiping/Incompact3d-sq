# Stage 2.2 supplement: minimal Lagrangian-to-Euler spreading conservation test

## What is implemented
- A minimal spreading module (`src/fiber_spread.f90`) for **operator-only** conservation checks.
- Static single-fiber quadrature weights (trapezoidal rule): endpoint `ds/2`, interior `ds`.
- Minimal test Lagrangian force field `fiber_test_force(3, fiber_nl)`.
- Independent spreading test mode (`spread_test_active`) with summary outputs.

## What is intentionally not implemented
- No IBM forcing in momentum equations.
- No no-slip closure.
- No rigid/flexible dynamics.
- No collision/wall-contact/multi-fiber.
- No feedback to real solver fields.

## Discrete form used
For each Lagrangian point `l`:

`f(i,j,k) += F_l * w_l * delta_h(x_i-X_l, y_j-Y_l, z_k-Z_l)`

where:
- `F_l` is test force,
- `w_l` is line quadrature weight,
- `delta_h` is regularized kernel.

Euler total is computed by volume integration:

`sum_{i,j,k} f(i,j,k) * dV`

Conservation is verified against:

`sum_l F_l * w_l`.

## Grid/periodicity assumptions in this stage
- Test spreading grid is uniform in x/y/z (operator-level test grid).
- x/z distances use minimum-image periodic handling.
- y is non-periodic.

## How to run
```bash
./xcompact3d baseline/channel_re180_single_phase/input_fiber_spread_inner_test.i3d
./xcompact3d baseline/channel_re180_single_phase/input_fiber_spread_near_x0_test.i3d
./xcompact3d baseline/channel_re180_single_phase/input_fiber_spread_near_z0_test.i3d
```

## Outputs
- `fiber_spread_lagrangian.dat`: Lagrangian points, test force, quadrature weights.
- `fiber_spread_summary.dat`: Lagrangian totals, Euler totals, abs/rel errors, spread sumw min/max.

## Why near-x0 / near-z0 tests are included
These cases verify that x/z periodic minimum-image handling preserves conservation near periodic boundaries and avoids boundary-cutoff loss.
