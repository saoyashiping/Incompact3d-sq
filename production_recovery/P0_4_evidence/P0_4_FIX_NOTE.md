# P0.4 Fix Note

This fix resolves the P0.4 validation failure caused by an incomplete velocity-bridge integration.

## Fixed blockers

1. `src/xcompact3d.f90` now imports `fibre_prod_velocity_bridge` explicitly and declares/initializes the velocity-sampling guard variables used by the existing one-way call path.
2. `src/CMakeLists.txt` now includes `fibre_prod_velocity_bridge.f90` in the `xcompact3d` target.
3. `src/CMakeLists.txt` now defines the missing `fibre_prod_velocity_bridge_check` target.

## Scope

- No pressure/projection/RK3/channel-forcing logic was changed.
- No structure advance was added.
- No reaction-force or force-buffer-to-RHS call was added to the velocity bridge.
- P0.4 remains diagnostic-only and one-way.

## Local direct checks

The following direct micro-checks were compiled and run with `gfortran` in the repair environment:

- `fibre_prod_main_hook_check`: PASS
- `fibre_prod_force_buffer_rhs_path_check`: PASS
- `fibre_prod_runtime_bridge_check`: PASS
- `fibre_prod_velocity_bridge_check`: PASS

Full validation still must be run on the target Ubuntu environment with the real 2decomp-fft install.
