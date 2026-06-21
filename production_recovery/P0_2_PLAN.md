# P0.2 Physical Force-Buffer-to-RHS Integration Plan

## Scope

P0.2 verifies the controlled micro-path from an IBM-generated `fibre_prod_force_buffer_type` Eulerian force-density buffer through the production main hook into DNS RHS arrays.

## Non-goals

- P0.2 does not run long DNS.
- P0.2 does not modify pressure/projection/RK3/channel-forcing logic.
- P0.2 does not claim production DNS-FSI readiness.
- P0.2 does not bridge runtime xcompact3d grid/state into the fibre production state.

## Objectives

1. Add a buffer-level public main-hook API for `fibre_prod_force_buffer_type`.
2. Preserve P0.1 lambda=0 no-op and lambda>0 missing-buffer blocking behavior.
3. Verify IBM spreading creates a nonzero finite force-density buffer.
4. Verify force-buffer total force is conserved by volume integration.
5. Verify RHS increments equal `lambda_fsi * penalty_beta * force_buffer%fx/fy/fz`.
6. Verify the RHS increment is not a uniform whole-domain contribution.

## Next stage

P0.3 must implement the runtime grid/state bridge and lambda=0 no-contamination in the `xcompact3d` call path. Production-run status remains blocked until later stages close the full DNS-FSI loop.
