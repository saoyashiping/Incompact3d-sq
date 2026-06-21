# P0.5 Production Fibre State Velocity Attachment Plan

## Scope

P0.5 initializes a minimal production fibre state, samples one-way fluid velocity at fibre node coordinates, and attaches sampled velocity to diagnostic state storage.

## Non-goals

- No long production DNS.
- No structure advance.
- No structure-to-fluid reaction force.
- No RHS feedback from sampled velocity.
- No two-way coupling or production-ready claim.

## Objectives

1. Add sampled velocity storage to production fibre state.
2. Add state velocity attachment APIs.
3. Use fibre node coordinates as Lagrangian velocity sampling points.
4. Attach sampled fluid velocity to state diagnostic storage only.
5. Preserve fibre coordinates, state velocity/acceleration, force buffer, velocity fields, and RHS arrays.
6. Preserve P0.2/P0.3/P0.4 validation paths.

## Next stage

P0.6 must add sampled-velocity-to-fluid-force-input candidate diagnostics with no RHS feedback. Production-run status remains blocked.
