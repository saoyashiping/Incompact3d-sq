# P0.4 Runtime Velocity Interpolation Bridge Plan

## Scope

P0.4 adds one-way diagnostic velocity sampling from runtime velocity fields to Lagrangian fibre sample points using the production IBM interpolation kernel.

## Non-goals

- No long production DNS.
- No fibre structure advance.
- No structure-to-fluid reaction force.
- No force-buffer-to-RHS feedback as a result of velocity sampling.
- No two-way coupling or production-ready claim.

## Objectives

1. Add a velocity sample container and allocation/finalization APIs.
2. Use production IBM interpolation to sample ux/uy/uz at Lagrangian points.
3. Preserve ux/uy/uz no-contamination.
4. Preserve RHS no-contamination.
5. Keep lambda=0 runtime bridge no-contamination intact.
6. Keep lambda>0 without feedback force from modifying RHS.
7. Preserve P0.2 and P0.3 validation paths.

## Next stage

P0.5 must initialize production fibre state and attach one-way sampled velocities. Production-run status remains blocked.
