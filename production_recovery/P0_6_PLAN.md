# P0.6 Plan: Sampled-Velocity-to-Hydrodynamic-Input Candidate Diagnostics

P0.6 converts one-way sampled fluid velocity already attached to the production fibre state into a structure-side hydrodynamic input candidate diagnostic.

## Scope

- Read sampled fluid velocity from production fibre state storage.
- Compare it with deterministic structure node velocity to form relative velocity.
- Compute `candidate_force = beta_hydro * relative_u`.
- Attach the candidate only to production fibre state diagnostic storage.
- Preserve P0.2/P0.3/P0.4/P0.5 safety checks.

## Explicit non-goals

- No production DNS run.
- No two-way coupling.
- No structure advance.
- No Eulerian force-buffer write.
- No RHS adapter call.
- No main-hook force-buffer-to-RHS call.
- No production-ready claim.

## Next stage

P0.7 is required for controlled structure-input candidate handoff, still with no structure advance and no RHS feedback.
