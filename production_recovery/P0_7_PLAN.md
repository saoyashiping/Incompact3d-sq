# P0.7 Plan: Controlled Structure-Input Candidate Handoff

P0.7 performs a diagnostic-only handoff of the P0.6 hydrodynamic input candidate into a controlled structure-side input slot.

## Scope

- Read `hydro_force_candidate` from production fibre state storage.
- Copy it to `structure_input_force` handoff storage.
- Attach the handoff buffer back to production fibre state structure-input storage.
- Preserve P0.2/P0.3/P0.4/P0.5/P0.6 safety checks.

## Explicit non-goals

- No production DNS run.
- No two-way coupling.
- No structure advance.
- No structure state time integration.
- No Eulerian force-buffer write.
- No RHS adapter call.
- No main-hook force-buffer-to-RHS call.
- No production-ready claim.

## Next stage

P0.8 is required for controlled structure dry-step preflight, still with no fluid feedback.
