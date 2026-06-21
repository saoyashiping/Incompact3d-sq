# P0.8 Plan: Controlled Structure Dry-Step Preflight

P0.8 consumes the P0.7 `structure_input_force` in a controlled dry-step predictor that writes only scratch/trial storage.

## Scope

- Read structure coordinates, structure velocity or zero, and structure input force.
- Compute scratch `dx_trial`, `x_trial`, and `u_trial` using the dry-step predictor.
- Check trial response finiteness and boundedness.
- Preserve P0.2/P0.3/P0.4/P0.5/P0.6/P0.7 safety checks.

## Explicit non-goals

- No production DNS run.
- No two-way coupling.
- No committed structure state write-back.
- No structure-to-fluid reaction force.
- No Eulerian force-buffer write.
- No RHS adapter call.
- No main-hook force-buffer-to-RHS call.
- No production-ready claim.

## Next stage

P0.9 is required for controlled structure dry-step commit gate, still with no fluid feedback.
