# P0.1 Plan: Production Module Build-Chain Integration and RHS Safety Gate

## Purpose

P0.1 repairs the blocker identified by P0.0: the `xcompact3d` target was wired only to a weak production hook path, and the active nonzero-lambda RHS adapter used a whole-field uniform scalar perturbation.

## Scope

P0.1 makes only build-chain and RHS-safety changes:

1. Add the complete production fibre helper chain to the `xcompact3d` source list.
2. Add the same chain to `fibre_prod_main_hook_check` so the check is compiled in the same production-module context.
3. Remove the fallback whole-domain uniform RHS perturbation.
4. Require an explicit Eulerian force-density buffer for any nonzero-lambda RHS modification.
5. Preserve lambda=0 and disabled-mode no-contamination behavior.

## Non-goals

P0.1 does not claim full paper-production DNS-FSI readiness. It does not yet construct the force-density buffer from runtime Lagrangian fibres inside `xcompact3d`. That is the next P0 stage.

## Pass meaning

P0.1 PASS means the unsafe uniform RHS smoke path has been removed and the full production fibre module chain is visible to the main build target. It does not mean long-time DNS statistics, MPI production FSI, or paper-level validation are complete.
