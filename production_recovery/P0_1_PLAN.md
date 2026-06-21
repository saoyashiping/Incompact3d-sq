# P0.1 Production Module Build-Chain Integration Plan

## Scope

P0.1 implements production module build-chain integration and lambda=0 no-op RHS safety gating for the main hook path.

## Non-goals

- P0.1 does not run long DNS cases.
- P0.1 does not claim production DNS-FSI readiness.
- P0.1 does not complete the physical FSI closed loop.
- P0.1 does not advance fibre state from the `xcompact3d` time integrator.

## Objectives

1. Add the complete production fibre module chain to the `xcompact3d` target.
2. Preserve strict lambda=0 no-op behavior.
3. Block nonzero-lambda RHS modification unless explicit Eulerian force-density buffers are supplied.
4. Remove the old small-lambda uniform RHS constant injection path.
5. Reserve the explicit force-buffer-to-RHS interface for P0.2 physical integration.

## Next stage

P0.2 is the stage that must connect physical force buffers into the production RHS path. Until then, production-run status remains blocked.
