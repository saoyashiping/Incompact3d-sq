# P0.0 Standalone Helper vs Production DNS-FSI Boundary Audit

## What standalone helper checks can establish

The `fibre_prod_*_check` executables can serve as module-level self-checks. They are useful for testing isolated production-module behavior such as state handling, grid adaptation, IBM interpolation/spreading helpers, structure solver behavior, closed-loop surrogate behavior, wall contact, fibre collision, and main-hook smoke behavior.

## What standalone helper checks cannot establish

Standalone helper PASS cannot be treated as equivalent to `xcompact3d` main-program DNS-FSI closed-loop PASS. A helper can validate an isolated module or surrogate path without proving that the same modules are compiled into, called by, synchronized with, and numerically coupled inside the production `xcompact3d` time-integration path.

## Current boundary finding

The current `xcompact3d` target directly includes only the production runtime config, main diagnostics, RHS adapter, and main hook. The complete production FSI module chain is not integrated into the main target. The active RHS adapter small-lambda path is a uniform RHS injection smoke test rather than a physical IBM/Lagrangian fibre reaction-force spreading path.

## Audit conclusion

The current uploaded source version cannot directly enter paper production DNS-FSI runs. The correct next step is production module build-chain integration with no-op/lambda0 safety, not production physics runs.
