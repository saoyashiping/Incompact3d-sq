# P0.0 Production-Readiness Hard Audit Plan

## Scope

P0.0 is a source/evidence hard audit of the uploaded code state. It does not repair source, modify physics, or run long DNS/FSI production cases.

## Restrictions honored

- Do not modify `src/*.f90`.
- Do not modify `src/CMakeLists.txt`.
- Do not modify `stage10_checks` through `stage22_checks`.
- Do not delete, overwrite, or rewrite existing `production_recovery/R0` through `production_recovery/R12` files.
- Add only P0.0 audit files under `production_recovery/` and `production_recovery/P0_0_evidence/`.

## Audit steps

1. Record R10/R11/R12 PASS/BLOCKED/FAIL status evidence.
2. Inspect `xcompact3d` target wiring in `src/CMakeLists.txt`.
3. Inspect `src/fibre_prod_rhs_adapter.f90` small-lambda RHS behavior.
4. Document the boundary between standalone helper PASS and production DNS-FSI main-program readiness.
5. Produce a P0.0 PASS/FAIL decision where PASS means the blocker is identified, not that production DNS-FSI is ready.
