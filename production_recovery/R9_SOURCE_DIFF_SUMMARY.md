# Production Recovery R9 Source Diff Summary

## New source files

- `src/fibre_prod_collision_geometry.f90` — finite point/segment distance, segment-segment clearance, effective gap, and geometry validation utilities.
- `src/fibre_prod_fibre_collision.f90` — collision state, deterministic pair metadata, near-contact and overlap detection, and equal-and-opposite force-candidate routines.
- `src/fibre_prod_collision_diagnostics.f90` — min-gap, max-penetration, force norm, action-reaction residual, energy candidate, finite diagnostics, and summary utilities.
- `src/fibre_prod_fibre_collision_check.f90` — standalone R9 fibre-fibre collision safety check driver.

## New evidence/documentation files

- `production_recovery/R9_PLAN.md`
- `production_recovery/R9_BUILD_LOG.txt`
- `production_recovery/R9_RUN_LOG.txt`
- `production_recovery/R9_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R9_PASS_FAIL.md`
- `production_recovery/R9_evidence/README.md`

## Modified files

- `src/CMakeLists.txt` — adds only the standalone `fibre_prod_fibre_collision_check` target.

## `src/xcompact3d.f90` modified?

No.  `src/xcompact3d.f90` was not modified in R9.

## `src/CMakeLists.txt` modified?

Yes.  The modification is limited to the standalone `fibre_prod_fibre_collision_check` target.

## Connected to `xcompact3d` executable?

No.  R9 modules are not added to the `xcompact3d` executable source list and no main-loop hook is added.

## Real RHS coupling implemented?

No.  R9 does not write or couple to a real Navier-Stokes RHS.

## Wall-contact logic modified?

No.  R9 does not modify wall-contact logic.

## Production DNS-FSI validation executed?

No.  R9 does not execute production DNS-FSI validation or np=1/2/4 MPI consistency checks.
