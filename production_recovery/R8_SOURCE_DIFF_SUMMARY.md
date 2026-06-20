# Production Recovery R8 Source Diff Summary

## New source files

- `src/fibre_prod_wall_geometry.f90` — channel wall geometry, validation, point and segment signed-gap utilities.
- `src/fibre_prod_wall_contact.f90` — near-wall warning, penetration detection, and wall-contact penalty-force candidate routines.
- `src/fibre_prod_wall_contact_diagnostics.f90` — wall-contact diagnostic counts, norms, energy candidate, finite checks, and status summary.
- `src/fibre_prod_wall_contact_check.f90` — standalone R8 wall-contact safety check driver.

## New evidence/documentation files

- `production_recovery/R8_PLAN.md`
- `production_recovery/R8_BUILD_LOG.txt`
- `production_recovery/R8_RUN_LOG.txt`
- `production_recovery/R8_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R8_PASS_FAIL.md`
- `production_recovery/R8_evidence/README.md`

## Modified files

- `src/CMakeLists.txt` — adds only the standalone `fibre_prod_wall_contact_check` target.

## `src/xcompact3d.f90` modified?

No.  `src/xcompact3d.f90` was not modified in R8.

## `src/CMakeLists.txt` modified?

Yes.  The modification is limited to the standalone `fibre_prod_wall_contact_check` target.

## Connected to `xcompact3d` executable?

No.  R8 modules are not added to the `xcompact3d` executable source list and no main-loop hook is added.

## Real RHS coupling implemented?

No.  R8 does not write or couple to a real Navier-Stokes RHS.

## Fibre-fibre collision implemented?

No.  R8 does not implement fibre-fibre collision.

## Production DNS-FSI validation executed?

No.  R8 does not execute production DNS-FSI validation or np=1/2/4 MPI consistency checks.
