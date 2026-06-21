# Production Recovery R9 Plan — Fibre-Fibre Collision Production Safety

## R9 scope

R9 adds standalone production fibre-fibre collision safety support.  The scope is limited to collision geometry, near-contact warning, overlap/penetration detection, deterministic pair metadata, penalty-force candidates, equal-and-opposite force consistency, damping-direction safety, and dry-structure response checks.

## R9 non-goals

R9 does not modify `src/xcompact3d.f90`, does not add an Xcompact3D main-loop hook, does not write a real Navier-Stokes RHS, does not modify wall-contact logic, does not run restart/statistics/visualization validation, does not run np=1/2/4 MPI production consistency, does not claim production DNS-FSI closure, and does not enter R10.

## Files created/modified

Created source files:

- `src/fibre_prod_collision_geometry.f90`
- `src/fibre_prod_fibre_collision.f90`
- `src/fibre_prod_collision_diagnostics.f90`
- `src/fibre_prod_fibre_collision_check.f90`

Created evidence files:

- `production_recovery/R9_PLAN.md`
- `production_recovery/R9_BUILD_LOG.txt`
- `production_recovery/R9_RUN_LOG.txt`
- `production_recovery/R9_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R9_PASS_FAIL.md`
- `production_recovery/R9_evidence/README.md`

Modified build file:

- `src/CMakeLists.txt` only to add the standalone `fibre_prod_fibre_collision_check` target.

## Build strategy

Configure a dedicated R9 build directory and build only the standalone `fibre_prod_fibre_collision_check` target.  The target compiles collision modules together with existing standalone fibre state, diagnostics, boundary, bending, tension, and structure solver modules.  It is not linked into `xcompact3d`.

## Run strategy

Run only the standalone check executable.  The check constructs safe separated, near-contact, node-node overlap, segment-segment overlap, closing/separating velocity, dry-structure response, invalid parameter, NaN coordinate, and degenerate-geometry cases without invoking DNS or MPI production validation.

## Pass/fail criteria

R9 PASS requires a successful build and a real run log containing `R9_FIBRE_PROD_FIBRE_COLLISION_CHECK PASS`.

R9 FAIL applies when the standalone target builds/runs but any required collision safety assertion fails.

R9 BLOCKED applies when configure, build, or execution cannot complete because required compiler, MPI wrapper, dependency path, or executable is unavailable.

## Evidence boundary

R9 evidence is limited to standalone fibre-fibre collision safety and force-candidate validation.  R9 PASS would not mean xcompact3d main-loop production coupling, real DNS RHS coupling, restart/statistics/visualization validation, np=1/2/4 MPI production consistency, or production DNS-FSI final closure.
