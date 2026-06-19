# Production Recovery R8 Plan — Wall Contact Production Safety

## R8 scope

R8 adds a standalone production wall-contact safety layer for fibre/channel-wall interactions.  It covers wall geometry, signed gap diagnostics, near-wall warnings, penetration detection, wall-contact penalty-force candidates, and a dry-structure safety check using the existing standalone structure solver.

## R8 non-goals

R8 does not modify `src/xcompact3d.f90`, does not add an Xcompact3D main-loop hook, does not write the Navier-Stokes RHS, does not implement fibre-fibre collision, does not validate restart/statistics/visualization, does not run np=1/2/4 production MPI consistency, and does not claim production DNS-FSI final closure.

## Files created/modified

Created source files:

- `src/fibre_prod_wall_geometry.f90`
- `src/fibre_prod_wall_contact.f90`
- `src/fibre_prod_wall_contact_diagnostics.f90`
- `src/fibre_prod_wall_contact_check.f90`

Created evidence files:

- `production_recovery/R8_PLAN.md`
- `production_recovery/R8_BUILD_LOG.txt`
- `production_recovery/R8_RUN_LOG.txt`
- `production_recovery/R8_SOURCE_DIFF_SUMMARY.md`
- `production_recovery/R8_PASS_FAIL.md`
- `production_recovery/R8_evidence/README.md`

Modified build file:

- `src/CMakeLists.txt` only to add the standalone `fibre_prod_wall_contact_check` target.

## Build strategy

Configure a dedicated R8 build directory and build only the standalone `fibre_prod_wall_contact_check` target.  The R8 target compiles wall-contact modules with the existing standalone fibre state, diagnostics, boundary, bending, tension, and structure solver modules.  It is not linked into the `xcompact3d` executable.

## Run strategy

Run the standalone check executable only.  The check constructs safe, near-wall, lower-penetration, upper-penetration, damping-direction, and fail-closed cases without invoking a DNS run or MPI production validation.

## Pass/fail criteria

R8 PASS requires a successful build and a real run log containing `R8_FIBRE_PROD_WALL_CONTACT_CHECK PASS`.

R8 FAIL applies when the standalone target builds/runs but the check fails or omits required wall-contact safety assertions.

R8 BLOCKED applies when the environment prevents configure, build, or execution, or when the required executable is not produced.

## Evidence boundary

R8 evidence is limited to standalone wall-contact safety and force-candidate validation.  R8 PASS would not mean xcompact3d main-loop coupling, real DNS RHS coupling, fibre-fibre collision, restart/statistics/visualization, np=1/2/4 MPI production consistency, or production DNS-FSI final closure.
