# R3 Source Diff Summary

## Root-cause repair

The R3 failure was caused by a CMake target registration omission.

The files `src/fibre_prod_grid_adapter.f90` and `src/fibre_prod_grid_adapter_check.f90` existed, but `src/CMakeLists.txt` did not contain:

```cmake
add_executable(fibre_prod_grid_adapter_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_grid_adapter_check.f90)
```

As a result, this command failed:

```bash
cmake --build build_r3_grid --target fibre_prod_grid_adapter_check -j 2
```

with:

```text
gmake: *** No rule to make target 'fibre_prod_grid_adapter_check'. Stop.
```

## Files added by R3

```text
src/fibre_prod_grid_adapter.f90
src/fibre_prod_grid_adapter_check.f90
production_recovery/R3_PLAN.md
production_recovery/R3_BUILD_LOG.txt
production_recovery/R3_RUN_LOG.txt
production_recovery/R3_SOURCE_DIFF_SUMMARY.md
production_recovery/R3_PASS_FAIL.md
production_recovery/R3_evidence/README.md
```

## Files modified by this repair

```text
src/CMakeLists.txt
production_recovery/R3_RUN_LOG.txt
production_recovery/R3_PASS_FAIL.md
production_recovery/R3_SOURCE_DIFF_SUMMARY.md
production_recovery/PRODUCTION_RECOVERY_R3_CLOSED.md
```

## xcompact3d main-loop status

```text
src/xcompact3d.f90 was not modified.
No fibre_prod_* module is connected to the production time loop.
```

## CMake status

`src/CMakeLists.txt` was modified only to add the standalone R3 target:

```text
fibre_prod_grid_adapter_check
```

The R3 module was not added to the `xcompact3d` executable.

## R3 implementation scope

R3 implements only the standalone production grid adapter foundation:

1. coordinate storage;
2. local range metadata;
3. periodic flags;
4. spacing calculation;
5. cell-volume calculation;
6. point-to-cell lookup;
7. validation and destroy/deallocation checks.

## Forbidden scope not implemented

R3 does not implement IBM interpolation.
R3 does not implement IBM spreading.
R3 does not implement RHS coupling.
R3 does not implement structure advancement.
R3 does not implement wall contact.
R3 does not implement fibre-fibre collision.
R3 does not modify Stage 20/21/22 source-only checks.

## Additional code-quality cleanup during repair

The standalone grid-adapter module was also checked with direct `gfortran -Wall -Wextra -fcheck=all` compilation. Two minor warnings were cleaned:

1. removed an unused `last_index` dummy argument from the internal spacing helper;
2. replaced real equality at the upper coordinate endpoint with an interval-safe `>=` endpoint check.

After cleanup, the direct standalone check still produced:

```text
R3_FIBRE_PROD_GRID_ADAPTER_CHECK PASS
```
