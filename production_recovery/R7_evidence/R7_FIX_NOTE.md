# R7 Fix Note

## Problem

The R7 CMake build failed while compiling `src/fibre_prod_tension_solver.f90` because two functions were declared as `pure` while assigning an `intent(out)` status argument.

GNU Fortran rejects this pattern for pure functions:

```text
Error: Argument 'status' of pure function ... must be INTENT(IN) or VALUE
Error: Variable 'status' cannot appear in a variable definition context
```

## Fix

Removed the `pure` attribute from:

```text
fibre_prod_tension_segment_length_residual
fibre_prod_tension_max_stretch_error
```

The numerical logic and public API were otherwise preserved.

## Verification performed here

A direct standalone `gfortran` compilation/link/run of the R7 check source set was performed. The executable printed:

```text
R7_FIBRE_PROD_FSI_CLOSED_LOOP_CHECK PASS
```

## Required user-side verification

After applying this patch to the Ubuntu project, rerun the normal R7 CMake validation command to regenerate authoritative `R7_BUILD_LOG.txt`, `R7_RUN_LOG.txt`, and `R7_PASS_FAIL.md`.

## Boundary

This fix does not connect R7 modules to `xcompact3d.f90`, does not write real DNS RHS, does not implement wall contact, and does not implement fibre-fibre collision.
