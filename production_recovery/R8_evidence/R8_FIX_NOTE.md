# R8 Fix Note — Tension Pure-Function Status Output

## Cause

GNU Fortran rejects `pure function` procedures that modify an `integer, intent(out) :: status` argument.  The R8 wall-contact standalone target depends on `src/fibre_prod_tension_solver.f90`, where `fibre_prod_tension_segment_length_residual` and `fibre_prod_tension_max_stretch_error` previously had the `pure` attribute while preserving status-output semantics.

## Fix

The fix removes the `pure` attribute from exactly these two functions:

- `fibre_prod_tension_segment_length_residual`
- `fibre_prod_tension_max_stretch_error`

The function names, arguments, return values, numeric logic, and status semantics are unchanged.  They were not converted to subroutines and no caller interface was changed.

## Explicit non-actions

- `src/xcompact3d.f90` was not modified.
- No R8 code was connected to the xcompact3d main loop.
- No fibre-fibre collision was implemented.
- No real Navier-Stokes RHS write or coupling was implemented.
- No R9 work was performed.
