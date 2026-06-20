# R8 Fix Note

## Failure

R8 failed while compiling `src/fibre_prod_tension_solver.f90` because two functions were declared as `pure` while also using an `integer, intent(out) :: status` argument.

GNU Fortran 15.x correctly rejects this pattern: a `pure function` cannot define an `intent(out)` dummy argument or pass it to another procedure in a variable-definition context.

Affected routines:

- `fibre_prod_tension_segment_length_residual`
- `fibre_prod_tension_max_stretch_error`

## Fix

Removed only the `pure` attribute from the two affected functions.

The following were kept unchanged:

- function names;
- argument lists;
- return values;
- status-code semantics;
- segment-length residual logic;
- max-stretch-error logic.

## Verification performed in this environment

A standalone direct `gfortran` compile/link/run was performed for the R8 wall-contact check dependency chain:

- `fibre_prod_config.f90`
- `fibre_prod_state.f90`
- `fibre_prod_diagnostics.f90`
- `fibre_prod_boundary_conditions.f90`
- `fibre_prod_bending_solver.f90`
- `fibre_prod_tension_solver.f90`
- `fibre_prod_structure_solver.f90`
- `fibre_prod_wall_geometry.f90`
- `fibre_prod_wall_contact.f90`
- `fibre_prod_wall_contact_diagnostics.f90`
- `fibre_prod_wall_contact_check.f90`

The direct run printed:

```text
R8_FIBRE_PROD_WALL_CONTACT_CHECK PASS
```

## Boundary

- `src/xcompact3d.f90` was not modified.
- R8 was not connected to the Xcompact3D main loop.
- No real Navier-Stokes RHS coupling was implemented.
- No fibre-fibre collision implementation was added.
- `stage20_checks/`, `stage21_checks/`, and `stage22_checks/` were not modified.
