# R5 Fix Note

The R5 build failure was caused by declaring two functions in `src/fibre_prod_tension_solver.f90` as `pure` while they assign an `intent(out)` status argument.

Fortran pure functions cannot define non-local state through `intent(out)`/`intent(inout)` dummy arguments in this way. GNU Fortran 15.2 correctly rejected the functions.

Fix applied:

- removed `pure` from `fibre_prod_tension_segment_length_residual`;
- removed `pure` from `fibre_prod_tension_max_stretch_error`.

The routines are still side-effect limited in practice, but are no longer declared `pure`, allowing their explicit status-return API to compile.

Standalone direct gfortran validation in this patch environment produced:

```text
R5_FIBRE_PROD_STRUCTURE_SOLVER_CHECK PASS
```

Official R5 status should still be confirmed by rerunning the Ubuntu CMake/mpif90 R5 validation command in the project environment.
