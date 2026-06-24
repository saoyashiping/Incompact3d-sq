# P1_4 identifier-length / MPI buffer fix

This fix addresses a build-stage P1_4 failure in `fibre_prod_p1_np_consistency_closure_case.f90`.

Root cause:
- Several public Fortran procedure names exceeded the standard 63-character identifier limit, e.g. `fibre_prod_p1_np_consistency_closure_case_record_structure_response` and `fibre_prod_p1_np_consistency_closure_case_record_global_signatures`.
- With gfortran this causes `Name ... is too long` / `Syntax error in PUBLIC statement`, so P1_4 fails before any real DNS-FSI run.
- The MPI global-count reduction also used scalar actual arguments with an explicit MPI interface. It is now converted to one-element arrays for safer `MPI_Allreduce` compatibility.

Scope:
- Renamed only the P1_4 internal API symbols to short `p1_4_np_*` names.
- Updated `xcompact3d.f90` and the P1_4 check program to use the new names.
- Kept the module filename, check target name, runtime tokens, physical logic, lambda gate, RHS formula, pressure/projection/RK3/channel forcing unchanged.
- P1_4 remains self-contained and does not read P1_0-P1_3 evidence as pass criteria.
