# P1_4 source diff summary

## Fixed build-stage issues

1. `src/fibre_prod_p1_np_consistency_closure_case.f90`
   - Shortened P1_4 public procedure names from `fibre_prod_p1_np_consistency_closure_case_*` to `p1_4_np_*` where needed.
   - This removes Fortran identifier-length violations such as `record_structure_response` and `record_global_signatures` exceeding 63 characters.
   - Changed the MPI global-count reduction from scalar actual arguments to one-element array buffers.

2. `src/fibre_prod_p1_np_consistency_closure_case_check.f90`
   - Updated calls to the shortened `p1_4_np_*` API.

3. `src/xcompact3d.f90`
   - Updated `use ... only:` and runtime calls to the shortened `p1_4_np_*` API.
   - No changes to pressure/projection/RK3/channel-forcing logic.

4. `src/CMakeLists.txt`
   - Keeps the P1_4 check target linked with `decomp2d` and `MPI::MPI_Fortran`.

5. `production_recovery/P1_4_evidence/P1_4_VALIDATION_COMMAND.sh`
   - Still self-contained.
   - Does not read P1_0/P1_1/P1_2/P1_3 PASS files or logs as pass criteria.

## Production status

P1_4 must still be rerun. P1 is not closed until `production_recovery/P1_4_PASS_FAIL.md` reports `Result: PASS`.

Production-run status: STILL BLOCKED FOR PAPER-SCALE DNS
