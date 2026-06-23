# P1_4 build-target/link fix

Reason:
- The P1_4 validation stopped during the CMake build step and reported:
  `Reason: command failed; see .../P1_4_BUILD_LOG.txt`.
- Code inspection showed that the newly added target
  `fibre_prod_p1_np_consistency_closure_case_check`
  was declared in `src/CMakeLists.txt` but was not linked against `decomp2d` or `MPI::MPI_Fortran`.
- The target compiles modules that use MPI/decomp-linked infrastructure, especially
  `fibre_prod_p1_np_consistency_closure_case.f90`, so the missing target link can cause the P1_4 build step to fail before any real DNS-FSI run is reached.

Fix:
- Add `target_link_libraries(fibre_prod_p1_np_consistency_closure_case_check PRIVATE decomp2d)`.
- Add MPI link guard for the same target.
- Add install rule for the target.
- Improve P1_4 failure reporting by writing `P1_4_LAST_FAILED_LOG_TAIL.txt` and by printing build-log tails from `run_P1_4.sh`.
- Mark unfinished P1_4 evidence audits as `Result: SKIPPED` instead of leaving stale `Result: PENDING` when the validation fails before those audits are reached.

No physics change:
- No pressure/projection/RK3/channel-forcing logic was modified.
- No P1_0/P1_1/P1_2/P1_3 pass/fail file is read by P1_4.
- P1_4 remains self-contained.
