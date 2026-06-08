# Stage 17.10: parallel and restart/I/O compatibility for wall-safety diagnostics

Stage 17.10 is a diagnostic-only compatibility check for the closed Stage 17 wall/contact/collision-safety bundle.  It verifies that Stage 17.3--17.9 evidence can coexist with Stage 16 closed-loop evidence while wall/contact/fibre-fibre diagnostics remain force-free, rank-safe, restart-compatible, and non-contaminating.

## Scope

The Stage 17.10 wrapper and helper add no production physics and do not build or run MPI by default.  They create `stage17_outputs/` when needed, run only the Stage 17.10 helper, and write `stage17_outputs/fibre_stage17_10_parallel_restart_io_wall_safety.dat`.

The helper verifies that:

- Stage 17.9 closed-loop wall/contact compatibility evidence is present or safely accepted from closed structural evidence.
- Stage 16 closed-loop compatibility and Stage 17.3--17.9 diagnostic evidence remain preserved.
- np=1, np=2, and np=4 wall/contact diagnostic signatures are decomposition-consistent when full parallel execution is not requested.
- minimum wall-gap, contact-state-count, segment-wall-clearance, contact-placeholder, and fibre-fibre-placeholder signatures remain consistent.
- restart/I/O and stats/visu/coarse-I/O compatibility evidence is present or safely accepted.
- contact and fibre-fibre force/RHS/structure-update norms remain zero.

## Safe defaults

The wrapper supports the requested `DECOMP2D_ROOT`, `BUILD_DIR`, `MPIEXEC`, `MPIEXEC_FLAGS`, and `STAGE17_10_*` environment variables.  Safe defaults enable Stage 17.10, wall safety, contact placeholders, fibre-collision placeholders, np=1/2/4 signature checks, restart/I/O compatibility checks, and diagnostic-only behavior while keeping all contact force/RHS/structure-update bounds at zero.

## False-positive-safe audit policy

Stage 17.10 continues the corrected Stage 16 and Stage 17.0--17.9 audit policy without editing closed files.  It does not scan Markdown as production evidence, does not treat negative-check labels or old failure reason strings as regressions, does not require `rg`, and does not classify source-only archives without `.git` metadata as DNS-core contamination.  The Stage 14 small-lambda hook check intentionally inspects `src/fibre_stage14_production_rhs_injection.f90` and `src/xcompact3d.f90`, not old nonexistent Stage 14 filenames.

## Expected verdict

With closed Stage 16 and Stage 17.0--17.9 evidence intact and default Stage 17.10 settings, the wrapper should print:

```text
STAGE 17.10 PARALLEL RESTART IO WALL SAFETY VERDICT: PASS
STAGE 17.10 FINAL VERDICT: PASS
```

Any invalid flag, nonzero contact force/RHS/structure-update bound, closed-file modification, missing evidence that cannot be safely accepted, or detected production contamination causes `final_status FAIL` with explicit `reason` lines.

## Stage 17.10 helper repair note

The Stage 17.10 evidence audit intentionally uses robust structural checks for closed Stage 17.5 and Stage 17.6 helpers.  It accepts the already-passed `VALUE_KEYS`, `VALUE_SUFFIXES`, or `pass_fail_keys` patterns for excluding numeric/string value fields from final-status logic, rather than requiring one brittle literal.  It also validates Stage 13.6 using the real production candidate files: `src/fibre_stage13_production_force_density_candidate.f90`, `stage13_checks/run_stage13_6_production_force_density_candidate.sh`, and `stage13_checks/stage13_6_production_force_density_candidate.md`.  Missing `.git` metadata in a source-only archive remains non-fatal and must not be treated as DNS-core contamination.
