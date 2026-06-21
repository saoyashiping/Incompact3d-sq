# P0.9 Source Diff Summary

- Added `fibre_prod_structure_commit_gate` for finite/bounded/shape-gated trial-state evaluation and optional structure-side commit.
- Extended `fibre_prod_state` with committed structure velocity storage and helpers for structure trial commits.
- Added `fibre_prod_structure_commit_gate_check` and CMake target coverage.
- Added a default-off `xcompact3d` diagnostic call path guarded by `FIBRE_PROD_STRUCTURE_DRY_COMMIT_ENABLE`.
- Added P0.9 validation and audit evidence scaffolding.


## P0.9 validation portability repair

The P0.9 source implementation already contains the commit gate module, state commit helpers, xcompact3d import/call path, and CMake check target. The validation script, however, used `rg` without a fallback. On systems without ripgrep this caused false static-audit failures.

Changed file:
- `production_recovery/P0_9_evidence/P0_9_VALIDATION_COMMAND.sh`

Change summary:
- added `search_quiet_regex()` and `search_lines_regex()` wrappers;
- use `rg` if available, otherwise use POSIX `grep -E`;
- token checks use `grep -Fq`;
- no production Fortran source or DNS/RHS logic changed by this repair.

## P0.9 executable-path repair

The latest P0.9 validation reached the build stage successfully, but failed when running checks because the validation script assumed check executables were located under `build_p0_9/src/<target>`. In this build layout, CMake links them under `build_p0_9/bin/<target>`.

Changed file:
- `production_recovery/P0_9_evidence/P0_9_VALIDATION_COMMAND.sh`

Change summary:
- added `find_built_exe()` to resolve check executables from `build_p0_9/bin`, `build_p0_9/src`, or a fallback `find` search;
- changed the run loop to execute the resolved path and log it;
- no production Fortran source, physics, RHS coupling, pressure/projection/RK3/channel-forcing logic changed by this repair.
