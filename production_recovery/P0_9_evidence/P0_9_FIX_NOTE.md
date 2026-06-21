# P0.9 validation script repair note

The latest P0.9 run passed the static audits and built the requested targets, including `fibre_prod_main_hook_check` and `fibre_prod_structure_commit_gate_check`, but failed at the run phase because the validation script tried to execute checks from:

```text
build_p0_9/src/<target>
```

In the actual CMake layout, the linked executables are placed under:

```text
build_p0_9/bin/<target>
```

Fix:
- add `find_built_exe()` to resolve executable paths robustly;
- search in `build_p0_9/bin`, `build_p0_9/src`, and `build_p0_9` before falling back to `find`;
- update the run loop to execute the resolved path and report it in `P0_9_VALIDATION_RESULT.txt`;
- keep all P0.9 production Fortran source unchanged.

No physics, RHS coupling, pressure/projection/RK3/channel-forcing, force-buffer, IBM spreading, or production DNS logic was changed by this repair. P0.9 remains blocked until the repaired validation script is rerun and reports `Result: PASS`.
