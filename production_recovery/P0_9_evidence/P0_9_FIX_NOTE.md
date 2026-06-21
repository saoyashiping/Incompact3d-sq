# P0.9 validation script repair note

The first P0.9 failure on the user machine was caused by the validation script depending on `rg` (ripgrep). The user environment did not have `rg`, so static `contains()` checks returned FAIL even though the corresponding source code, xcompact3d imports, state commit helper, and CMake target were present.

Fix:
- add a portable grep fallback for all static regex checks;
- replace evidence-producing `rg -n` calls with a wrapper that uses `rg` when available and `grep -En` otherwise;
- replace required-token checks with `grep -Fq`;
- keep production source behavior unchanged;
- keep P0.9 production status as STILL BLOCKED until validation is rerun and passes.

No physics, RHS coupling, pressure/projection/RK3/channel-forcing logic was changed by this repair.
