# Stage 9.8 restart I/O regression

Stage 9.8 validates real production restart write/read/continue in no-fibre/no-coupling mode.

- Unlike Stage 9.7 (stats/visu/coarse I/O), Stage 9.8 is restart-focused.
- Compile portability fix: avoid character-valued `MERGE` with unequal string lengths; use explicit `if/else` branch prints.
- Two phases per np: `cold` (write restart) and `restart` (read + continue).
- Cold phase uses original input file.
- Restart phase uses generated input `stage9_outputs/stage9_8_input_restart_np<np>.i3d` with `irestart = 1`.
- Restart-read status is recorded only in real production `irestart /= 0` branch around `restart(...,0)` in `xcompact3d`.
- Restart file pattern checks: `restart*` files must exist and be non-empty before restart phase starts.
- Restart signature persists across runs through `X3D_STAGE9_8_SIGNATURE_FILE` and compares global sums/max-abs of `ux/uy/uz` with tolerance `X3D_STAGE9_8_RESTART_SIGNATURE_TOL`.
- Finite checks include velocity/pressure/divergence fields.
- Timeout per phase: `STAGE9_8_TIMEOUT_SEC`.
- Stage 9.8 remains no-fibre/no-coupling.

Manual run:
`bash stage9_checks/run_stage9_8_restart_io_regression.sh`

Expected lines:
- `STAGE 9.8 RESTART IO REGRESSION VERDICT: PASS`
- `STAGE 9.8 FINAL VERDICT: PASS`


- Stage 9.8 gate generates per-np temporary inputs under `stage9_outputs/`.
- Cold input sets `irestart=0` and `icheckpoint=1`; restart input sets `irestart=1` and `icheckpoint=1`.
- Stale restart/checkpoint files are removed before each np cold phase.
- Restart/checkpoint detection uses a real matched file path from `find ... -print -quit`, not bare `find` exit status.
