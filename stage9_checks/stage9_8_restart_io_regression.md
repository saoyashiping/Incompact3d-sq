# Stage 9.8 restart I/O regression

Stage 9.8 validates real production restart write/read/continue in no-fibre/no-coupling mode.

- Unlike Stage 9.7 (stats/visu/coarse I/O), Stage 9.8 is restart-focused.
- Compile portability fix: avoid character-valued `MERGE` with unequal string lengths; use explicit `if/else` branch prints.
- Two phases per np: `cold` (write restart) and `restart` (read + continue).
- Stage 9.8 gate generates per-np temporary inputs under `stage9_outputs/`.
- Cold phase uses generated cold input with `irestart=0` and `icheckpoint=1`.
- Restart phase uses generated restart input with `irestart=1`, `icheckpoint=1`, and `ifirst=STAGE9_8_MAX_STEPS_BEFORE_RESTART+1` (default `ifirst=4` when max steps before restart is 3).
- Cold phase writes checkpoint/restart data after `STAGE9_8_MAX_STEPS_BEFORE_RESTART` complete outer steps; the restart input `ifirst` aligns the real restart-read branch with that written checkpoint state.
- Stale restart/checkpoint files are removed before each np cold phase.
- Restart/checkpoint detection uses a real matched file path from `find ... -print -quit`, not bare `find` exit status.
- Restart-read status is recorded only in real production `irestart /= 0` branch around `restart(...,0)` in `xcompact3d`.
- Restart file pattern checks require non-empty `checkpoint*`/`restart*` files before restart phase starts.
- Restart signature persists across runs through `X3D_STAGE9_8_SIGNATURE_FILE` and compares global sums/max-abs of `ux/uy/uz` with tolerance `X3D_STAGE9_8_RESTART_SIGNATURE_TOL`.
- Restart signature comparison is performed at the restart **read point** (immediately after real production `restart(...,0)`), before additional continuation steps advance the state.
- Post-restart continuation steps do not overwrite the read-point signature decision; read-point comparison is one-shot and locked.
- On signature failure, dat/log include absolute differences `d_sum_*`, `d_max_*`, and `stage9_8_signature_tolerance`.
- Finite checks include velocity/pressure/divergence fields.
- Timeout per phase: `STAGE9_8_TIMEOUT_SEC`.
- Stage 9.8 remains no-fibre/no-coupling.

Manual run:
`bash stage9_checks/run_stage9_8_restart_io_regression.sh`

Expected lines:
- `STAGE 9.8 RESTART IO REGRESSION VERDICT: PASS`
- `STAGE 9.8 FINAL VERDICT: PASS`
