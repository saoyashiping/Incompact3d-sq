# Stage 9.7 stats/visu/coarse-io smoke

Stage 9.7 verifies real production channel DNS output plumbing for statistics, visualization, and coarse/decomp I/O while staying no-fibre/no-coupling.

- Distinct from Stage 9.5 (projection/divergence) and Stage 9.6 (RK3/RHS/mass-flux).
- Uses real production calls (`postprocessing` -> `run_postprocessing` -> `write_snapshot`/`overall_statistic`, plus normal decomp IO init/finalise).
- No Stage 8 fibre/IBM coupling is enabled.

## Checked paths
- Statistics path reached and finite-status tracking.
- Visualization/snapshot path reached and finite-status tracking.
- Coarse/decomp I/O descriptor/open/write/close status tracking.
- Output file patterns checked by script: non-empty files under `data/` and `statistics/`.

## Finite-field criteria
- Output-field finite status checks real `ux1,uy1,uz1,pp3(:,:,:,1),divu3` before end-of-step audit.

## Not tested yet
- Detailed scientific correctness of statistics values.
- Fibre coupling, IBM force spreading, structure advancement.

## Pass criteria
- Stage 9.7 diagnostic status line reports `stage9_7_stats_visu_io_smoke_status 1`.
- Program prints `STAGE 9.7 STATS VISU IO SMOKE VERDICT: PASS`.
- Gate script prints `STAGE 9.7 FINAL VERDICT: PASS`.

## Manual run
`bash stage9_checks/run_stage9_7_stats_visu_io_smoke.sh`


## Gate orchestration
- Stage 9.7 calls earlier Stage 9 gates as prerequisites by default.
- Set `STAGE9_SKIP_PREREQS=1` to skip nested prerequisite calls when Stage 9.7 is invoked from a higher-level orchestrator.
- This flag only skips prerequisite gate invocations; each stage still executes its own build/run/log/dat validation.
- Stage 9.7 remains no-fibre/no-coupling/no-IBM-injection.

## Bounded smoke execution

- Stage 9.7 smoke runs are bounded by `X3D_STAGE9_7_MAX_STEPS` and exit after complete outer steps.
- The gate script enforces `STAGE9_7_TIMEOUT_SEC` per np run to prevent hangs.
- The executable prints per-step Stage 9.7 progress and a final-audit start message.
- If timeout occurs, inspect the printed log tail.
- Stage 9.7 remains no-fibre/no-coupling/no-IBM-injection.


## Environment variable parsing guard
- Stage 9.7 env parsing in Fortran must call `get_environment_variable` with keywords (`value=...`, `status=...`).
- Do not use positional argument 3 for status: positional arg 3 is **length**, not status, and can silently disable smoke-mode detection.
