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
