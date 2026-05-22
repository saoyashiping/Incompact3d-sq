# Stage 9.8 restart I/O regression

Stage 9.8 validates real production restart write/read/continue in no-fibre/no-coupling mode.

- Unlike Stage 9.7 (stats/visu/coarse I/O), Stage 9.8 is restart-focused.
- Two phases per np: `cold` (write restart) and `restart` (read + continue).
- Restart file pattern checks: `restart*` files must exist and be non-empty.
- Restart signature uses global sums and max-abs of `ux/uy/uz` with tolerance `X3D_STAGE9_8_RESTART_SIGNATURE_TOL`.
- Finite checks include velocity/pressure/divergence fields.
- Timeout per phase: `STAGE9_8_TIMEOUT_SEC`.
- Pass requires cold and restart phase PASS lines and dat status pass.

Manual run:
`bash stage9_checks/run_stage9_8_restart_io_regression.sh`

Expected lines:
- `STAGE 9.8 RESTART IO REGRESSION VERDICT: PASS`
- `STAGE 9.8 FINAL VERDICT: PASS`
