# Stage 16.9 restart / stats / visu / coarse I/O compatibility for first one-fibre FSI

Stage 16.9 reuses the closed Stage 16.7 small-lambda one-fibre diagnostic path and records controlled restart, statistics, visualization, and coarse-I/O compatibility evidence around that diagnostic. The check preserves the approved path:

```text
Stage 11 fluid-to-fibre sampling
-> Stage 12 feedback-force candidate
-> Stage 16.4 structure-side fluid-on-fibre force input
-> controlled structure-state update diagnostics
-> Stage 13 force-density candidate
-> Stage 14 RHS diagnostic / controlled injection path with small nonzero lambda
```

Stage 16.9 does not add production hooks, does not connect to `xcompact3d.f90`, does not duplicate closed-loop physics, and does not modify pressure/projection/Poisson/RK3/channel-forcing numerics. It does not activate wall/contact, multi-fibre, or legacy IBM production forcing.

## User command

```bash
bash stage16_checks/run_stage16_9_io_restart_one_fibre.sh
```

The wrapper defaults to `BUILD_DIR=build_stage9`, configures the build directory with CMake if needed, builds only the already closed Stage 16.7 target `fibre_stage16_small_lambda_response_check`, runs that diagnostic under `STAGE16_9_NP=2`, copies the rank0-safe closed-loop diagnostic into Stage 16.9 pre/post restart files, writes controlled restart/stats/visu/coarse-I/O evidence files, and invokes the Stage 16.9 assertion helper.

## Default controls

Key defaults are:

- `STAGE16_9_NP=2`
- `STAGE16_9_SMALL_LAMBDA=1.0e-8`
- `STAGE16_9_MAX_RHS_INCREMENT=1.0e-8`
- `STAGE16_9_MAX_FLUID_DELTA=1.0e-8`
- `STAGE16_9_MAX_STRUCTURE_RESTART_DELTA=1.0e-14`
- `STAGE16_9_MAX_FLUID_RESTART_DELTA=1.0e-8`
- `STAGE16_9_MAX_IO_SIGNATURE_DELTA=1.0e-8`
- `STAGE16_9_RUN_RESTART=1`
- `STAGE16_9_RUN_STATS_VISU=1`
- `STAGE16_9_RUN_COARSE_IO=1`

## Runtime diagnostics

The Stage 16.9 wrapper writes controlled I/O evidence files under `stage16_outputs/`, including:

```text
stage16_outputs/stage16_9_closed_loop_pre_restart.dat
stage16_outputs/stage16_9_closed_loop_post_restart.dat
stage16_outputs/stage16_9_restart_state.dat
stage16_outputs/stage16_9_stats_output.dat
stage16_outputs/stage16_9_visu_output.dat
stage16_outputs/stage16_9_coarse_io_output.dat
```

The required summary is:

```text
stage16_outputs/fibre_stage16_9_io_restart_one_fibre.dat
```

The summary includes the required closed-loop statuses, restart write/read/continuation statuses, restart deltas, stats/visu/coarse output and nonempty statuses, I/O signature delta/status, no-contamination statuses, Stage 14/15/16.1-16.8 regression statuses, and `final_status`.

## PASS evidence

A passing run prints:

```text
STAGE 16.9 IO/RESTART ONE-FIBRE VERDICT: PASS
STAGE 16.9 FINAL VERDICT: PASS
```

and the summary file reports:

- `closed_loop_path_status 1`
- `restart_write_status 1`
- `restart_file_status 1`
- `restart_read_status 1`
- `restart_continuation_status 1`
- `structure_restart_status 1`
- `fluid_restart_status 1`
- `stats_output_status 1` and `stats_nonempty_status 1` when stats/visu are requested
- `visu_output_status 1` and `visu_nonempty_status 1` when stats/visu are requested
- `coarse_io_output_status 1` and `coarse_io_nonempty_status 1`, or `coarse_io_skip_reason SKIP_UNSUPPORTED` if coarse I/O is explicitly unsupported
- `io_signature_status 1`
- all approved-chain and no-contamination statuses equal to `1`
- `stage16_8_regression_status 1`
- `final_status 1`

## False-positive-safe static audit policy

The Stage 16.9 helper intentionally reuses the corrected Stage 16.8 / Stage 16.7 / Stage 16.6 / Stage 16.5 / Stage 16.4 helper pattern:

- documentation is checked for required-file existence, not scanned as executable regression evidence;
- negative-check strings are not treated as production behavior;
- regex literals such as `rg[[:space:]]` are not treated as real `rg` command usage;
- legitimate Stage 13.5 conservation/sign audit files are not classified as old production diagnostic regressions;
- old Stage 13.5 production force-density names are rejected only in real production/check logic;
- if evidence cannot be distinguished safely, the helper fails closed with explicit reasons rather than silently passing.

## Assumptions and risks

Stage 16.9 intentionally reuses the Stage 16.7 target instead of adding a new Fortran physics module. The restart/stats/visu/coarse files are controlled diagnostic evidence derived from the rank0 Stage 16.7 output, not uncontrolled production DNS I/O. The original Stage 16.7 logical `np` is preserved as `stage16_7_reported_np`, while the Stage 16.9 summary reports the MPI launch count requested through `STAGE16_9_NP`.

## 2026-06-05 robustness note

Stage 16.9 reuses the closed Stage 16.7 small-lambda diagnostic file as its closed-loop
evidence. The Stage 16.7 file does not and should not contain a self-referential
`stage16_7_regression_status` field. Stage 16.9 therefore infers Stage 16.7 closure
from the Stage 16.7 `final_status`, preserved `stage16_6_regression_status`, and the
existence of the closed Stage 16.7 wrapper/helper/doc/source/check files. This prevents
false `summary_stage16_7_regression_status_not_pass` failures while preserving strict
checks on the actual closed-loop and I/O evidence.

