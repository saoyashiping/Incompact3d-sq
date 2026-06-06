# Stage 16.11 short-time stability / bounded-energy smoke for first one-fibre FSI

Stage 16.11 is a diagnostic-only short-time boundedness smoke after the Stage 16.10 contamination audit. It reuses the closed Stage 16.7 small-lambda one-fibre diagnostic path over a short controlled sequence and summarizes bounded structure/force/RHS/fluid/work/energy proxy evidence. The approved path remains:

```text
Stage 11 fluid-to-fibre sampling
-> Stage 12 feedback-force candidate
-> Stage 16.4 structure-side fluid-on-fibre force input
-> controlled structure-state diagnostics
-> Stage 13 force-density candidate
-> Stage 14 RHS diagnostic / controlled injection path
```

Stage 16.11 does not add production hooks, does not connect to `xcompact3d.f90`, does not duplicate closed-loop physics, and does not modify pressure/projection/Poisson/RK3/channel-forcing numerics. It does not activate wall/contact, multi-fibre, legacy IBM, bending, tension, or full-structure production solves.

## User command

```bash
bash stage16_checks/run_stage16_11_short_time_stability_smoke.sh
```

The wrapper defaults to `BUILD_DIR=build_stage9`, configures the build directory with CMake if needed, builds only the already closed Stage 16.7 target `fibre_stage16_small_lambda_response_check`, runs that diagnostic for `STAGE16_11_NSTEPS=5`, writes Stage 16.11 runtime evidence, and invokes the Stage 16.11 assertion helper.

## Default controls

Key defaults are:

- `STAGE16_11_RUN_STAGE16_10=0`
- `STAGE16_11_REQUIRE_STAGE14_CLOSED=1`
- `STAGE16_11_REQUIRE_STAGE15_CLOSED=1`
- `STAGE16_11_REQUIRE_STAGE16_10=1`
- `STAGE16_11_ACCEPT_STAGE16_10_CLOSED_EVIDENCE=1`
- `STAGE16_11_NP=2`
- `STAGE16_11_NPTS=8`
- `STAGE16_11_NSTEPS=5`
- `STAGE16_11_SMALL_LAMBDA=1.0e-8`
- `STAGE16_11_DT=1.0e-4`
- `STAGE16_11_MAX_FORCE_INPUT=1.0e-6`
- `STAGE16_11_MAX_STRUCTURE_UPDATE=1.0e-12`
- `STAGE16_11_MAX_VELOCITY_UPDATE=1.0e-10`
- `STAGE16_11_MAX_ACCELERATION_UPDATE=1.0e-6`
- `STAGE16_11_MAX_RHS_INCREMENT=1.0e-8`
- `STAGE16_11_MAX_FLUID_DELTA=1.0e-8`
- `STAGE16_11_MAX_WORK_PROXY=1.0e-14`
- `STAGE16_11_MAX_ENERGY_PROXY=1.0e-14`
- `STAGE16_11_MAX_GROWTH_RATIO=10.0`

## Runtime diagnostics

The wrapper writes short-time runtime evidence to:

```text
stage16_outputs/stage16_11_short_time_stability_evidence.dat
```

The required Stage 16.11 summary is:

```text
stage16_outputs/fibre_stage16_11_short_time_stability_smoke.dat
```

The summary contains the required closed-loop status fields, bounded position/velocity/acceleration/force/Stage13/RHS/fluid statuses, finite and bounded work/energy proxy statuses, growth ratio and no-runaway statuses, no-contamination statuses, rank0/rank-corruption statuses, Stage 14/15/16.1-16.10 regression statuses, and `final_status`.

## PASS evidence

A passing run prints:

```text
STAGE 16.11 SHORT-TIME STABILITY SMOKE VERDICT: PASS
STAGE 16.11 FINAL VERDICT: PASS
```

and the summary file reports:

- `closed_loop_path_status 1`
- `position_update_bounded_status 1`
- `velocity_update_bounded_status 1`
- `acceleration_update_bounded_status 1`
- `force_input_bounded_status 1`
- `stage13_force_density_bounded_status 1`
- `rhs_increment_bounded_status 1`
- `fluid_signature_bounded_status 1`
- `work_proxy_finite_status 1` and `work_proxy_bounded_status 1`
- `energy_proxy_finite_status 1` and `energy_proxy_bounded_status 1`
- `growth_ratio_bounded_status 1`
- `no_runaway_growth_status 1`
- all approved-chain and no-contamination statuses equal to `1`
- all Stage 14/15/16.1-16.10 regression statuses equal to `1`
- `final_status 1`

## False-positive-safe static audit policy

The Stage 16.11 helper intentionally reuses the corrected Stage 16.10 / Stage 16.9 / Stage 16.8 / Stage 16.7 / Stage 16.6 / Stage 16.5 / Stage 16.4 helper pattern:

- documentation is checked for required-file existence, not scanned as executable regression evidence;
- negative-check strings are not treated as production behavior;
- regex literals such as `rg[[:space:]]` are not treated as real `rg` command usage;
- legitimate Stage 13.5 conservation/sign audit files are not classified as old production diagnostic regressions;
- old Stage 13.5 production force-density names are rejected only in real production/check logic;
- if evidence cannot be distinguished safely, the helper fails closed with explicit reasons rather than silently passing.

## Assumptions and risks

Stage 16.11 intentionally reuses closed Stage 16.7/16.8/16.9/16.10 evidence instead of adding a new Fortran physics module. The work and energy values are simple diagnostic proxies derived from the controlled Stage 16.7 force/update signatures, not production physical energy accounting.

## Stage 16.7 reused-runtime evidence note

The Stage 16.11 wrapper reuses the passed Stage 16.7 small-lambda runtime output as the per-step kernel. That Stage 16.7 file is produced by Stage 16.7 itself and intentionally does **not** contain a self-referential `stage16_7_regression_status` field. Stage 16.11 must infer Stage 16.7 closed evidence from `final_status = 1`, preserved `stage16_6_regression_status = 1`, and the existence of the Stage 16.7 wrapper/helper/doc/source/check files. Missing `stage16_7_regression_status` in the reused Stage 16.7 runtime file is not a regression and must not be treated as a failed status.

