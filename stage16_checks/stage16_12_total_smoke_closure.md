# Stage 16.12 total smoke and closure

Stage 16.12 is the final Stage 16 closure aggregation for the first controlled
one-fibre FSI diagnostic path. It does not add physics, Fortran sources, CMake
targets, or production hooks. The check reads existing Stage 16.0--16.11 runtime
or closed-file evidence, confirms the approved Stage 11 -> Stage 12 -> Stage
16.4 -> Stage 13 -> Stage 14 chain, and writes the Stage 16.12 summary.

## Manual command

```bash
bash stage16_checks/run_stage16_12_total_smoke_closure.sh
```

By default the wrapper runs no physics and builds nothing. It creates
`stage16_outputs/`, invokes the Stage 16.12 Python helper, writes
`stage16_outputs/fibre_stage16_12_total_smoke_closure.dat`, and generates
`stage16_checks/STAGE16_CLOSED.md` only after a full Stage 16.12 PASS.

## Defaults and environment controls

- `BUILD_DIR=build_stage9`
- `MPIEXEC=mpirun`
- `MPIEXEC_FLAGS` may be empty
- `STAGE16_12_RUN_STAGE16_11=0`
- `STAGE16_12_REQUIRE_STAGE14_CLOSED=1`
- `STAGE16_12_REQUIRE_STAGE15_CLOSED=1`
- `STAGE16_12_REQUIRE_STAGE16_11=1`
- `STAGE16_12_ACCEPT_STAGE16_11_CLOSED_EVIDENCE=1`
- `STAGE16_12_ENABLE=1`
- `STAGE16_12_DIAGNOSTIC_ONLY=1`
- `STAGE16_12_GENERATE_CLOSURE_FILE=1`

`DECOMP2D_ROOT` is accepted for consistency with earlier wrappers but is not used
unless a future maintainer explicitly opts into a build path. Stage 16.12 should
remain evidence aggregation and closure, not new runtime physics.

## Evidence policy

The helper accepts a prior stage through its runtime `.dat` file when the file is
present and has `final_status 1`. If prior runtime files are absent in a fresh
checkout, it accepts closed evidence only when the prior stage wrapper, helper,
documentation, and required source/check files exist, and when Stage 14 and Stage
15 closure records are present.

The helper intentionally does not require old stages to contain self-referential
regression fields they never wrote. For example, Stage 16.7 can be accepted from
`final_status` plus prior-stage evidence or from closed-file evidence.

## False-positive protections

The Stage 16.12 helper reuses the corrected Stage 16.11 / 16.10 / 16.9 / 16.8 /
16.7 / 16.6 / 16.5 / 16.4 false-positive-safe audit pattern:

- Markdown is used for required-file presence only, not as executable behavior.
- Negative-check strings are not treated as real regressions.
- Regex literals such as `rg[[:space:]]` are not treated as actual `rg` command
  usage.
- Legitimate Stage 13.5 conservation/sign audit files are not classified as old
  production diagnostic regressions.
- Old Stage 13.5 production force-density names are a regression only in real
  production/check logic where they replace Stage 13.6 evidence.
- Any actual `rg` usage in the Stage 16.12 wrapper has a `grep` fallback.

## Required PASS prints

A successful run prints:

```text
STAGE 16.12 TOTAL SMOKE CLOSURE VERDICT: PASS
STAGE 16.12 FINAL VERDICT: PASS
```

and produces:

- `stage16_outputs/fibre_stage16_12_total_smoke_closure.dat`
- `stage16_checks/STAGE16_CLOSED.md` when `STAGE16_12_GENERATE_CLOSURE_FILE=1`
