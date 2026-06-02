# Stage 14.0 — RHS Injection Configuration Gate

## Target

Stage 14.0 introduces only the configuration layer for future controlled
RHS injection of the Stage 13 Eulerian force-density candidate.  It defines:

- global Stage 14 RHS-injection switches;
- the future injection gain `lambda_14`;
- the requirement that Stage 13 diagnostics remain available before later
  injection stages;
- diagnostic-only protection that remains forced on in Stage 14.0.

## Mathematical / physical meaning

Stage 13 provides the diagnostic Eulerian force-density candidate `f_i_cand`.
Later Stage 14 steps may use the formula:

```text
RHS_new = RHS_old + lambda_14 * f_i_cand
```

Stage 14.0 does **not** perform that formula.  At Stage 14.0:

```text
f_fsi = 0
RHS_stage14.0 = RHS_stage13 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

No RHS addition is performed, no production force is applied, and no fibre
structure equation is advanced.

## Environment variables

The Stage 14.0 configuration module parses:

```text
X3D_STAGE14_RHS_INJECTION
X3D_STAGE14_INJECTION_GAIN
X3D_STAGE14_MAX_STEPS
X3D_STAGE14_REQUIRE_STAGE13
X3D_STAGE14_DIAGNOSTIC_ONLY
```

Defaults keep Stage 14 disabled: RHS injection is off, `lambda_14 = 0`,
`max_steps = 3`, Stage 13 is required, and diagnostic-only mode is forced on.
Invalid gain, step-count, or dependency values fall back to safe defaults.

## Required principles

- RHS injection is default off.
- Stage 14 requires the Stage 13 Eulerian force-density candidate path.
- Later stages must validate `lambda_14 = 0` before any nonzero gain.
- Stage 14.0 does not advance the flexible-fibre structure.
- Closed Stage 10, Stage 11, Stage 12, and Stage 13 files are not edited.

## Intentionally not done

Stage 14.0 intentionally does not add:

- an `xcompact3d` hook call;
- a production main-loop insertion;
- RHS-injection buffer allocation;
- RHS addition;
- production IBM forcing application;
- feedback force application beyond the already-closed Stage 13 diagnostics;
- two-way force application;
- fibre structure advance.

## Pass criteria

The Stage 14.0 gate passes only when:

1. `xcompact3d` and all required closed-stage check targets still build.
2. `fibre_stage14_config_check` builds.
3. The default-disabled run reports Stage 14 unrequested and no RHS injection.
4. The requested zero-gain run reports the request and keeps `lambda_14 = 0`.
5. The invalid-input run falls back safely while keeping diagnostic-only mode.
6. All no-modification, no-hook, no-RHS-buffer, no-RHS-injection, no-IBM,
   no-feedback, no-two-way, and no-structure statuses remain `1`.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE14_0_RUN_STAGE13_CLOSURE=0 \
bash stage14_checks/run_stage14_0_config.sh
```
