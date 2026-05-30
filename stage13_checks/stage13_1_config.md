# Stage 13.1 Eulerian Force-Density Candidate Configuration

## Target

Stage 13.1 adds the configuration layer for the future Stage 13 production Eulerian force-density candidate. It defines global switches, readonly controls, point-count limits, and spreading-normalization selection while remaining diagnostic-only.

Stage 13.1 does not connect to `xcompact3d`, does not add a production hook, and does not introduce force-density physics.

## Mathematical and physical meaning

Stage 13 will later introduce the Eulerian force-density candidate

```text
f_i_cand = (1 / DeltaV_i) * sum_k F_fibre_to_fluid_cand(k) * W_ik * Delta_s_k
```

Stage 13.1 does not compute `f_i_cand`, does not allocate an Eulerian force-density buffer, and does not spread Lagrangian force. The production equations remain unchanged:

```text
f_fsi = 0
RHS_stage13.1 = RHS_stage13.0 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

No spreading is performed and RHS injection remains forbidden.

## Environment variables

The Stage 13.1 configuration module reads:

```text
X3D_STAGE13_FORCE_DENSITY_CANDIDATE
X3D_STAGE13_FORCE_READONLY
X3D_STAGE13_SPREADING_READONLY
X3D_STAGE13_MAX_POINTS
X3D_STAGE13_MAX_EULERIAN_POINTS
X3D_STAGE13_SPREADING_NORMALIZATION
```

Defaults are safe and diagnostic-only: Stage 13 disabled, force readonly enabled, spreading readonly enabled, `X3D_STAGE13_MAX_POINTS=8`, `X3D_STAGE13_MAX_EULERIAN_POINTS=64`, and conservative normalization.

## What is intentionally not done

Stage 13.1 intentionally does not add:

- an `xcompact3d` hook call;
- production main-loop insertion;
- force-density buffer allocation;
- Lagrangian-to-Eulerian spreading;
- Eulerian force-density computation;
- RHS injection;
- IBM spreading;
- feedback force application to the fluid;
- two-way force;
- fibre structure advance.

## Pass criteria

The Stage 13.1 gate passes only if:

1. `xcompact3d`, required Stage 11 checks, required Stage 12 checks, and `fibre_stage13_config_check` build.
2. The standalone Stage 13.1 config check prints `STAGE 13.1 CONFIG VERDICT: PASS`.
3. `stage13_outputs/fibre_stage13_1_config.dat` reports requested, readonly, spreading-readonly, valid max-point parsing, valid max-Eulerian-point parsing, valid normalization parsing, no force-density allocation, no spreading, no RHS injection, no IBM spreading, no feedback application, no two-way force, no structure advance, no fluid-field access, no fluid-field modification, and config status all as pass.
4. No closed Stage 10, Stage 11, Stage 12, or Stage 13.0 scripts are run by default.

The gate writes `stage13_outputs/stage13_1_config_gate.dat` and prints:

```text
STAGE 13.1 FINAL VERDICT: PASS
```

Failures print explicit non-empty reasons.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE13_1_RUN_STAGE13_0=0 \
bash stage13_checks/run_stage13_1_config.sh
```
