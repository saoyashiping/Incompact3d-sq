# Stage 20.0 preflight boundary

Stage 20.0 is a diagnostic-only boundary before Stage 20 begins. It accepts
Stage 19 closure using available Stage 19.11 closure evidence, including
source-only archive evidence, and it does **not** require rerunning any previous
stage.

## Safety boundary

Stage 20.0 declares these conceptual gates without activating production
coupling:

- `STAGE20_ENABLE`: true only for preflight diagnostic
- `STAGE20_TWOWAY_COUPLING_ENABLE`: false
- `STAGE20_FLUID_TO_STRUCTURE_ENABLE`: false
- `STAGE20_STRUCTURE_TO_FLUID_ENABLE`: false
- `STAGE20_RHS_COUPLING_ENABLE`: false
- `STAGE20_LAMBDA_COUPLING`: 0.0
- `STAGE20_DIAGNOSTIC_ONLY`: true
- `STAGE20_FAIL_CLOSED`: true

Stage 20.0 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 file. It does not run MPI, DNS, previous stages, builds, or
production validation.

## Source-only compatibility

Missing old stage outputs and missing old closure files are allowed when Stage
19 closure is accepted from Stage 19.11 PASS evidence, Stage 19 closure
metadata, source-level Stage 19.11 helper/wrapper/documentation, or manual/user
accepted source-only closure. A missing `.git` directory is not an unknown
failure.

## Next stage

Stage 20.1: two-way coupling config and lambda gate.

## Manual command

```bash
stage20_checks/run_stage20_0_preflight_boundary.sh
```

Expected output includes:

```text
STAGE 20.0 PREFLIGHT BOUNDARY VERDICT: PASS
STAGE 20.0 FINAL VERDICT: PASS
```
