# Stage 15.0 Configuration / Global Switches Only

## Objective

Stage 15.0 introduces only the guarded Stage 15 configuration layer for future production flexible-fibre structure advance. It does not allocate production structure state, does not advance the fibre, does not update fibre position or velocity, does not compute bending or tension, does not enable wall/contact or multi-fibre logic, and does not change fluid/RHS behavior.

## Mathematical context

Future Stage 15 work will eventually connect the structure equation

```text
rho_tilde * X_tt = d_s(T d_s X) - gamma * d_ssss X - F_fs
```

with inextensibility

```text
d_s X · d_s X = 1
```

Stage 15.0 does not solve either equation. It only exposes guarded configuration controls.

## Configuration controls

The Stage 15.0 wrapper supports:

- `DECOMP2D_ROOT`
- `BUILD_DIR`
- `MPIEXEC`
- `MPIEXEC_FLAGS`
- `STAGE15_0_RUN_STAGE14_11`
- `STAGE15_0_REQUIRE_STAGE14_CLOSED`
- `STAGE15_0_ENABLE`
- `STAGE15_0_STRUCTURE_ADVANCE_ENABLE`
- `STAGE15_0_DIAGNOSTIC_ONLY`

Safe defaults are:

```text
BUILD_DIR=build_stage9
MPIEXEC=mpirun
MPIEXEC_FLAGS=
STAGE15_0_RUN_STAGE14_11=0
STAGE15_0_REQUIRE_STAGE14_CLOSED=1
STAGE15_0_ENABLE=0
STAGE15_0_STRUCTURE_ADVANCE_ENABLE=0
STAGE15_0_DIAGNOSTIC_ONLY=1
```

The Fortran configuration module maps these wrapper controls to Stage 15 configuration values and forces Stage 15.0 to remain diagnostic-only with structure advance disabled.

## Static audit

The Stage 15.0 wrapper verifies that:

- Stage 15 source files do not call structure advance.
- Stage 15 source files do not call bending solve.
- Stage 15 source files do not call tension solve.
- Stage 15 source files do not update fibre position.
- Stage 15 source files do not update fibre velocity.
- Stage 15 source files do not call wall/contact or multi-fibre logic.
- Stage 15 source files do not modify pressure/projection/Poisson/RK3/channel forcing.
- Stage 15 source files do not call production IBM forcing outside the approved Stage 11–14 path.
- The forbidden Stage 14 `stage14_get_injection_gain() == 0.0` hook-registration gate is absent.
- Stage 11.5, Stage 13, and Stage 14 production diagnostics still exist in source.
- Rank0-safe diagnostic writing protections are preserved.
- Stage 13 force-density diagnostic sampling repair is preserved.

## Build/runtime check

The wrapper builds and runs only the Stage 15.0 config check target:

```text
fibre_stage15_config_check
```

The check writes:

```text
stage15_outputs/fibre_stage15_0_config.dat
```

and verifies safe defaults, diagnostic-only behavior, structure-advance-disabled behavior, closure-gating controls, and no production behavior change.

## Pass criteria

Stage 15.0 passes only if:

- Stage 14 is closed when `STAGE15_0_REQUIRE_STAGE14_CLOSED=1`.
- The static regression audit passes.
- `fibre_stage15_config_check` builds and runs.
- Stage 15 defaults are safe.
- Stage 15 structure advance remains disabled.
- Stage 15 diagnostic-only mode remains active.
- No structure state allocation, structure advance, bending/tension solve, fibre position update, fibre velocity update, wall/contact, or multi-fibre logic is detected.
- No pressure/projection/Poisson/RK3/channel-forcing changes are detected.
- No Stage 11–14 diagnostic regression is detected.

The wrapper prints:

```text
STAGE 15.0 CONFIG VERDICT: PASS
STAGE 15.0 FINAL VERDICT: PASS
```

or FAIL with explicit reasons.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE15_0_RUN_STAGE14_11=0 \
STAGE15_0_REQUIRE_STAGE14_CLOSED=1 \
STAGE15_0_ENABLE=0 \
STAGE15_0_STRUCTURE_ADVANCE_ENABLE=0 \
STAGE15_0_DIAGNOSTIC_ONLY=1 \
bash stage15_checks/run_stage15_0_config.sh
```

If `stage14_checks/STAGE14_CLOSED.md` is not present, rerun Stage 14.11 manually or set `STAGE15_0_RUN_STAGE14_11=1` to request the prerequisite from this wrapper.


## Stage 15.0 repair notes

The static call scanner checks only the called routine name. It must not flag ordinary configuration calls such as `get_environment_variable` merely because an environment-variable string contains words like `STRUCTURE_ADVANCE`.

Stage 15.0 requires the versioned Stage 14 closure record `stage14_checks/STAGE14_CLOSED.md` unless `STAGE15_0_RUN_STAGE14_11=1` is explicitly used to regenerate closure evidence.
