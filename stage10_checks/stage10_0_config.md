# Stage 10.0 config layer

## Target
Stage 10.0 introduces only production-coupling configuration switches and a standalone check target. It does not connect any production coupling hook.

## Mathematical/physical meaning
The no-fibre DNS remains unchanged. `f_fsi` exists only formally, and in Stage 10.0:

- `f_fsi = 0`

No fluid equation term is added or modified.

## Environment variables
- `X3D_STAGE10_PRODUCTION_HOOK`
- `X3D_STAGE10_FORCE_NOOP`
- `X3D_STAGE10_MAX_STEPS`

## Intentionally not done in Stage 10.0
- no `xcompact3d` main-loop hook call
- no production coupling insertion
- no IBM interpolation/spreading
- no fibre state allocation
- no force injection

## Pass criteria
- `xcompact3d` builds
- `fibre_stage10_config_check` builds
- check executable prints `STAGE 10.0 CONFIG VERDICT: PASS`
- dat statuses report requested/noop/no-fibre/no-force/no-rhs-injection/config all passing

## Manual command
```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
bash stage10_checks/run_stage10_0_config.sh
```
