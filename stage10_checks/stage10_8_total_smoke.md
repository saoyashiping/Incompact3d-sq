# Stage 10.8 Total Smoke / Closure

## Stage 10.8 target

- Stage 10 total smoke
- Stage 10 closure
- production no-op hook skeleton closure

## Mathematical / physical meaning

- `f_fsi = 0`
- `RHS_stage10 = RHS_stage9`
- hook active but no-op

## Closure evidence

- Stage 10.4 np=1 no-op invariance
- Stage 10.5 np=1/2/4 parallel no-op invariance
- Stage 10.6 RHS contamination audit
- Stage 10.7 restart/stats/visu/coarse I/O no-op validation

## Explicit old-logic warning

- do not run Stage 10.2 / Stage 10.3 by default
- do not forbid valid guarded hook calls
- do not broad-grep `ibm` / `structure` / `fibre_`
- do not classify `no_*` diagnostic fields as contamination
- do not bare-run `xcompact3d` for Stage 10.3 semantics

## What is intentionally not done

- no real IBM
- no real fibre force
- no RHS injection
- no structure advance
- no two-way coupling

## Pass criteria

Stage 10.8 passes only if:

1. required build targets pass,
2. Stage 10.4/10.5/10.6/10.7 pass,
3. Stage 10 hook no-op safety dat keys pass,
4. `stage10_outputs/stage10_8_total_smoke.dat` final status is `1`.

`STAGE10_CLOSED.md` is generated only on full pass.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE10_SKIP_PREREQS=1 \
bash stage10_checks/run_stage10_8_total_smoke.sh
```
