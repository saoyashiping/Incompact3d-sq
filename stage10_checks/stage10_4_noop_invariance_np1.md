# Stage 10.4 No-Op Invariance (np=1)

## Stage 10.4 target

Stage 10.4 verifies np=1 invariance between:

1. baseline deterministic no-fibre run, and
2. Stage 10 production-hook-enabled no-op run.

## Mathematical/physical meaning

- `f_fsi = 0`
- hook is active but no-op
- final DNS signatures should remain unchanged within strict mixed tolerance

## Compared signatures

- `stage9_9_signature_sum_ux`
- `stage9_9_signature_sum_uy`
- `stage9_9_signature_sum_uz`
- `stage9_9_signature_max_ux`
- `stage9_9_signature_max_uy`
- `stage9_9_signature_max_uz`
- `stage9_9_signature_l2_ux`
- `stage9_9_signature_l2_uy`
- `stage9_9_signature_l2_uz`

## Tolerance formula

A metric passes if:

`abs(delta) <= max(abs_tol, rel_tol * max(1, abs(reference)))`

Defaults:

- `STAGE10_4_INVARIANCE_ABS_TOL=1.0e-12`
- `STAGE10_4_INVARIANCE_REL_TOL=1.0e-14`

## Intentionally not done

- no np=2/4 comparison (reserved for later stage)
- no IBM interpolation/spreading
- no real fibre force
- no RHS injection
- no structure advance

## Pass criteria

- required build targets pass
- Stage 10.0–10.3 prerequisites pass (unless skipped)
- baseline and hook runs both produce valid `np=1` Stage 9.9 dat
- Stage 10.3 hook dat confirms active no-op safety
- all signature comparisons pass mixed tolerance
- `stage10_outputs/stage10_4_noop_invariance_np1.dat` reports final status `1`

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE10_4_TIMEOUT_SEC=240 \
STAGE10_4_INVARIANCE_ABS_TOL=1.0e-12 \
STAGE10_4_INVARIANCE_REL_TOL=1.0e-14 \
bash stage10_checks/run_stage10_4_noop_invariance_np1.sh
```
