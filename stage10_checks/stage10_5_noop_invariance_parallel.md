# Stage 10.5 Parallel No-Op Invariance

## Stage 10.5 target

Stage 10.5 validates `np=1/2/4` baseline vs Stage 10 hook-enabled no-op invariance.

## Mathematical/physical meaning

- `f_fsi = 0`
- Stage 10 hook is active in hook runs, but no-op
- no rank-dependent side effect
- final DNS signatures should remain unchanged for each matching `np`

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

A metric passes when:

`abs(delta) <= max(abs_tol, rel_tol * max(1, abs(reference)))`

Defaults:

- `STAGE10_5_INVARIANCE_ABS_TOL=1.0e-12`
- `STAGE10_5_INVARIANCE_REL_TOL=1.0e-14`

## Explicit warning

- Stage 10.5 does **not** rerun Stage 10.2/10.3 by default.
- This avoids re-entering old false-positive audit logic.
- Stage 10.5 validates runtime evidence directly (dat outputs and comparisons).

## Intentionally not done

- no real IBM interpolation/spreading
- no real fibre force
- no RHS injection
- no structure advance
- no two-way coupling

## Pass criteria

- required build targets pass
- baseline and hook-enabled Stage 9.9 runs complete
- Stage 10.3 hook no-op evidence dat passes safety keys
- all `np=1/2/4` baseline-vs-hook signature metrics pass mixed tolerance
- `stage10_outputs/stage10_5_noop_invariance_parallel.dat` final status is `1`

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE10_5_INVARIANCE_ABS_TOL=1.0e-12 \
STAGE10_5_INVARIANCE_REL_TOL=1.0e-14 \
STAGE10_5_RUN_PREREQS=0 \
bash stage10_checks/run_stage10_5_noop_invariance_parallel.sh
```
