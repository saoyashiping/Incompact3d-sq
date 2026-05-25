# Stage 10.3 Main No-Op Hook Connection

Stage 10.3 is the first stage where the Stage 10 hook API is connected in `xcompact3d` main-loop control flow, while keeping hook operators no-op.

## Physical meaning

- `f_fsi = 0`
- no-fibre DNS equations remain unchanged
- hook operators are present but return no physical action

## Hook placement

- `stage10_hook_init`: after init and before time loop
- `stage10_hook_pre_step`: at outer-step start
- `stage10_hook_pre_rhs`: before RK/RHS work in main loop
- `stage10_hook_post_projection`: after substep block completion at main-loop level
- `stage10_hook_post_step`: after completed outer step
- `stage10_hook_finalize`: after loop exit / before shutdown

## Forbidden placements

- Poisson internals
- projection-math internals
- derivative kernels
- RHS math-expression internals
- restart I/O internals
- Stage 9 diagnostics internals

## Safety guarantees

- no fluid or pressure arrays passed to Stage 10 hooks in this stage
- no force injection
- no IBM interpolation/spreading
- no fibre state allocation
- no fibre structure advance

## Pass criteria

- required build targets succeed
- Stage 10.0/10.1/10.2 prerequisites pass (unless skipped)
- `fibre_stage10_3_main_noop_hook.dat` reports required status keys
- `stage10_3_field_modified_status = 0`

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE10_3_TIMEOUT_SEC=240 \
bash stage10_checks/run_stage10_3_main_noop_hook.sh
```
