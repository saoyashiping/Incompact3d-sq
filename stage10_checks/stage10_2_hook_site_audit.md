# Stage 10.2 Hook Site Audit

## Stage 10.2 target

Stage 10.2 performs a static production hook placement audit and records safe future insertion points for Stage 10 hook operators, without connecting any hook call in production code.

## Mathematical/physical meaning

Production no-fibre DNS remains unchanged:

- `div(u) = 0`
- `du/dt + u·grad(u) = -grad(p) + nu laplacian(u) + f_drive + f_fsi`
- Stage 10.2 constraint: `f_fsi = 0`

No production force path is activated in Stage 10.2.

## Candidate hook sites (future Stage 10.3 connection points)

- `hook_init`: after field allocation/initialization and before entering outer time loop.
- `hook_pre_step`: at outer-step entry before RK substep work.
- `hook_pre_rhs`: before RHS assembly boundary, outside RHS math implementation.
- `hook_post_projection`: after projection/velocity correction synchronization point, outside Poisson internals.
- `hook_post_step`: after full outer-step completion at a safe diagnostics point.
- `hook_finalize`: after time-loop exit and before shutdown/final diagnostics.

## Unsafe / forbidden insertion sites

- inside Poisson internals
- inside projection mathematics
- inside derivative kernels
- inside restart I/O internals
- inside Stage 9 diagnostic logic
- inside closed Stage 9 gate internals

## Intentionally not done

- no `xcompact3d` hook call insertion
- no RHS injection
- no IBM interpolation/spreading
- no fibre state allocation
- no fibre structure advance

## Pass criteria

Stage 10.2 passes only when:

1. Build targets pass:
   - `xcompact3d`
   - `fibre_stage10_config_check`
   - `fibre_stage10_noop_hook_check`
2. Stage 10.0 and Stage 10.1 prerequisite checks pass.
3. Candidate future insertion points are identified.
4. No active `stage10_hook_*` production calls are present in audited production files.
5. No coupling contamination indicators are found.
6. Dat file `stage10_outputs/stage10_2_hook_site_audit.dat` reports all required statuses as `1`.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
bash stage10_checks/run_stage10_2_hook_site_audit.sh
```
