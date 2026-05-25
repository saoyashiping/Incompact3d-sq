# Stage 10.6 RHS Contamination Audit

## Stage 10.6 target

Stage 10.6 performs a production RHS contamination audit to prove that no FSI/IBM/fibre force has been injected into the Navier-Stokes RHS.

## Mathematical / physical meaning

- `f_fsi = 0`
- `RHS_stage10 = RHS_stage9`

## Allowed production hook state

Guarded Stage 10 hook calls in `xcompact3d.f90` are allowed after Stage 10.3:

- `if (stage10_reg) call stage10_hook_...`

These guarded calls are control-flow hooks and are **not** RHS force injection.

## Explicit old-logic warning

- Do not forbid valid guarded hook calls.
- Do not classify `no_ibm_call_status`, `no_structure_advance_status`, `no_fibre_state_status` as contamination.
- Do not broad-grep `ibm`, `structure`, or `fibre_`.
- Audit only active Fortran `use` / `call` / assignment statements after stripping comments.

## Files audited

- `src/xcompact3d.f90`
- `src/navier.f90`
- `src/time_integrators.f90`
- `src/derive.f90`
- `src/poisson.f90`
- `src/Case-Channel.f90`

## Intentionally not done

- no real IBM
- no real fibre force
- no RHS injection
- no structure advance
- no two-way coupling

## Pass criteria

- required build targets pass
- no `stage10_hook_*` calls in RHS/projection/Poisson/channel files
- guarded hook call pattern remains present in `xcompact3d.f90`
- no active contamination import/call/assignment patterns in audited RHS-related files
- dat file `stage10_outputs/stage10_6_rhs_contamination_audit.dat` reports final status `1`

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
STAGE10_6_RUN_PREREQS=0 \
bash stage10_checks/run_stage10_6_rhs_contamination_audit.sh
```
