# Stage 10.7 IO/Restart/Stats/Visu No-Op Validation

## Stage 10.7 target

Validate restart, statistics, visualization, and coarse I/O behavior under Stage 10 no-op hook environment.

## Mathematical / physical meaning

- `f_fsi = 0`
- `RHS_stage10 = RHS_stage9`
- no-op hook must not alter restart or I/O behavior

## Reused gates

- Stage 9.7 stats/visu/coarse I/O smoke
- Stage 9.8 restart I/O regression

## Explicit old-logic warning

- do not run Stage 10.2 / Stage 10.3 by default
- do not forbid valid guarded hook calls
- do not broad-grep `ibm` / `structure` / `fibre_`
- do not classify `no_*` diagnostic fields as contamination

## Intentionally not done

- no real IBM
- no real fibre force
- no RHS injection
- no structure advance
- no two-way coupling

## Pass criteria

- required build targets pass
- Stage 9.7 passes under Stage 10 no-op hook environment
- Stage 9.8 passes under Stage 10 no-op hook environment
- Stage 10.3 hook dat safety keys remain valid
- dat file `stage10_outputs/stage10_7_io_restart_stats_visu_noop.dat` reports final status `1`

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE10_7_RUN_PREREQS=0 \
bash stage10_checks/run_stage10_7_io_restart_stats_visu_noop.sh
```
