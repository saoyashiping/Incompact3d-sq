# Stage 14.5 — Guarded Production RHS Hook Skeleton

## Target

Stage 14.5 installs the first guarded production RHS-injection hook skeleton for the Stage 14 pathway.

The production hook is connected only under explicit Stage 14 configuration and is validated only with:

```text
lambda_14 = 0
```

No small-lambda production response is tested in Stage 14.5. Small nonzero response remains reserved for Stage 14.7.

## Mathematical / physical meaning

Stage 13 provides the Eulerian force-density candidate:

```text
f_i_cand = S_h[F_fibre_to_fluid_cand]
```

Stage 14.1 provides the diagnostic RHS accumulator:

```text
f_rhs_stage14 = lambda_14 * f_i_cand
```

Stage 14.2 verifies the controlled standalone formula:

```text
RHS_new = RHS_old + f_rhs_stage14
```

Stage 14.3 verifies the fibre-on-fluid sign and lambda scaling. Stage 14.4 verifies the planned RK-substep timing.

Stage 14.5 installs the guarded production hook skeleton. For the Stage 14.5 validation smoke:

```text
lambda_14 = 0
RHS_new = RHS_old
f_fsi = 0
```

The hook records that it was called and that the computed Stage 14 increment is zero, but it does not alter the RHS for Stage 14.5 validation.

## Production policy

Stage 14.5 follows these production-safety rules:

- default off;
- requires explicit Stage 14 request;
- requires RHS injection flag enabled;
- requires Stage 13 dependency to be requested;
- validates production behavior only with `lambda_14=0`;
- blocks nonzero lambda response in this stage;
- writes diagnostics showing hook activity and zero RHS increment.

## What is intentionally not done

Stage 14.5 intentionally does not:

- test small-lambda production response;
- modify pressure, projection, Poisson, RK3, or channel forcing source;
- call production IBM forcing application;
- apply feedback force beyond already-closed Stage 13 diagnostics;
- activate real two-way force coupling;
- advance the fibre structure.

## Pass criteria

The Stage 14.5 gate passes only if:

- all required standalone hook-check keys are `1`;
- the production smoke reports `STAGE 9.9 FINAL VERDICT: PASS`;
- the production hook diagnostics report the hook initialized and applied;
- `lambda_14=0` is confirmed;
- nonzero lambda is marked blocked for this stage;
- RHS signature delta, increment L2, and increment max-abs are all `<= 1.0e-12`;
- pressure/projection/Poisson/RK3/channel-forcing and structure-advance guard statuses are all `1`.

The production hook writes:

```text
stage14_outputs/fibre_stage14_5_production_rhs_hook.dat
```

The standalone check writes:

```text
stage14_outputs/fibre_stage14_5_production_rhs_hook_check.dat
```

The gate writes:

```text
stage14_outputs/stage14_5_production_rhs_hook_gate.dat
```

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE14_5_RUN_STAGE14_4=0 \
bash stage14_checks/run_stage14_5_production_rhs_hook.sh
```
