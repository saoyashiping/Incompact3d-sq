# Stage 11.9 RHS / Coupling Contamination Audit

## Stage 11.9 target

Perform a strict RHS / coupling contamination audit for Stage 11 one-way sampling.

## Mathematical / physical meaning

- `U_f = I_h[u](X_f)` is read-only sampling.
- `f_fsi = 0`.
- `RHS_stage11.9 = RHS_stage10 = RHS_stage9`.
- no spreading
- no feedback
- no two-way force
- no structure advance

## Static audit scope

- `src/xcompact3d.f90`
- `src/fibre_stage11_production_oneway_hook.f90`

## Runtime evidence

- Stage 9.9 deterministic no-fibre smoke with Stage 11 hook enabled.
- Hook diagnostics from `stage11_outputs/fibre_stage11_5_production_oneway_hook.dat`.

## Hook diagnostic evidence

- sampled velocity finite
- sample performed
- no field modification
- no RHS modification
- no IBM / feedback / two-way / structure advance

## What is intentionally not done

- no RHS injection
- no IBM spreading
- no feedback force
- no two-way force
- no fibre structure advance
- no Stage 10 script reruns
- no closed-stage edits

## Pass criteria

Pass requires:

- required build targets succeed,
- static audit finds no forbidden active imports/calls/writes,
- Stage 9.9 runtime smoke passes under Stage 11 hook,
- Stage 11.5 hook diagnostics are present and all required keys match expected values.

## Manual command

```bash
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE11_9_RUN_STAGE11_8=0 \
bash stage11_checks/run_stage11_9_rhs_coupling_contamination_audit.sh
```
