# Stage 13.0 Preflight Closure Integrity

## Target

Stage 13.0 is a preflight integrity stage before real Stage 13 Eulerian force-density candidate work begins. It checks Stage 12 closure-file integrity and repairs the Stage 11 finalize chain in the Stage 9.7, Stage 9.8, and Stage 9.9 early-stop branches of `src/xcompact3d.f90`.

Stage 13.0 does not introduce Stage 13 physics.

## Mathematical and physical meaning

The no-fibre production DNS remains unchanged:

```text
f_fsi = 0
RHS_stage13.0 = RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

Stage 13.0 does not create an Eulerian force-density candidate and does not perform RHS injection:

```text
f_IBM_candidate is not created
RHS <- RHS + f_IBM_candidate is forbidden
```

There is no spreading, no feedback application to the fluid, no two-way coupling, and no fibre structure advance.

## Repair scope

The allowed repair scope is intentionally narrow:

- add missing `stage11_production_oneway_finalize()` calls in the Stage 9.7 / 9.8 / 9.9 early-stop finalize branches in `src/xcompact3d.f90`;
- preserve existing Stage 10 finalize calls;
- preserve existing Stage 12 finalize calls;
- create `stage12_checks/STAGE12_CLOSED.md` only if it is missing;
- add this Stage 13.0 gate and documentation.

No closed Stage 10, Stage 11, or Stage 12 scripts/modules are modified.

## What is intentionally not done

Stage 13.0 intentionally does not add:

- Stage 13 force-density candidate physics;
- a new xcompact3d physics hook;
- RHS modification;
- IBM spreading;
- feedback force application to fluid;
- two-way force density;
- fibre structure advance.

## Pass criteria

The preflight gate passes only if all of the following are true:

1. `stage12_checks/STAGE12_CLOSED.md` exists.
2. `xcompact3d` and all required Stage 11 / Stage 12 check targets through Stage 12.6 build.
3. Normal exit contains Stage 11 finalize.
4. Stage 9.7 early-stop contains Stage 11 finalize.
5. Stage 9.8 early-stop contains Stage 11 finalize.
6. Stage 9.9 early-stop contains Stage 11 finalize.
7. Stage 10 finalize remains present.
8. Stage 12 finalize remains present.
9. Static checks find no active Eulerian force-density, RHS-injection, IBM-spreading, feedback-application, two-way-force, or structure-advance signs in `src/xcompact3d.f90`.

The gate writes `stage13_outputs/stage13_0_preflight_closure_integrity.dat` and prints:

```text
STAGE 13.0 PREFLIGHT CLOSURE INTEGRITY VERDICT: PASS
STAGE 13.0 FINAL VERDICT: PASS
```

Failures print explicit non-empty reasons, including branch-specific Stage 11 finalize failures.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
STAGE13_0_RUN_STAGE12_CLOSURE=0 \
bash stage13_checks/run_stage13_0_preflight_closure_integrity.sh
```
