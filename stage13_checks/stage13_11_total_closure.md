# Stage 13.11 total closure gate

## Target

Stage 13.11 is the total closure gate for the Stage 13 production Eulerian
force-density candidate path. It closes Stage 13 only after every required build,
smoke, invariance, parallel-consistency, I/O-compatibility, contamination-audit,
and diagnostic check passes.

## Required closure evidence

The gate requires all of the following:

1. automatic CMake configure when `build_stage9` is absent;
2. `xcompact3d` build PASS;
3. all Stage 11 check targets build PASS;
4. all Stage 12 check targets through Stage 12.6 build PASS;
5. all Stage 13.1–13.6 check targets build PASS;
6. Stage 13.6 production force-density candidate hook smoke PASS;
7. Stage 13.7 np=1 force-density no-contamination invariance PASS;
8. Stage 13.8 np=1/2/4 parallel force-density consistency PASS;
9. Stage 13.9 restart / stats / visualization / coarse I/O compatibility PASS;
10. Stage 13.10 RHS injection / production IBM forcing / structure contamination audit PASS;
11. `stage13_outputs/fibre_stage13_6_production_force_density_candidate.dat` exists;
12. force-density candidate computed and finite statuses are PASS;
13. force-density norm and integrated-force diagnostics are finite and conservative;
14. spreading input sign and wrong-sign rejection statuses are PASS;
15. field and RHS modified statuses remain zero;
16. no RHS injection, production IBM forcing, feedback application, two-way force, or structure advance is reported.

## Closure file policy

`stage13_checks/STAGE13_CLOSED.md` is generated only after the full Stage 13.11
closure gate passes. A failed or incomplete Stage 13.11 run must not create the
closure marker.

## Closed-stage protection

Stage 13.11 does not modify Stage 10, Stage 11, Stage 12, Stage 13.0–13.10,
production solver, RHS/projection/Poisson/RK3/channel-forcing, production IBM, or
structure-solver source files.

## Manual command

```sh
DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4 \
BUILD_DIR=build_stage9 \
MPIEXEC=mpirun \
MPIEXEC_FLAGS="--mca btl self,vader,tcp" \
bash stage13_checks/run_stage13_11_total_closure.sh
```

## Expected final verdict

A complete successful closure run prints:

```text
STAGE 13.11 FINAL VERDICT: PASS
```
