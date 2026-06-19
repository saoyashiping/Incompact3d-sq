# Production Recovery Status

## Authoritative current status

**Stage 22 source-only closure exists; production DNS-FSI closure is not complete.**

This document supersedes unconditional final-closure language in historical Stage 22 closure artifacts. Stage 22 artifacts are preserved as historical source-only audit evidence, but they must not be used as proof of production-ready DNS-FSI. True production closure requires completion of Production Recovery R1-R12.

## R0 scope statement

R0 is documentation/status correction only.

No source solver modification is performed in R0.

No DNS/MPI/build validation is performed in R0.

## Current source facts

* The `xcompact3d` production target does not contain production Stage 20-22 FSI closure.
* The Stage 15 hook is diagnostic/no-op or non-advancing unless later production code is added.
* The Stage 14 RHS injection is not a physical IBM spreading implementation.
* Existing Stage 22 checks are source-only and must not replace actual build/MPI/DNS evidence.

## Evidence boundary

Historical Stage 22 evidence may document helper logic, static checks, synthetic metrics, and source-only audit conclusions. It does not certify:

* production DNS-FSI closure,
* real MPI `np=1/2/4` execution,
* real build/run validation,
* production wall-contact coupling,
* production fibre-fibre collision coupling, or
* a complete production restart/statistics/visualization loop for fibre state.

## Production Recovery route

Production DNS-FSI closure must be re-established through the following staged route:

1. **R1 baseline build/run cleanup**
2. **R2 production fibre state**
3. **R3 real Xcompact3D grid adapter**
4. **R4 production IBM interpolation**
5. **R5 production structure solver**
6. **R6 IBM spreading into real RHS**
7. **R7 two-way FSI closure**
8. **R8 wall contact**
9. **R9 fibre-fibre collision**
10. **R10 restart/statistics/visualization**
11. **R11 real MPI np=1/2/4 consistency**
12. **R12 paper-grade validation matrix**

## Status-use rule

When this document conflicts with historical closure files, this document is authoritative. Historical files must be read as source-only closure artifacts unless a later Production Recovery stage provides real build/MPI/DNS evidence.
