# Production Recovery

This directory is the entry point for Production Recovery R0-R12.

## Current entry status

R0 records the status correction: Stage 22 source-only closure exists, but production DNS-FSI closure is not complete. The authoritative project status is maintained in `../PRODUCTION_RECOVERY_STATUS.md`.

## Stage dependencies

Production Recovery is sequential. Later stages depend on the evidence and production integration established by earlier stages:

1. **R0 status correction and evidence-boundary reset**: documentation-only correction; no source solver changes and no build/DNS/MPI validation.
2. **R1 baseline build/run cleanup**: establish a clean baseline compile/run path before physics recovery.
3. **R2 production fibre state**: introduce production-managed fibre state.
4. **R3 real Xcompact3D grid adapter**: connect fibre/IBM logic to the real solver grid.
5. **R4 production IBM interpolation**: implement production interpolation from Eulerian flow to Lagrangian fibre state.
6. **R5 production structure solver**: advance fibre structure in the production loop.
7. **R6 IBM spreading into real RHS**: spread IBM forces into the real RHS rather than a placeholder or diagnostic path.
8. **R7 two-way FSI closure**: close the production two-way coupling loop.
9. **R8 wall contact**: integrate wall-contact coupling into the production path.
10. **R9 fibre-fibre collision**: integrate fibre-fibre collision coupling into the production path.
11. **R10 restart/statistics/visualization**: persist, report, and visualize production fibre state.
12. **R11 real MPI np=1/2/4 consistency**: validate real MPI consistency with actual runs.
13. **R12 paper-grade validation matrix**: complete publication-grade validation evidence.

## Evidence rule

Source-only checks and historical closure files are not substitutes for real production build, MPI, DNS, or FSI validation. Each recovery stage must clearly label its evidence boundary.
