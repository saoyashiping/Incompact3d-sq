# R1 Known Warnings

This file records R1 residual warnings without rerunning R1 and without modifying the real R1 run logs.

## Git metadata warning

R1 evidence may include transient `git status --short` metadata from files that were being created during evidence capture. This is a documentation/evidence-timing warning only.

## `airfoils.f90` warning

Any `airfoils.f90` build warning observed in downstream environments should be treated as a baseline build warning to triage separately. R1 warnings do not certify or invalidate fibre/FSI physics.

## OpenMPI warning

OpenMPI launcher warnings, including root-execution or environment warning banners, should be preserved in raw run logs when present. They must not be rewritten as PASS evidence.
