# Stage 19.0 preflight boundary

Stage 19.0 is a diagnostic-only boundary gate between the fully closed Stage 18
physical-structure helper sequence and future Stage 19 production-side physical
structure integration work.

## Required evidence inspected read-only

The helper checks that `stage18_checks/STAGE18_CLOSED.md` exists and has valid
closure content.  Stage 18 is accepted as closed when that marker is valid and
either the Stage 18.12 closure output reports PASS or the closure marker itself
contains strong text proving Stage 18.0 through Stage 18.12 were closed after the
total audit.

Individual Stage 18.0 through Stage 18.11 helper outputs are informational in
Stage 19.0.  They are reported as `PASS` when present and passing, or as
`ACCEPTED_BY_STAGE18_CLOSURE` when Stage 18 closure evidence supersedes them.
Stage 19.0 must not force users to rerun Stage 18.0 through Stage 18.11 checks.
Helper-local Stage 18.11 JSON snapshots under `stage18_outputs` are explicitly
accepted as helper evidence, not production restart/statistics/visualization I/O.

## Boundary definitions

* **PREFLIGHT BOUNDARY** means read-only evidence inspection and Stage 19
  boundary definition.
* **PRODUCTION STRUCTURE INTEGRATION** means actual production X/V/A state,
  production structure hook, production advance API, or production commit. Stage
  19.0 does not add this yet.
* **HELPER OUTPUT** means `stage19_outputs` only.
* **PRODUCTION I/O** means runtime restart/statistics/visualization output or
  production Fortran I/O. Stage 19.0 does not modify this.
* **PRODUCTION FSI COUPLING** means RHS/IBM/DNS-core coupling. Stage 19.0 does
  not activate this.

## Future Stage 19 scope, not implemented in Stage 19.0

Stage 19 will later introduce production-side physical structure state,
production-side candidate structure advance API, production-side structure hook,
no-op invariance, controlled production-side single-fibre response, parallel
consistency, and restart/I/O boundary checks.

## Explicit Stage 19.0 prohibitions

Stage 19.0 forbids production X/V/A state creation, production structure hook
insertion, production structure advance activation, Stage 14 RHS injection,
force spreading to Eulerian RHS, IBM modification, DNS-core modification,
pressure projection / Poisson / RK3 channel forcing modification, production
restart/statistics/visualization I/O modification, MPI execution, production DNS
execution, real wall contact force, real fibre-fibre collision force, penalty /
repulsive / lubrication / friction / adhesion / contact damping force,
collision-induced RHS or collision-induced structure update, and production
multifibre logic.

## Manual command

```bash
stage19_checks/run_stage19_0_preflight_boundary.sh
```

Expected PASS evidence in a repository with preserved Stage 18 closure marker and Stage 18.12 closure evidence:

```text
STAGE 19.0 PREFLIGHT BOUNDARY VERDICT: PASS
STAGE 19.0 FINAL VERDICT: PASS
final_status PASS
```
