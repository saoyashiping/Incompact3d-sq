# Stage 18.12: total contamination audit and closure

Stage 18.12 is the final diagnostic-only total audit for Stage 18.  It performs
read-only evidence inspection for Stage 18.0--18.11, checks syntax/compile health
of Stage 18 helpers and wrappers, verifies false-positive-safe preservation
signals, and creates `stage18_checks/STAGE18_CLOSED.md` only after every status
field passes.

## Scope

**TOTAL AUDIT** means read-only Stage 18 evidence inspection and targeted
source/path contamination checks.  **CLOSURE** means writing
`stage18_checks/STAGE18_CLOSED.md` only after all checks pass.  **HELPER-LOCAL
OUTPUT** means files under `stage18_outputs` only.  **PRODUCTION I/O** means
runtime restart/statistics/visu output or Fortran production I/O modifications;
Stage 18.12 does not modify those.  **PRODUCTION PHYSICS** means runtime `X/V/A`
updates or RHS/IBM/DNS-core coupling; Stage 18.12 does not activate those.

## Evidence checked

Stage 18.12 verifies the presence of Stage 18.0--18.11 wrappers, helpers,
documentation, and existing `.dat` evidence files.  By default it requires prior
outputs and does not rerun prior stages:

```text
STAGE18_12_REQUIRE_PRIOR_OUTPUTS=1
STAGE18_12_RERUN_PRIOR_STAGES=0
```

Stage 18.11 helper-local JSON snapshots under `stage18_outputs` are accepted as
helper-local evidence and are not production restart/statistics/visu I/O.

## Closure marker

If and only if all Stage 18.12 status checks pass, the helper creates
`stage18_checks/STAGE18_CLOSED.md` containing:

* Stage 18 closed
* generated date/time when available
* list of closed Stage 18.0--18.12
* a diagnostic-only single-fibre physical structure dynamics statement
* a no-production RHS/IBM/DNS-core/stats/visu/restart I/O contamination statement
* a real contact/collision/multifibre production logic exclusion statement
* a note that the next stage should be Stage 19 or another user-defined milestone

## False-positive-safe policy

Stage 18.12 continues the corrected Stage 18.11 / Stage 18.10 / Stage 18.9 /
Stage 18.8 / Stage 18.7 / Stage 18.6 / Stage 18.5 / Stage 18.0 / Stage 17.11 /
Stage 17.10 / Stage 17.6 / Stage 16 audit pattern: no broad scans, no Markdown as
activation evidence, no mandatory `rg`, source-only archives accepted, helper
outputs under `stage18_outputs` are not production I/O, MPI compatibility text is
not MPI execution, and only `*_status` fields determine `final_status`.

## Manual command and expected PASS evidence

Run:

```bash
stage18_checks/run_stage18_12_total_contamination_audit_closure.sh
```

Expected terminal evidence when all prior Stage 18 outputs are present and PASS:

```text
STAGE 18.12 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: PASS
STAGE 18.12 FINAL VERDICT: PASS
```

The wrapper writes:

```text
stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat
```

and creates `stage18_checks/STAGE18_CLOSED.md` only after every check passes.
