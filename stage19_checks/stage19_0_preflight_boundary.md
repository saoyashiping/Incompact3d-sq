# Stage 19.0 — Stage 18 closure and Stage 19 preflight boundary

Stage 19.0 is a diagnostic-only preflight gate for entering Stage 19. It verifies
that Stage 18 is closed or, for a source-only archive, that the Stage 18.12
closure gate source evidence is present and sufficient to accept the already
completed Stage 18 closure without rerunning Stage 18.0--18.11.

Important behavior:

- Stage 19.0 does **not** rerun Stage 18.0--18.11 by default.
- Stage 19.0 does **not** require individual `stage18_outputs/*.dat` files when
  Stage 18 closure evidence or source-only Stage 18.12 closure-gate evidence is
  valid.
- Stage 19.0 remains preflight-only: no production structure state, no production
  hook, no production advance API, no RHS/IBM/DNS-core changes, no MPI, no DNS,
  no production restart/statistics/visualization I/O changes.
- Python byte-code syntax checks use a temporary `.pyc` path and do not write to
  `/dev/null`, which is required for Python 3.13 compatibility.

Expected verdict:

```text
STAGE 19.0 PREFLIGHT BOUNDARY VERDICT: PASS
STAGE 19.0 FINAL VERDICT: PASS
```
