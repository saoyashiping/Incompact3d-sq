# Stage 21.9: contact diagnostic integration

Stage 21.9 is a diagnostic-only integration layer for Stage 21 wall-contact and fibre-fibre collision safety diagnostics. It integrates wall gap diagnostics, fibre-fibre distance/gap diagnostics, warning and fail-closed states, the contact candidate registry, pair ownership, deterministic ordering, metadata consistency, persistence/reload readiness, and reduction-ready `np=1,2,4` summaries into one diagnostic report.

## Source-only and no-rerun policy

Stage 21.9 accepts Stage 21.8 PASS evidence when present, does not rerun Stage 21.8 or any earlier stage, and accepts source-only archives with missing old outputs by default.

## Integrated diagnostic chain

```text
Stage 21.1 wall distance and signed-gap audit
  -> Stage 21.2 fibre-fibre point/segment distance audit
  -> Stage 21.3 near-contact warning and fail-closed gate
  -> Stage 21.4 contact candidate registry
  -> Stage 21.5 contact pair ownership audit
  -> Stage 21.6 deterministic pair ordering
  -> Stage 21.7 contact metadata consistency
  -> Stage 21.8 contact candidate persistence audit
  -> Stage 21.9 integrated diagnostic summary
```

## Integrated groups

The report contains wall diagnostics, fibre-fibre diagnostics, warning/fail-closed counts, candidate registry counts, ownership validity, ordering determinism, metadata consistency, persistence/reload hash equality, production isolation, and next-stage readiness.

## Safety boundary

Stage 21.9 does not compute contact force, collision force, or wall force; does not apply contact/collision force; does not modify structure advance or RHS; does not call Stage 14 RHS injection; does not run MPI or production DNS; does not activate production multifibre logic; and does not modify or write production restart, statistics, visualization, checkpoint, flow-field, source, CMake, DNS, RHS, or IBM files.

## Safe defaults

The helper defaults are diagnostic-only and fail-closed with `STAGE21_9_ENABLE=1`, `STAGE21_9_CONTACT_DIAGNOSTIC_INTEGRATION_ENABLE=1`, `STAGE21_9_DIAGNOSTIC_ONLY=1`, `STAGE21_9_FAIL_CLOSED=1`, `STAGE21_9_NP_VALUES=1,2,4`, and all force, RHS, structure, Stage 14 RHS injection, production DNS, MPI, production multifibre, restart I/O, statistics I/O, and visualization I/O gates disabled.

## Manual command

```bash
stage21_checks/run_stage21_9_contact_diagnostic_integration.sh
```

## Expected PASS evidence

The wrapper writes `stage21_outputs/fibre_stage21_9_contact_diagnostic_integration.dat` and prints:

```text
STAGE 21.9 CONTACT DIAGNOSTIC INTEGRATION VERDICT: PASS
STAGE 21.9 FINAL VERDICT: PASS
```

## Next stage

Stage 21.10: collision-force-disabled proof
