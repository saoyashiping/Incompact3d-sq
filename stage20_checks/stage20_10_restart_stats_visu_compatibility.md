# Stage 20.10: restart/statistics/visualization compatibility for active coupling

Stage 20.10 is a diagnostic-only helper-level compatibility audit. It summarizes the Stage 20 active-coupling diagnostic metadata that restart, statistics, and visualization reporting would need, while keeping every payload strictly inside the Stage 20.10 audit output.

This stage writes only `stage20_outputs/fibre_stage20_10_restart_stats_visu_compatibility.dat`. It does not write production restart files, checkpoint files, statistics files, visualization files, flow-field files, or production I/O schemas.

It does not run actual MPI, production DNS, production validation, production RHS updates, Stage 14 RHS injection, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, restart writers, statistics writers, or visualization writers.

## Source-only and no-rerun policy

Stage 20.10 accepts Stage 20.9 PASS evidence when present and preserves source-only acceptance for prior Stage 20 closure behavior. It does not require old closure files, does not require all old stage outputs, and does not rerun Stage 10-19 or Stage 20.0-20.9.

## Audit-only summary groups

The audit constructs helper-local metadata summaries only:

* coupling gate summary;
* lambda summary;
* force summary;
* structure summary;
* RHS diagnostic summary;
* np=1/2/4 parallel consistency summary;
* restart compatibility markers;
* statistics compatibility markers;
* visualization compatibility markers.

These summaries are compatibility payload candidates, not production I/O payloads.

## Safety boundary

`STAGE20_10_RESTART_AUDIT_ONLY`, `STAGE20_10_STATISTICS_AUDIT_ONLY`, and `STAGE20_10_VISUALIZATION_AUDIT_ONLY` must remain enabled. Production restart, statistics, and visualization I/O gates must remain disabled. Compatibility payloads are written only to the Stage 20.10 `.dat` audit output.

## Next stage

Stage 20.11: Stage 20 total contamination audit and closure.

## Manual command

```bash
stage20_checks/run_stage20_10_restart_stats_visu_compatibility.sh
```

## Expected PASS evidence

```text
STAGE 20.10 RESTART STATS VISU COMPATIBILITY VERDICT: PASS
STAGE 20.10 FINAL VERDICT: PASS
```
