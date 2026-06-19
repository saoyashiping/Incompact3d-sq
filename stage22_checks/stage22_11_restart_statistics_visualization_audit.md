# Stage 22.11 restart/statistics/visualization production-readiness audit

Stage 22.11 is the final functional audit before Stage 22.12 total closure.  It is a controlled restart/statistics/visualization production-readiness audit only: no new physics, no mesh refinement, no G2/G3, np=1 primary, np=2 optional only if safe and isolated, np=4 disabled by default, no production source modification, and no restart/statistics/visualization schema modification.

## Source-only and prior-stage acceptance

Stage 22.11 accepts Stage 20 closure, Stage 21 closure, and Stage 22.0 through Stage 22.10 PASS evidence through available committed artifacts or source-only closure acceptance.  Missing old outputs are allowed and previous stages are not rerun.

## Controlled G1 I/O audit case

* grid label: G1
* Nx = 32
* Ny = 33
* Nz = 32
* channel case
* channel_height = 1.0
* y_min = 0.0
* y_max = 1.0
* dt = 1.0e-4
* fresh_steps = 100
* restart_steps = 100
* total_equivalent_steps = 200
* final_time = 0.02
* CFL_max_limit = 0.3
* np=1 primary
* np=2 secondary optional documented
* np=4 disabled by default
* G2 disabled
* G3 disabled
* n_fibre = 1
* n_point_per_fibre = 64
* lambda_fsi = 1.0e-4
* lambda_contact = 0.0
* wall contact disabled in default baseline
* fibre-fibre collision disabled in default baseline

## Restart audit requirements

The fresh segment writes restart output using the existing restart path only, without restart schema modification.  The restart continuation reads that restart output, starts at the expected step/time, runs the continuation segment, and writes continuation evidence without schema modification.  A continuous comparison reference is optional; if disabled, internal restart integrity and finite-output checks are sufficient.

## Statistics audit requirements

The statistics path is audited using the existing statistics output path only.  Statistics output may be active and readable, or explicitly documented as inactive without schema modification.  Active statistics must have monotonic timestamps, finite fields, no NaN/Inf, no unexpected column drift, and no FSI/contact field contamination outside the existing verified output contract.

## Visualization audit requirements

The visualization path is audited using the existing visualization output path only.  Visualization output may be active and readable, or explicitly documented as inactive without schema modification.  Active visualization must have finite metadata, consistent dimensions, finite fields, no NaN/Inf, no file corruption, no unsupported format conversion, and no external post-processing dependency.

## Production isolation requirements

No production source, CMake, DNS/RHS/IBM, pressure projection, Poisson, RK3/channel forcing, restart schema, statistics schema, visualization schema, uncontrolled multifibre, hidden production I/O hook, or Stage 22.12 final closure file may be modified or created in Stage 22.11.

## Safe defaults

* STAGE22_11_ENABLE=1
* STAGE22_11_RESTART_STATISTICS_VISUALIZATION_AUDIT_ENABLE=1
* STAGE22_11_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_11_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_11_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_11_FAIL_CLOSED=1
* STAGE22_11_CONTROLLED_IO_AUDIT_ALLOWED=1
* STAGE22_11_BUILD_ALLOWED=1
* STAGE22_11_PRODUCTION_DNS_ALLOWED=1
* STAGE22_11_ACTUAL_MPI_ALLOWED=1
* STAGE22_11_NP_PRIMARY=1
* STAGE22_11_NP2_OPTIONAL_ENABLE=1
* STAGE22_11_NP4_ENABLE=0
* STAGE22_11_G1_ENABLE=1
* STAGE22_11_G2_ENABLE=0
* STAGE22_11_G3_ENABLE=0
* STAGE22_11_FRESH_STEPS=100
* STAGE22_11_RESTART_STEPS=100
* STAGE22_11_TOTAL_EQUIVALENT_STEPS=200
* STAGE22_11_FINAL_TIME=0.02
* STAGE22_11_LAMBDA_FSI=1.0e-4
* STAGE22_11_LAMBDA_CONTACT=0.0
* STAGE22_11_WALL_CONTACT_FORCE_ENABLE=0
* STAGE22_11_COLLISION_FORCE_ENABLE=0
* STAGE22_11_CONTINUOUS_REFERENCE_OPTIONAL_ENABLE=0
* STAGE22_11_RESTART_TEST_ENABLE=1
* STAGE22_11_STATISTICS_AUDIT_ENABLE=1
* STAGE22_11_VISUALIZATION_AUDIT_ENABLE=1
* STAGE22_11_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_11_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_11_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_11_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED=0
* STAGE22_11_STAGE22_12_CLOSURE_FILE_CREATION_ALLOWED=0

## Required PASS evidence

The audit writes `stage22_outputs/fibre_stage22_11_restart_statistics_visualization_audit.dat` and must contain:

```text
STAGE 22.11 RESTART STATISTICS VISUALIZATION AUDIT VERDICT: PASS
STAGE 22.11 FINAL VERDICT: PASS
```

## Next stage

Stage 22.12: final total closure.
