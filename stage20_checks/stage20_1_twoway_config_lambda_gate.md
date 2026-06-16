# Stage 20.1 two-way coupling config and lambda gate

Stage 20.1 is a definition-only, diagnostic-only configuration boundary for
future real two-way FSI activation. It defines and audits Stage 20 conceptual
gates, but it does not activate production coupling or runtime hooks.

## Source-only and previous-stage policy

Stage 20.1 accepts Stage 20.0 PASS evidence when present and preserves Stage
20.0 source-only acceptance behavior. It does not rerun Stage 10-19 or Stage
20.0, does not require all old closure files, and does not require all old stage
output directories in source-only archives.

## Conceptual gates and safe defaults

- `STAGE20_ENABLE`: true only for Stage 20 diagnostic/config audit
- `STAGE20_TWOWAY_COUPLING_ENABLE`: false
- `STAGE20_FLUID_TO_STRUCTURE_ENABLE`: false
- `STAGE20_STRUCTURE_TO_FLUID_ENABLE`: false
- `STAGE20_RHS_COUPLING_ENABLE`: false
- `STAGE20_LAMBDA_COUPLING`: 0.0
- `STAGE20_DIAGNOSTIC_ONLY`: true
- `STAGE20_FAIL_CLOSED`: true
- `STAGE20_SINGLE_FIBRE_ONLY`: true
- `STAGE20_CONTACT_ENABLE`: false
- `STAGE20_COLLISION_ENABLE`: false
- `STAGE20_MULTIFIBRE_ENABLE`: false

## Consistency rules

- If two-way coupling is disabled, fluid-to-structure, structure-to-fluid, and
  RHS coupling remain disabled.
- If lambda coupling is zero, effective coupling is zero.
- If RHS coupling is disabled, Stage 14 RHS injection is not callable.
- Diagnostic-only and fail-closed mode prevent production runtime activation.
- Contact, collision, and multifibre gates are disabled by default and introduce
  no wall contact force, fibre-fibre collision force, or production multifibre
  logic.

## Safety boundary

Stage 20.1 does not modify production Fortran, CMake, IBM, DNS-core,
projection, Poisson, RK3/channel forcing, restart/statistics/visualization I/O,
or any Stage 10-19 or Stage 20.0 file. It does not run MPI, DNS, previous
stages, builds, or production validation.

## Next stage

Stage 20.2: fluid-to-structure force input adapter.

## Manual command

```bash
stage20_checks/run_stage20_1_twoway_config_lambda_gate.sh
```

Expected output includes:

```text
STAGE 20.1 TWOWAY CONFIG LAMBDA GATE VERDICT: PASS
STAGE 20.1 FINAL VERDICT: PASS
```
