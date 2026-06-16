# Stage 21.4: contact candidate registry

Stage 21.4 creates a diagnostic-only registry/schema for future wall-contact and fibre-fibre collision candidates. It stores metadata derived from existing helper-local distance, signed-gap, and warning logic, but it does not compute or apply real contact force, collision force, wall force, structure advance updates, RHS coupling, MPI, DNS, production hooks, or production multifibre logic.

Stage 21.4 accepts Stage 21.3 PASS evidence when present and preserves source-only closure acceptance. Missing old outputs and closure files are allowed when source-only acceptance is enabled, and no previous stage is rerun.

## Registry schema

Each candidate record contains:

* `candidate_id`
* `candidate_type`: `wall_lower`, `wall_upper`, or `fibre_fibre`
* `fibre_i`, `fibre_j`
* `point_i`, `point_j`
* `segment_i`, `segment_j`
* `gap_value`
* `penetration_depth`
* `warning_trigger`
* `fail_closed_trigger`
* `risk_level`
* `risk_label`
* `nearest_wall_side`
* `closest_pair_valid`
* `candidate_active`
* `force_computation_allowed`
* `force_application_allowed`
* `diagnostic_only`

## Validation policy

The registry verifies unique candidate IDs, valid candidate types, valid wall nearest-side metadata, ordered fibre-fibre pairs, risk-label/risk-level consistency, warning/fail flag consistency, diagnostic-only operation, disabled force computation/application, disabled structure advance, disabled RHS coupling, disabled DNS/MPI, and no closed-stage/source/CMake modification.

## Output

```text
stage21_outputs/fibre_stage21_4_contact_candidate_registry.dat
```

## Manual command

```bash
stage21_checks/run_stage21_4_contact_candidate_registry.sh
```

## Expected PASS evidence

```text
STAGE 21.4 CONTACT CANDIDATE REGISTRY VERDICT: PASS
STAGE 21.4 FINAL VERDICT: PASS
```

## Next stage

Stage 21.5: contact pair ownership audit
