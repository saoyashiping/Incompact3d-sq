# Stage 21.8: contact candidate persistence audit

Stage 21.8 is a diagnostic-only persistence audit for Stage 21 contact/collision candidate metadata. It verifies that gap, warning/risk/fail, registry, ownership, deterministic ordering, and reduction-ready summary metadata can be serialized to Stage 21 diagnostic output, reloaded, reconstructed, hashed, and compared deterministically.

## Source-only and no-rerun policy

Stage 21.8 accepts Stage 21.7 PASS evidence when available, but it does not rerun Stage 21.7 or any earlier stage. Source-only archives are accepted by default with `STAGE21_8_ALLOW_SOURCE_ONLY_ARCHIVE=1` and missing old outputs are allowed with `STAGE21_8_ALLOW_MISSING_OLD_OUTPUTS=1`.

## Persistence scope

The persisted diagnostic metadata chain is:

```text
gap metadata
warning/risk/fail metadata
candidate registry metadata
ownership metadata
deterministic ordering metadata
reduction-ready global summary
```

## Persistence rules

The helper creates a deterministic text serialization using sorted canonical candidate order. The reloaded and reconstructed payloads must preserve candidate count, candidate IDs, canonical keys, canonical sort keys, global order indices, gap metadata, risk metadata, registry metadata, owner ranks for `np=1,2,4`, deterministic ordering metadata, and disabled diagnostic gate values.

Hash equality is required across the serialization, reload, reconstruction, and final roundtrip hashes. The diagnostic schema is named `stage21_contact_candidate_metadata` with schema version `1` by default.

## Production I/O isolation

Stage 21.8 may write only `stage21_outputs/fibre_stage21_8_contact_candidate_persistence_audit.dat`. It does not modify production restart, checkpoint, statistics, visualization, flow-field, DNS, RHS, IBM, CMake, or source files. Restart/statistics/visualization payloads are documented only as diagnostic candidates and are explicitly marked as not production payloads.

## Safety boundary

No contact force, collision force, wall force, force application, structure advance update, RHS coupling, Stage 14 RHS injection, MPI execution, production DNS execution, production multifibre activation, production hook, or production I/O schema modification is introduced.

## Manual command

```bash
stage21_checks/run_stage21_8_contact_candidate_persistence_audit.sh
```

## Expected PASS evidence

The wrapper writes `stage21_outputs/fibre_stage21_8_contact_candidate_persistence_audit.dat` and prints:

```text
STAGE 21.8 CONTACT CANDIDATE PERSISTENCE AUDIT VERDICT: PASS
STAGE 21.8 FINAL VERDICT: PASS
```

## Next stage

Stage 21.9: contact diagnostic integration
