# Stage 21.7: contact metadata consistency

Stage 21.7 is a diagnostic-only metadata consistency audit for the Stage 21 wall-contact / fibre-fibre collision safety chain. It verifies that helper-local gap metadata, warning and fail-closed metadata, candidate registry metadata, ownership metadata, deterministic ordering metadata, and reduction-ready summary metadata agree with each other.

## Source-only and no-rerun policy

This stage accepts Stage 21.6 PASS evidence when available, but it does not rerun Stage 21.6 or any earlier stage. Source-only archives are accepted by default through `STAGE21_7_ALLOW_SOURCE_ONLY_ARCHIVE=1` and `STAGE21_7_ALLOW_MISSING_OLD_OUTPUTS=1`.

## Metadata chain

```text
distance/gap metadata
  -> warning/risk/fail metadata
  -> candidate registry metadata
  -> ownership metadata
  -> deterministic ordering metadata
  -> reduction-ready global summary
```

## Consistency rules

For every diagnostic candidate, penetration metadata must satisfy `penetration_depth = max(0, -gap_value)`. Risk metadata follows the fail-closed classification: `SAFE` has risk level `0`, `NEAR_CONTACT` has risk level `1`, `OVERLAP` has risk level `2`, and `FAIL_CLOSED` has risk level `3`.

Allowed candidate types are `wall_lower`, `wall_upper`, and `fibre_fibre`. Wall candidates use `nearest_wall_side` values `lower` or `upper` and do not use fibre-fibre pair semantics. Fibre-fibre candidates use an ordered canonical fibre pair and `none` as the nearest-wall sentinel.

Registry keys, canonical pair keys, canonical sort keys, global candidate IDs, owner ranks for `np=1,2,4`, local IDs, rank counts, and global ordering indices are checked for uniqueness, contiguity, deterministic ordering, and reduction readiness.

## Safety boundary

Stage 21.7 does not compute contact force, collision force, wall force, or contact/collision force application. It does not modify structure advance, RHS coupling, Stage 14 RHS injection, MPI, production DNS, production multifibre logic, Fortran sources, CMake files, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, restart, statistics, or visualization paths.

## Safe defaults

The helper defaults are diagnostic-only and fail-closed: `STAGE21_7_ENABLE=1`, `STAGE21_7_CONTACT_METADATA_CONSISTENCY_ENABLE=1`, `STAGE21_7_DIAGNOSTIC_ONLY=1`, `STAGE21_7_FAIL_CLOSED=1`, `STAGE21_7_NP_VALUES=1,2,4`, `STAGE21_7_METADATA_CHAIN=gap_warning_registry_ownership_ordering`, and all force, structure, RHS, production DNS, MPI, and production multifibre gates disabled.

## Manual command

```bash
stage21_checks/run_stage21_7_contact_metadata_consistency.sh
```

## Expected PASS evidence

The wrapper writes `stage21_outputs/fibre_stage21_7_contact_metadata_consistency.dat` and prints:

```text
STAGE 21.7 CONTACT METADATA CONSISTENCY VERDICT: PASS
STAGE 21.7 FINAL VERDICT: PASS
```

## Next stage

Stage 21.8: contact candidate persistence audit
