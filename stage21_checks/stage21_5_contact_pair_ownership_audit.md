# Stage 21.5: contact pair ownership audit

Stage 21.5 is a diagnostic-only ownership audit for wall-contact and fibre-fibre contact candidate records. It verifies deterministic, unique, and reduction-ready ownership metadata under helper-level `np=1`, `np=2`, and `np=4` semantics without running MPI or production DNS.

Stage 21.5 accepts Stage 21.4 PASS evidence when present and preserves source-only archive acceptance. Missing old outputs are allowed when accepted evidence or source-only mode is available, and no previous stage is rerun.

## Ownership keys

Wall candidate canonical key:

```text
(candidate_type, fibre_i, point_i, nearest_wall_side)
```

Fibre-fibre candidate canonical key:

```text
(candidate_type, min(fibre_i,fibre_j), max(fibre_i,fibre_j), min(segment_i,segment_j), max(segment_i,segment_j))
```

The documented ownership rule is `owner_rank = stable_hash(candidate_key) mod np`, where `stable_hash` is a deterministic SHA-256 integer hash of the canonical key string.

## Required helper semantics

The audit evaluates `np_values = 1,2,4`, assigns owner ranks, assigns contiguous local candidate IDs per owner rank, records rank candidate counts, checks unique global IDs and candidate IDs, rejects duplicate/self/unordered pair contamination, and verifies reduction-ready ownership metadata.

## No production activation

Stage 21.5 does not compute contact force, collision force, wall force, contact/collision force application, structure advance updates, RHS updates, Stage 14 RHS injection, MPI, production DNS, production multifibre logic, source changes, CMake changes, or production hooks.

## Output

```text
stage21_outputs/fibre_stage21_5_contact_pair_ownership_audit.dat
```

## Manual command

```bash
stage21_checks/run_stage21_5_contact_pair_ownership_audit.sh
```

## Expected PASS evidence

```text
STAGE 21.5 CONTACT PAIR OWNERSHIP AUDIT VERDICT: PASS
STAGE 21.5 FINAL VERDICT: PASS
```

## Next stage

Stage 21.6: deterministic pair ordering
