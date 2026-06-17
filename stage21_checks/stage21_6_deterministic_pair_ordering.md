# Stage 21.6: deterministic pair ordering

Stage 21.6 is a diagnostic-only deterministic ordering audit for wall-contact and fibre-fibre contact candidate records. It verifies that canonical global ordering is independent of discovery order, helper owner rank, helper `np=1/2/4` grouping, and repeated evaluations.

Stage 21.6 accepts Stage 21.5 PASS evidence when present and preserves source-only archive acceptance. Missing old outputs are allowed when accepted evidence or source-only mode is available, and no previous stage is rerun.

## Canonical sort key

Wall candidate canonical sort key:

```text
(candidate_type_priority, fibre_i, point_i, nearest_wall_side_priority, candidate_id)
```

Fibre-fibre candidate canonical sort key:

```text
(candidate_type_priority, min(fibre_i,fibre_j), max(fibre_i,fibre_j), min(segment_i,segment_j), max(segment_i,segment_j), min(point_i,point_j), max(point_i,point_j), candidate_id)
```

Candidate type priority is `wall_lower=0`, `wall_upper=1`, and `fibre_fibre=2`. Nearest-wall side priority is `lower=0` and `upper=1`.

## Ordering invariance checks

The audit verifies that sorting by `canonical_sort_key` produces the same reference ordering from original discovery order, reversed order, fixed-seed shuffled order, owner-rank grouped order for `np=1`, `np=2`, and `np=4`, and repeated helper evaluations.

## No production activation

Stage 21.6 does not compute contact force, collision force, wall force, contact/collision force application, structure advance updates, RHS updates, Stage 14 RHS injection, MPI, production DNS, production multifibre logic, source changes, CMake changes, or production hooks.

## Output

```text
stage21_outputs/fibre_stage21_6_deterministic_pair_ordering.dat
```

## Manual command

```bash
stage21_checks/run_stage21_6_deterministic_pair_ordering.sh
```

## Expected PASS evidence

```text
STAGE 21.6 DETERMINISTIC PAIR ORDERING VERDICT: PASS
STAGE 21.6 FINAL VERDICT: PASS
```

## Next stage

Stage 21.7: contact metadata consistency
