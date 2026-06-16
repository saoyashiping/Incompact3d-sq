# Stage 21.1: wall distance and signed-gap audit

Stage 21.1 is a helper-local wall distance and signed-gap audit for one flexible fibre in a channel. It computes wall distances, signed gaps, nearest-wall side, penetration depth, and safe/near-wall/penetration classifications for deterministic audit geometry only.

Stage 21.1 accepts Stage 21.0 PASS evidence when present and preserves Stage 20/Stage 21.0 source-only closure acceptance. Missing old outputs and closure files are allowed when source-only acceptance is enabled, and no previous stage is rerun.

## Wall distance and signed-gap formulas

For channel walls `y_min` and `y_max` and point `X_q = (x_q, y_q, z_q)`:

* Lower-wall distance: `d_lower(q) = y_q - y_min`.
* Upper-wall distance: `d_upper(q) = y_max - y_q`.
* Lower signed gap: `g_lower(q) = d_lower(q) - a_f`.
* Upper signed gap: `g_upper(q) = d_upper(q) - a_f`.
* Pointwise wall gap: `g_wall(q) = min(g_lower(q), g_upper(q))`.
* Global minimum wall gap: `g_wall_min = min_q g_wall(q)`.
* Penetration depth: `penetration_depth(q) = max(0, -g_wall(q))`.
* Maximum penetration depth: `penetration_depth_max = max_q penetration_depth(q)`.
* Nearest wall side: lower when `g_lower(q) <= g_upper(q)`, otherwise upper.

## Non-penetration and classification

The non-penetration requirement is `g_lower(q) >= 0`, `g_upper(q) >= 0`, and `g_wall_min >= 0`.

Classifications are helper-local only:

* `safe` when `g_wall(q) > warning_gap`.
* `near_wall` when `0 <= g_wall(q) <= warning_gap`.
* `penetration` when `g_wall(q) < 0`.
* `fail` for NaN/Inf or `penetration_depth_max > penetration_fail_limit`.

## Safe default gates

Wall contact, wall contact force candidate, contact force application, contact in structure advance, contact to RHS, fibre collision, fibre-fibre gap audit, fibre collision force candidate, production multifibre, production DNS, and actual MPI are disabled. Diagnostic-only and fail-closed modes are enabled.

## No production activation

Stage 21.1 computes wall distances and wall signed gaps only. It does not compute wall contact force, fibre-fibre distance, collision force, penalty/repulsive/lubrication/friction/adhesion/contact damping force, contact-augmented structure advance, contact-to-RHS forcing, Stage 14 RHS injection, MPI, production DNS, IBM/DNS-core/projection/Poisson/RK3/channel-forcing changes, or restart/statistics/visualization production I/O changes.

## Output

```text
stage21_outputs/fibre_stage21_1_wall_distance_signed_gap_audit.dat
```

## Manual command

```bash
stage21_checks/run_stage21_1_wall_distance_signed_gap_audit.sh
```

## Expected PASS evidence

```text
STAGE 21.1 WALL DISTANCE SIGNED GAP AUDIT VERDICT: PASS
STAGE 21.1 FINAL VERDICT: PASS
```

## Next stage

Stage 21.2: fibre-fibre point/segment distance audit
