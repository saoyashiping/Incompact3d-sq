# Stage 21.2: fibre-fibre point/segment distance audit

Stage 21.2 is a helper-local fibre-fibre point/segment distance and signed-gap audit for controlled two-fibre geometry. It computes point-point distances, segment-segment distances, closest point/segment IDs, fibre-fibre signed gaps, finite/bounded distance metadata, ownership metadata, no self-pair/duplicate-pair contamination, and safe/near-contact/overlap classification.

Stage 21.2 accepts Stage 21.1 PASS evidence when present and preserves Stage 20/Stage 21.0/Stage 21.1 source-only closure acceptance. Missing old outputs and closure files are allowed when source-only acceptance is enabled, and no previous stage is rerun.

## Point-point distance formulas

For two fibre centerlines `X_i(s)` and `X_j(r)`, the point-point distance is `d_pq = |X_i,p - X_j,q|`, the closest point-point distance is `d_pp_min = min_{p,q} d_pq`, and the point-point signed gap is `g_pp = d_pp_min - (a_i + a_j)`.

## Segment-segment distance formulas

Each segment pair uses endpoints `A = X_i,p`, `B = X_i,p+1`, `C = X_j,q`, and `D = X_j,q+1`. Closest points are `P* = A + s_clamp * (B - A)` and `Q* = C + t_clamp * (D - C)` with `0 <= s_clamp <= 1` and `0 <= t_clamp <= 1`. The segment distance is `d_seg = |P* - Q*|`, the closest segment distance is `d_ss_min = min_{p,q} d_seg`, and the segment signed gap is `g_ss = d_ss_min - (a_i + a_j)`.

## Fibre-fibre signed gap and classification

The final fibre-fibre distance and signed gap are `d_ff = min(d_pp_min, d_ss_min)` and `g_ff = d_ff - (a_i + a_j)`. The non-penetration requirement is `g_ff >= 0`; penetration depth is `penetration_depth_ff = max(0, -g_ff)`.

Classifications are helper-local only: `safe` when `g_ff > warning_gap`, `near_contact` when `0 <= g_ff <= warning_gap`, `overlap_or_penetration` when `g_ff < 0`, and `fail` for NaN/Inf or `penetration_depth_ff > penetration_fail_limit`.

## Safe default gates

Fibre collision, fibre collision force candidate, wall contact, wall contact force candidate, contact/collision force application, contact in structure advance, contact to RHS, production multifibre, production DNS, and actual MPI are disabled. Diagnostic-only and fail-closed modes are enabled.

## No production activation

Stage 21.2 computes fibre-fibre distances and signed gaps only. It does not compute collision force, wall contact force, penalty/repulsive/lubrication/friction/adhesion/contact damping force, contact/collision force application, contact/collision structure advance, contact/collision RHS forcing, Stage 14 RHS injection, MPI, production DNS, IBM/DNS-core/projection/Poisson/RK3/channel-forcing changes, or restart/statistics/visualization production I/O changes.

The two-fibre geometry is helper-local diagnostic geometry only and is not production multifibre logic.

## Output

```text
stage21_outputs/fibre_stage21_2_fibre_fibre_distance_audit.dat
```

## Manual command

```bash
stage21_checks/run_stage21_2_fibre_fibre_distance_audit.sh
```

## Expected PASS evidence

```text
STAGE 21.2 FIBRE-FIBRE DISTANCE AUDIT VERDICT: PASS
STAGE 21.2 FINAL VERDICT: PASS
```

## Next stage

Stage 21.3: near-contact warning and fail-closed gate
