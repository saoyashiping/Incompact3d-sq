# Stage 17.6 segment-level wall clearance safety

Stage 17.6 extends the closed Stage 17 point-wise wall diagnostics to adjacent-control-point segment diagnostics.  It remains diagnostic-only: it does not apply wall/contact/fibre-fibre force, mutate fibre state, modify fluid RHS, add production hooks, or activate structure dynamics.

For each segment from `X_i` to `X_{i+1}`, the helper evaluates `y_a`, `y_b`, `y_mid = 0.5 * (y_a + y_b)`, `y_seg_min`, and `y_seg_max`.  It then computes segment centerline wall distances and effective-radius gaps:

```text
seg_centerline_lower = y_seg_min - y_min
seg_centerline_upper = y_max - y_seg_max
seg_centerline_wall = min(seg_centerline_lower, seg_centerline_upper)

seg_gap_lower = y_seg_min - y_min - r_eff
seg_gap_upper = y_max - y_seg_max - r_eff
seg_gap_wall = min(seg_gap_lower, seg_gap_upper)
```

The wrapper supports safe defaults through `STAGE17_6_*` environment variables for channel bounds, effective fibre radius, wall-clearance thresholds, penetration tolerance, point count, test case, and zero tolerance.  It creates `stage17_outputs` and writes `stage17_outputs/fibre_stage17_6_segment_wall_clearance_safety.dat`.

Analytic Stage 17.6 test cases are:

* `all_segments_clear`
* `segment_near_lower_wall`
* `segment_contact_placeholder`
* `segment_lower_penetration`
* `segment_upper_penetration`
* `mixed_segment_states_priority`

Segment-level states may reuse the closed Stage 17.5 diagnostic state names (`CLEAR`, `NEAR_WALL_WARNING`, `CONTACT_PLACEHOLDER`, `PENETRATED_FAIL_CLOSED`) only as force-free diagnostic labels.  Stage 17.6 does not introduce a contact response.

The helper preserves the corrected false-positive-safe audit pattern from Stage 16 and Stage 17.0--17.5: it uses targeted evidence checks, avoids broad repository-wide scanning, does not scan Markdown as executable regression evidence, ignores diagnostic negative-check strings and failure labels, and treats regex literals such as `rg[[:space:]]` as text rather than real ripgrep command use.
