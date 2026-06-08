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

## Stage 17.6 helper evidence-audit correction

The Stage 17.6 helper must not regress to brittle checks for previously closed stages.  In fresh source archives, `.git` metadata may be absent and generated runtime outputs may be missing; this must not be treated as closed-file modification.  Stage 17.6 therefore accepts read-only structural evidence for closed Stage 17.0--17.5 files and preserves the earlier fixes:

* Stage 17.0 fresh-archive Stage 16 closure acceptance is recognized by its helper/wrapper/documentation structure rather than by old failure-label strings.
* Stage 17.1 evidence/final-status handling is recognized whether the accepted helper uses `VALUE_KEYS` or the earlier `pass_fail_keys` pattern to exclude numeric value fields.
* Git-status unavailability in source-only archives is not interpreted as RHS/IBM/DNS-core contamination.
* Stage 14 small-nonzero-lambda hook evidence is checked against `src/fibre_stage14_production_rhs_injection.f90` and `src/xcompact3d.f90`; the obsolete/nonexistent `fibre_stage14_rhs_apply.f90` name is not used.
