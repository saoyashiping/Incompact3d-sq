# Stage 17.7 contact placeholder interface, no force

Stage 17.7 adds a diagnostic-only placeholder contact-record interface for future wall/contact/collision models.  It consumes the closed Stage 17 wall-gap, containment, contact-state, and segment-safety semantics as read-only evidence and does not modify fibre state, fluid RHS, IBM forcing, production hooks, or DNS-core numerics.

Allowed diagnostic record fields include `contact_pair_type`, `contact_location_type`, entity/point/segment identifiers, `contact_gap`, wall normal components, `contact_state`, placeholder status, and zero-valued force/RHS/structure-update norms.

Allowed diagnostic values are:

* `contact_pair_type`: `WALL_LOWER`, `WALL_UPPER`, `FUTURE_FIBRE_FIBRE_PLACEHOLDER`
* `contact_location_type`: `POINT`, `SEGMENT`, `MOCK_PAIR`
* `contact_state`: `CLEAR`, `NEAR_WALL_WARNING`, `CONTACT_PLACEHOLDER`, `PENETRATED_FAIL_CLOSED`

Wall normals are diagnostic geometry only:

```text
lower wall normal = (0,+1,0)
upper wall normal = (0,-1,0)
```

All contact response channels remain inactive:

```text
contact_force_active_status = false
contact_force_norm = 0
contact_rhs_increment_norm = 0
contact_structure_update_norm = 0
```

Standalone analytic checks cover clear no-record behavior, lower/upper point placeholders, lower segment placeholders, penetrated fail-closed placeholders, and an inactive mock future fibre-fibre placeholder.  The future fibre-fibre placeholder is metadata only and must not activate production multi-fibre logic.

The wrapper writes `stage17_outputs/fibre_stage17_7_contact_placeholder_no_force.dat` and prints the Stage 17.7 PASS/FAIL verdict.  It performs no build, no MPI run, and no production physics.

The helper preserves the corrected false-positive-safe audit pattern from Stage 16 and Stage 17.0--17.6: targeted evidence checks only, no Markdown-as-code scanning, no broad repository-wide scan, no negative-check string regressions, no brittle old failure-label matching, no `rg`-only dependency, and source-only archives without `.git` metadata are not treated as DNS-core contamination.
