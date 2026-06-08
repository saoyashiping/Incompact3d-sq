# Stage 17.8 fibre-fibre collision placeholder as standalone geometry only

Stage 17.8 adds standalone mock two-fibre geometry diagnostics for future fibre-fibre collision placeholder review.  It does not activate production multi-fibre FSI, add production hooks, modify `xcompact3d.f90`, apply fibre-fibre/contact/wall forces, or change fibre/fluid state.

The diagnostic definitions are:

```text
d_ff = minimum centreline distance between selected mock entities
gap_ff = d_ff - 2*r_eff

d_segseg = minimum distance between mock segment A and mock segment B
gap_segseg = d_segseg - 2*r_eff
```

Allowed diagnostic states are:

* `CLEAR`
* `NEAR_FIBRE_WARNING`
* `COLLISION_PLACEHOLDER`
* `OVERLAPPED_FAIL_CLOSED`

The standalone analytic cases cover clear, near-warning, collision-placeholder, overlap fail-closed, segment-segment distance, and mixed-priority behaviour.  The expected priority is `OVERLAPPED_FAIL_CLOSED > COLLISION_PLACEHOLDER > NEAR_FIBRE_WARNING > CLEAR`.

All response channels remain inactive:

```text
fibre_fibre_force_active_status = false
fibre_fibre_force_norm = 0
fibre_fibre_rhs_increment_norm = 0
fibre_fibre_structure_update_norm = 0
```

The wrapper writes `stage17_outputs/fibre_stage17_8_fibre_fibre_placeholder_geometry.dat`, creates `stage17_outputs` if needed, and performs no project build, no MPI run, and no production physics.

The helper preserves the corrected false-positive-safe audit pattern from Stage 16 and Stage 17.0--17.7: targeted evidence checks only, no Markdown-as-code scanning, no broad repository-wide scan, no negative-check string regressions, no brittle old failure-label matching, no `rg`-only dependency, source-only archives without `.git` metadata are not treated as DNS-core contamination, and the Stage 14 small-lambda hook check uses `src/fibre_stage14_production_rhs_injection.f90` plus `src/xcompact3d.f90`.
