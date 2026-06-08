# Stage 17.4 boundary containment and wall-penetration fail-closed

Stage 17.4 adds guarded diagnostic-only boundary-containment checks for standalone analytic fibre control-point geometries. It uses the Stage 17.1 effective fibre radius concept, Stage 17.2 channel wall/domain-boundary metadata, and Stage 17.3 effective-radius wall-clearance gap formulas without applying forces or modifying production fibre/fluid state.

## Penetration formulas

For each control point with coordinate `y_i`:

- `lower_surface_gap(i) = y_i - r_eff - y_min`
- `upper_surface_gap(i) = y_max - (y_i + r_eff)`

These are equivalent to the Stage 17.3 effective-radius gaps:

- `gap_lower(i) = y_i - y_min - r_eff`
- `gap_upper(i) = y_max - y_i - r_eff`

Lower-wall penetration exists when `gap_lower(i) < -penetration_tolerance`. Upper-wall penetration exists when `gap_upper(i) < -penetration_tolerance`. A point is boundary-contained when both gaps are greater than or equal to `-penetration_tolerance`.

## Analytic test geometries

The helper supports standalone analytic cases:

- `contained_clear`: all points are safely inside the wall-bounded channel.
- `exact_radius_touch_lower`: one point has `y_i - r_eff = y_min` within tolerance and is non-penetrating.
- `exact_radius_touch_upper`: one point has `y_i + r_eff = y_max` within tolerance and is non-penetrating.
- `lower_wall_penetration`: one point has lower-wall penetration and must fail closed with an explicit lower-wall reason when selected.
- `upper_wall_penetration`: one point has upper-wall penetration and must fail closed with an explicit upper-wall reason when selected.

## Strict Stage 17.4 boundary

Stage 17.4 reports containment, penetration detection, penetration depths, and offending point indices. It does not classify contact states such as CLEAR / NEAR_WALL / CONTACT_PLACEHOLDER, does not issue near-wall warning states, does not apply wall/contact/collision forces, does not modify `X_f`, `V_f`, or `A_f`, does not modify fluid RHS, does not insert production hooks into `xcompact3d.f90`, does not activate production multi-fibre logic, and does not activate bending/tension/inextensibility.

Contact-state and near-wall classification begins in Stage 17.5. Segment-level wall safety begins in Stage 17.6. Real contact/collision force models belong to Stage 21. Bending, tension, and inextensibility belong to Stage 18.

## False-positive-safe audit policy

The Stage 17.4 helper reuses the corrected Stage 17.3 / Stage 17.2 / Stage 17.1 / Stage 17.0 / Stage 16 false-positive-safe audit pattern. Documentation is not scanned as executable regression evidence, negative-check strings are not treated as behavior, regex literals such as `rg[[:space:]]` are not treated as ripgrep usage, legitimate Stage 13.5 conservation/sign audits are not treated as production diagnostic regressions, and only targeted Stage 13.6 production/check evidence is inspected for diagnostic naming preservation. Stage 17.0, Stage 17.1, Stage 17.2, and Stage 17.3 diagnostic failure-label strings are labels, not rollback evidence.
