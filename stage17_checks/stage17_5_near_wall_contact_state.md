# Stage 17.5 near-wall warning and contact-state classification

Stage 17.5 adds guarded diagnostic-only wall-safety state classification for standalone analytic fibre control-point geometries. It uses Stage 17.1 effective-radius configuration, Stage 17.2 channel wall/domain-boundary metadata, Stage 17.3 wall-gap diagnostics, and Stage 17.4 penetration fail-closed meaning without applying forces or modifying production fibre/fluid state.

## Diagnostic states

Allowed point-wise states are:

- `CLEAR`
- `NEAR_WALL_WARNING`
- `CONTACT_PLACEHOLDER`
- `PENETRATED_FAIL_CLOSED`

For a point-wise effective-radius gap `gap_wall(i) = min(gap_lower(i), gap_upper(i))`, Stage 17.5 classifies:

1. `PENETRATED_FAIL_CLOSED` when `gap_wall(i) < -penetration_tolerance`.
2. `CONTACT_PLACEHOLDER` when `-penetration_tolerance <= gap_wall(i) <= min_wall_clearance`.
3. `NEAR_WALL_WARNING` when `min_wall_clearance < gap_wall(i) <= warning_wall_clearance`.
4. `CLEAR` when `gap_wall(i) > warning_wall_clearance`.

The global worst-state priority is `PENETRATED_FAIL_CLOSED > CONTACT_PLACEHOLDER > NEAR_WALL_WARNING > CLEAR`.

## Analytic test geometries

The helper supports standalone analytic cases:

- `all_clear`: all points classify as `CLEAR`.
- `near_wall_warning`: one point classifies as `NEAR_WALL_WARNING`.
- `contact_placeholder`: one point classifies as `CONTACT_PLACEHOLDER` with no force.
- `penetrated_fail_closed`: one point classifies as `PENETRATED_FAIL_CLOSED` and reports its index.
- `mixed_states_priority`: points include all four states and verify global priority.

## Strict Stage 17.5 boundary

Stage 17.5 classifies diagnostic states, reports counts, reports the global worst state, and reports penetrated point indices. It does not apply wall/contact/fibre-fibre forces, penalty/repulsive/lubrication/friction/adhesion/contact-damping forces, modify `X_f`, `V_f`, or `A_f`, modify fluid RHS, insert production hooks into `xcompact3d.f90`, activate production multi-fibre logic, or activate bending/tension/inextensibility.

Segment-level wall safety begins in Stage 17.6. The contact placeholder interface with no force begins in Stage 17.7. Real contact/collision force models belong to Stage 21. Bending, tension, and inextensibility belong to Stage 18.

## False-positive-safe audit policy

The Stage 17.5 helper reuses the corrected Stage 17.4 / Stage 17.3 / Stage 17.2 / Stage 17.1 / Stage 17.0 / Stage 16 false-positive-safe audit pattern. Documentation is not scanned as executable regression evidence, negative-check strings are not treated as behavior, regex literals such as `rg[[:space:]]` are not treated as ripgrep usage, legitimate Stage 13.5 conservation/sign audits are not treated as production diagnostic regressions, and only targeted Stage 13.6 production/check evidence is inspected for diagnostic naming preservation. Stage 17.0--17.4 diagnostic failure-label strings are labels, not rollback evidence.
