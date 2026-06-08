# Stage 17.3 effective-radius wall-clearance diagnostics

Stage 17.3 adds guarded diagnostic-only point-wise wall-clearance diagnostics for standalone analytic fibre control-point geometries. It uses the Stage 17.1 effective fibre radius concept and the Stage 17.2 channel wall/domain-boundary representation without modifying production fibre state or DNS-core numerics.

## Diagnostic formulas

For each control point with coordinate `y_i`:

- `d_lower(i) = y_i - y_min`
- `d_upper(i) = y_max - y_i`
- `d_center_wall(i) = min(d_lower(i), d_upper(i))`
- `gap_lower(i) = d_lower(i) - r_eff`
- `gap_upper(i) = d_upper(i) - r_eff`
- `gap_wall(i) = min(gap_lower(i), gap_upper(i))`
- `d_center_wall_min = min_i d_center_wall(i)`
- `gap_wall_min = min_i gap_wall(i)`

All computed distances and gaps must be finite. Numeric values are reported as diagnostics and are not treated as boolean pass/fail fields in `final_status`.

## Analytic test geometries

The helper supports standalone analytic cases:

- `centered_clear`: points safely away from both walls with positive centerline distance and positive effective-radius gap.
- `near_lower_wall_gap_positive`: one point close to the lower wall but with positive effective-radius gap.
- `exact_radius_touch_placeholder`: one point with lower-wall distance equal to `r_eff`, producing zero effective-radius gap within tolerance.
- `outside_gap_negative_diagnostic_only`: one point with lower-wall distance less than `r_eff`, producing a negative gap that is reported diagnostically only.

## Strict Stage 17.3 boundary

Stage 17.3 does not classify contact state, issue near-wall warning states, perform wall-penetration fail-closed checks, apply wall/contact/collision forces, modify `X_f`, `V_f`, or `A_f`, modify fluid RHS, insert production hooks into `xcompact3d.f90`, activate production multi-fibre logic, or activate bending/tension/inextensibility. Near-wall/contact-state classification begins in Stage 17.5, wall-penetration fail-closed checks begin in Stage 17.4, real contact/collision force models belong to Stage 21, and bending/tension/inextensibility belong to Stage 18.

## False-positive-safe audit policy

The Stage 17.3 helper reuses the corrected Stage 17.2 / Stage 17.1 / Stage 17.0 / Stage 16 false-positive-safe audit pattern. Documentation is not scanned as executable regression evidence, negative-check strings are not treated as behavior, regex literals such as `rg[[:space:]]` are not treated as ripgrep usage, legitimate Stage 13.5 conservation/sign audits are not treated as production diagnostic regressions, and only targeted Stage 13.6 production/check evidence is inspected for diagnostic naming preservation. Stage 17.0, Stage 17.1, and Stage 17.2 diagnostic failure-label strings are labels, not rollback evidence.
