# Stage 17.2 channel wall and domain-boundary representation

Stage 17.2 adds a guarded diagnostic-only representation for channel wall and domain-boundary metadata needed by later wall-clearance diagnostics. It records finite ordered channel wall bounds, identifies the wall-normal direction, and makes the x/y/z boundary-policy convention explicit without evaluating fibre positions against walls.

## Boundary metadata

The Stage 17.2 wrapper and helper audit metadata equivalent to:

- `stage17_enable`
- `stage17_wall_safety_enable`
- `stage17_boundary_check_enable`
- `stage17_diagnostic_only`
- `wall_normal_direction = y`
- lower wall identifier `lower_wall_y_min`
- upper wall identifier `upper_wall_y_max`
- `y_min`
- `y_max`
- `channel_height = y_max - y_min`
- `x_boundary_policy`
- `y_boundary_policy`
- `z_boundary_policy`
- `periodic_x_status`
- `periodic_z_status`
- `wall_y_status`
- finite and ordered domain-bound status

Safe defaults keep `diagnostic_only` true, set `wall_normal_direction=y`, use `y_min=-1.0`, `y_max=1.0`, keep x/z policies explicitly periodic, and set the y policy to `wall_bounded`.

## Fail-closed semantics

- `y_min` and `y_max` must be finite.
- `y_max` must be greater than `y_min`.
- `channel_height` must be finite and positive.
- The y boundary policy must be wall-bounded / wall-normal.
- x/z boundary policies must be explicit metadata.
- Invalid wall/domain metadata fails closed with explicit reasons in `stage17_outputs/fibre_stage17_2_channel_wall_domain_boundary.dat`.

## Strict Stage 17.2 boundary

Stage 17.2 must not compute point-wise wall clearance, effective-radius wall gap, or segment wall clearance. It must not classify near-wall/contact state, detect wall penetration, add real wall-contact force, add fibre-fibre collision force, add penalty/repulsive/lubrication/friction/adhesion/contact-damping force, inject collision-induced RHS, update structure from collision, add production multi-fibre collision dynamics, or activate bending/tension/inextensibility.

Point-wise effective-radius wall clearance begins in Stage 17.3. Wall penetration fail-closed checks begin in Stage 17.4. Contact-state classification begins in Stage 17.5. Real contact/collision force models belong to Stage 21. Bending, tension, and inextensibility belong to Stage 18.

## False-positive-safe audit policy

The Stage 17.2 helper reuses the corrected Stage 17.1 / Stage 17.0 / Stage 16 false-positive-safe audit pattern. Documentation is not scanned as executable regression evidence, negative-check strings are not treated as behavior, regex literals such as `rg[[:space:]]` are not treated as ripgrep usage, legitimate Stage 13.5 conservation/sign audits are not treated as production diagnostic regressions, and only targeted Stage 13.6 production/check evidence is inspected for diagnostic naming preservation. Stage 17.0 and Stage 17.1 diagnostic failure-label strings are labels, not rollback evidence.
