# Stage 17.1 wall/contact safety configuration

Stage 17.1 adds the guarded configuration layer for Stage 17 wall / boundary / contact-safety diagnostics. It is diagnostic-only and non-invasive: it does not compute wall distance, segment wall clearance, near-wall state, wall penetration, or contact-state classification, and it does not add real contact or collision physics.

## Configuration controls

The Stage 17.1 wrapper and helper audit controls equivalent to:

- `stage17_enable`
- `stage17_wall_safety_enable`
- `stage17_boundary_check_enable`
- `stage17_fail_closed_enable`
- `stage17_contact_placeholder_enable`
- `stage17_fibre_collision_placeholder_enable`
- `stage17_effective_fibre_radius`
- `stage17_min_wall_clearance`
- `stage17_warning_wall_clearance`
- `stage17_penetration_tolerance`
- `stage17_diagnostic_only`

The safe defaults are diagnostic-only and non-invasive. They leave real wall-contact physics, fibre-fibre collision physics, contact/collision force, `X_f`, `V_f`, `A_f`, RHS, Stage 13 force-density, Stage 14 RHS, Stage 16 closed-loop behavior, pressure/projection/Poisson/RK3 logic, and channel forcing untouched.

## Fail-closed numeric semantics

- `stage17_effective_fibre_radius` must be finite and nonnegative.
- `stage17_min_wall_clearance` must be finite and nonnegative.
- `stage17_warning_wall_clearance` must be finite and greater than or equal to `stage17_min_wall_clearance`.
- `stage17_penetration_tolerance` must be finite and nonnegative.
- `stage17_diagnostic_only` must remain true for Stage 17.1.

Invalid radius, clearance, or tolerance values fail closed with explicit reasons in `stage17_outputs/fibre_stage17_1_wall_contact_safety_config.dat`.

## Boundary for later stages

Wall distance diagnostics begin later in Stage 17.3. Contact-state classification begins later in Stage 17.5. Real wall-contact and fibre-fibre collision force models belong to Stage 21. Bending, tension, and inextensibility structure dynamics belong to Stage 18.

## False-positive-safe audit policy

The Stage 17.1 helper reuses the corrected Stage 17.0 / Stage 16 false-positive-safe audit pattern. Documentation is not scanned as executable regression evidence, negative-check strings are not treated as behavior, regex literals such as `rg[[:space:]]` are not treated as ripgrep usage, legitimate Stage 13.5 conservation/sign audits are not treated as production diagnostic regressions, and only targeted Stage 13.6 production/check evidence is inspected for diagnostic naming preservation.
