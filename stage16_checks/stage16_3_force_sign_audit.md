# Stage 16.3 structure-side force sign / action-reaction audit

Stage 16.3 is a sign-convention audit before any closed-loop one-fibre FSI
production update is activated. It verifies the Stage 12 feedback-force sign
used by the future structure-side input and the equal/opposite fluid-side sign
used by the Stage 13/14 force-density/RHS path.

Stage 16.3 does **not** advance the structure, insert a production hook, modify
fluid RHS, modify pressure/projection/Poisson/RK3/channel forcing, solve
bending/tension, handle wall/contact logic, or activate multi-fibre logic.

## Sign convention under audit

The controlled analytic relation is:

```text
slip = U_f - V_f
F_cand = alpha * slip
F_fluid_on_fibre = F_cand
F_fibre_on_fluid = -F_cand
F_fluid_on_fibre + F_fibre_on_fluid = 0
```

The future Stage 16 structure update must use `F_fluid_on_fibre`. The Stage
13/14 fluid-side force-density/RHS path must use `F_fibre_on_fluid`.

## Manual command

Run exactly:

```bash
bash stage16_checks/run_stage16_3_force_sign_audit.sh
```

The wrapper builds and runs only the standalone Stage 16.3 force-sign audit
check target, then writes:

```text
stage16_outputs/fibre_stage16_3_force_sign_audit.dat
```

A successful run prints:

```text
STAGE 16.3 FORCE SIGN AUDIT VERDICT: PASS
STAGE 16.3 FINAL VERDICT: PASS
```

A failed run prints `FAIL` and explicit failure reasons.

## Environment variables

| Variable | Default | Purpose |
| --- | ---: | --- |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `DECOMP2D_ROOT` | empty | Optional CMake prefix path when configuring a missing build directory. |
| `MPIEXEC` | `mpirun` | Preserved for wrapper consistency. |
| `MPIEXEC_FLAGS` | empty | Preserved for wrapper consistency. |
| `STAGE16_3_RUN_STAGE16_2` | `0` | Optional prerequisite Stage 16.2 rerun. |
| `STAGE16_3_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE16_3_REQUIRE_STAGE15_CLOSED` | `1` | Require `stage15_checks/STAGE15_CLOSED.md`. |
| `STAGE16_3_REQUIRE_STAGE16_2` | `1` | Require Stage 16.2 evidence. |
| `STAGE16_3_ACCEPT_STAGE16_2_CLOSED_EVIDENCE` | `1` | Accept existing Stage 16.2 files when fresh output evidence is absent. |
| `STAGE16_3_ENABLE` | `1` | Stage 16.3 sign-audit request gate. |
| `STAGE16_3_ONE_FIBRE_FSI_ENABLE` | `1` | One-fibre sign-audit request only; not production FSI activation. |
| `STAGE16_3_FEEDBACK_ALPHA` | `1.0` | Positive analytic feedback gain. |
| `STAGE16_3_MAX_ACTION_REACTION_ERROR` | `1.0e-14` | Action-reaction tolerance. |
| `STAGE16_3_MAX_SIGN_ERROR` | `1.0e-14` | Sign-check tolerance. |
| `STAGE16_3_DIAGNOSTIC_ONLY` | `1` | Stage 16.3 must remain diagnostic-only. |

## Added force-sign API

`src/fibre_stage16_force_sign_audit.f90` provides:

- `stage16_force_sign_audit_reset()`
- `stage16_force_sign_audit_load_from_environment()`
- `stage16_force_sign_audit_compute_reference(u_f, v_f, alpha)`
- `stage16_force_sign_audit_validate_action_reaction(tolerance)`
- `stage16_force_sign_audit_validate_wrong_sign_rejection(tolerance)`
- `stage16_force_sign_audit_write_diagnostics(unit)`

The API is standalone and is not connected to the production time loop.

## Audit coverage

Stage 16.3 verifies:

1. Stage 16.3 wrapper/doc/source/helper/check target and minimal build
   registration exist.
2. Stage 14 and Stage 15 closure files exist when required.
3. Stage 16.2 evidence is present or accepted via closed Stage 16.2 files.
4. The analytic slip and candidate force are finite.
5. The structure-side force is explicitly fluid-on-fibre.
6. The fluid-side force is explicitly fibre-on-fluid.
7. The two signs satisfy action-reaction within tolerance.
8. Wrong-sign fluid-side usage is detected and rejected.
9. Zero slip gives zero force within tolerance.
10. Reversing slip reverses both force signs consistently.
11. Stage 16.3 does not insert production hooks into `src/xcompact3d.f90`.
12. Stage 16.3 does not activate structure advance, RHS modification,
    bending/tension solves, wall/contact handling, multi-fibre logic,
    pressure/projection/Poisson/RK3/channel-forcing changes, or legacy IBM
    forcing outside the approved chain.
13. Stage 11.5, Stage 13.6, Stage 14, Stage 15.1-15.11, Stage 16.0, Stage
    16.1, and Stage 16.2 diagnostics remain present.
14. Rank0-safe diagnostic writing is preserved for Stage 11, Stage 13,
    Stage 14, Stage 15, and Stage 16 diagnostics.
15. No Stage 16.3 script requires ripgrep without a grep fallback.

## Required summary fields

`stage16_outputs/fibre_stage16_3_force_sign_audit.dat` contains at least:

```text
stage16_3_requested_status
feedback_alpha
slip_finite_status
feedback_force_finite_status
structure_side_force_sign_status
fluid_side_force_sign_status
action_reaction_error
action_reaction_status
wrong_sign_rejection_status
zero_slip_zero_force_status
slip_reversal_sign_status
approved_stage12_13_14_chain_status
no_production_hook_status
no_structure_advance_status
no_rhs_modification_status
no_bending_solve_status
no_tension_solve_status
no_wall_contact_status
no_multifibre_status
no_pressure_projection_modification_status
no_poisson_modification_status
no_rk3_channel_forcing_modification_status
stage14_regression_status
stage15_regression_status
stage16_1_regression_status
stage16_2_regression_status
final_status
```

`final_status 1` is required for PASS.

## Stage 16.3 boundary

Stage 16.3 intentionally adds only force-sign audit source/check code, a wrapper,
a parser/audit helper, documentation, and minimal standalone build registration.
It does not add one-fibre FSI physics, compute production feedback, advance the
structure, solve bending/tension, handle wall/contact interactions, enable
multi-fibre operation, modify fluid RHS, modify pressure/projection/Poisson/RK3
or channel forcing, or activate legacy IBM forcing outside the approved Stage
11-14 chain.

## 2026-06-05 audit-helper correction

Stage 16.3 static audit must not treat documentation or negative audit text as a real regression.  In particular, references to the forbidden Stage 14 `stage14_get_injection_gain()==0.0` gate inside Markdown or explanatory text are not production code, and Stage 13.5 remains a legitimate closed conservation/sign audit.  Only the old production force-density candidate naming pattern is forbidden.  The helper also distinguishes real shell `rg` command usage from quoted regular-expression text so the no-rg-only-dependency check does not produce false positives.
