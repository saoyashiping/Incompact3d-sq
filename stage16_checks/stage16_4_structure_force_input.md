# Stage 16.4 structure-side feedback force input hook

Stage 16.4 connects a controlled Stage 12-style feedback-force candidate into a
Stage 16 structure-side force input buffer using the Stage 16.3-audited
structure-side sign. It is an interface/buffer audit only: it does not advance
the structure, update position/velocity/acceleration, insert a production hook,
modify fluid RHS, or alter pressure/projection/Poisson/RK3/channel forcing.

## Sign convention

The controlled relation remains:

```text
F_cand = alpha * (U_f - V_f)
F_structure_input = F_fluid_on_fibre = F_cand
F_fibre_on_fluid = -F_structure_input
```

The Stage 13/14 fluid-side path continues to use the opposite fibre-on-fluid
sign. Stage 16.4 stores only the read-only structure-side force input for future
Stage 16 structure-update stages.

## Manual command

Run exactly:

```bash
bash stage16_checks/run_stage16_4_structure_force_input.sh
```

The wrapper builds and runs only the standalone Stage 16.4 structure force-input
check target, then writes:

```text
stage16_outputs/fibre_stage16_4_structure_force_input.dat
```

A successful run prints:

```text
STAGE 16.4 STRUCTURE FORCE INPUT VERDICT: PASS
STAGE 16.4 FINAL VERDICT: PASS
```

A failed run prints `FAIL` and explicit failure reasons.

## Environment variables

| Variable | Default | Purpose |
| --- | ---: | --- |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `DECOMP2D_ROOT` | empty | Optional CMake prefix path when configuring a missing build directory. |
| `MPIEXEC` | `mpirun` | Preserved for wrapper consistency. |
| `MPIEXEC_FLAGS` | empty | Preserved for wrapper consistency. |
| `STAGE16_4_RUN_STAGE16_3` | `0` | Optional prerequisite Stage 16.3 rerun. |
| `STAGE16_4_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE16_4_REQUIRE_STAGE15_CLOSED` | `1` | Require `stage15_checks/STAGE15_CLOSED.md`. |
| `STAGE16_4_REQUIRE_STAGE16_3` | `1` | Require Stage 16.3 evidence. |
| `STAGE16_4_ACCEPT_STAGE16_3_CLOSED_EVIDENCE` | `1` | Accept existing Stage 16.3 files when fresh output evidence is absent. |
| `STAGE16_4_ENABLE` | `1` | Stage 16.4 force-input audit request gate. |
| `STAGE16_4_ONE_FIBRE_FSI_ENABLE` | `1` | One-fibre force-input audit request only; not production FSI activation. |
| `STAGE16_4_FEEDBACK_ALPHA` | `1.0` | Positive controlled feedback gain. |
| `STAGE16_4_MAX_ACTION_REACTION_ERROR` | `1.0e-14` | Action-reaction tolerance. |
| `STAGE16_4_MAX_SIGN_ERROR` | `1.0e-14` | Sign/readback tolerance. |
| `STAGE16_4_MAX_FORCE_INPUT` | `1.0e-6` | Maximum allowed structure-side input magnitude. |
| `STAGE16_4_DIAGNOSTIC_ONLY` | `1` | Stage 16.4 must remain diagnostic-only. |

## Added force-input API

`src/fibre_stage16_structure_force_input.f90` provides:

- `stage16_structure_force_input_reset()`
- `stage16_structure_force_input_load_from_environment()`
- `stage16_structure_force_input_set_from_stage12_candidate(u_f, v_f, alpha, use_wrong_sign)`
- `stage16_structure_force_input_validate()`
- `stage16_structure_force_input_get_force(index, force)`
- `stage16_structure_force_input_write_diagnostics(unit)`

The API is standalone and exposes read-only force-input diagnostics for future
Stage 16.5/16.6 work. It is not connected to the production time loop.

## Audit coverage

Stage 16.4 verifies:

1. Stage 16.4 wrapper/doc/source/helper/check target and minimal build
   registration exist.
2. Stage 14 and Stage 15 closure files exist when required.
3. Stage 16.3 evidence is present or accepted via closed Stage 16.3 files.
4. The controlled Stage 12-style candidate is finite.
5. The stored structure input uses the fluid-on-fibre sign.
6. The fluid-side companion uses the equal/opposite fibre-on-fluid sign.
7. Action-reaction consistency is within tolerance.
8. Wrong-sign structure-side input is rejected.
9. Zero slip gives zero structure input within tolerance.
10. Structure force input is bounded by `STAGE16_4_MAX_FORCE_INPUT`.
11. Force input reset/clear and readback behavior are valid.
12. Stage 16.4 does not insert production hooks into `src/xcompact3d.f90`.
13. Stage 16.4 does not activate structure advance, position/velocity/
    acceleration updates, RHS modification, bending/tension solves, wall/contact
    handling, multi-fibre logic, pressure/projection/Poisson/RK3/channel-forcing
    changes, or legacy IBM forcing outside the approved chain.
14. Stage 11.5, Stage 13.6, Stage 14, Stage 15.1-15.11, and Stage 16.0-16.3
    diagnostics remain present.
15. Rank0-safe diagnostic writing is preserved for Stage 11, Stage 13, Stage 14,
    Stage 15, and Stage 16 diagnostics.
16. Static checks avoid documentation/negative-check/regex-literal false
    positives and do not treat the legitimate Stage 13.5 conservation/sign audit
    as an old production diagnostic regression.

## Required summary fields

`stage16_outputs/fibre_stage16_4_structure_force_input.dat` contains at least:

```text
stage16_4_requested_status
feedback_alpha
structure_force_input_finite_status
structure_force_input_bounded_status
structure_side_force_sign_status
fluid_side_force_sign_status
action_reaction_error
action_reaction_status
wrong_sign_rejection_status
zero_slip_zero_input_status
force_input_reset_status
force_input_readback_status
approved_stage12_13_14_chain_status
no_production_hook_status
no_structure_advance_status
no_position_update_status
no_velocity_update_status
no_acceleration_update_status
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
stage16_3_regression_status
final_status
```

`final_status 1` is required for PASS.

## Stage 16.4 boundary

Stage 16.4 intentionally adds only force-input interface source/check code, a
wrapper, a parser/audit helper, documentation, and minimal standalone build
registration. It does not add one-fibre FSI production physics, advance the
structure, update position/velocity/acceleration, solve bending/tension, handle
wall/contact interactions, enable multi-fibre operation, modify fluid RHS,
modify pressure/projection/Poisson/RK3 or channel forcing, or activate legacy IBM
forcing outside the approved Stage 11-14 chain.
