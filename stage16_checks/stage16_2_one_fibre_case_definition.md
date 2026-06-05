# Stage 16.2 single-fibre production case definition

Stage 16.2 defines a safe, standalone one-fibre channel-DNS FSI case description
for future Stage 16 closed-loop tests. It does **not** compute feedback force,
advance the structure, insert a production hook, modify RHS forcing, or alter any
pressure/projection/Poisson/RK3/channel-forcing behavior.

The future coupling route remains:

```text
Stage 15 controlled structure-state update
-> Stage 12 feedback-force candidate
-> Stage 13 Eulerian force-density candidate
-> Stage 14 controlled RHS diagnostic/injection path
```

Stage 16.2 only defines and validates initial placeholders:

```text
X_f(s,0), V_f(s,0), A_f(s,0)
```

## Manual command

Run exactly:

```bash
bash stage16_checks/run_stage16_2_one_fibre_case_definition.sh
```

The wrapper builds and runs only the standalone Stage 16.2 case-definition check
target, then writes:

```text
stage16_outputs/fibre_stage16_2_one_fibre_case_definition.dat
```

A successful run prints:

```text
STAGE 16.2 ONE-FIBRE CASE DEFINITION VERDICT: PASS
STAGE 16.2 FINAL VERDICT: PASS
```

A failed run prints `FAIL` and explicit failure reasons.

## Safe default case

The standalone default case uses a normalized safe box because Stage 16.2 is a
case-definition/check stage and does not claim production-grid insertion yet:

- exactly one fibre,
- `N_f = 8`,
- deterministic straight fibre from `(0.25, 0.50, 0.50)` to `(0.75, 0.50, 0.50)`,
- zero initial velocity placeholder `V_f = 0`,
- zero initial acceleration placeholder `A_f = 0`,
- positive point spacing,
- finite wall/domain clearance,
- wall/contact and multi-fibre logic inactive except for rejection diagnostics.

## Environment variables

| Variable | Default | Purpose |
| --- | ---: | --- |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `DECOMP2D_ROOT` | empty | Optional CMake prefix path when configuring a missing build directory. |
| `MPIEXEC` | `mpirun` | Preserved for wrapper consistency. |
| `MPIEXEC_FLAGS` | empty | Preserved for wrapper consistency. |
| `STAGE16_2_RUN_STAGE16_1` | `0` | Optional prerequisite Stage 16.1 rerun. |
| `STAGE16_2_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE16_2_REQUIRE_STAGE15_CLOSED` | `1` | Require `stage15_checks/STAGE15_CLOSED.md`. |
| `STAGE16_2_REQUIRE_STAGE16_1` | `1` | Require Stage 16.1 PASS evidence. |
| `STAGE16_2_ENABLE` | `1` | Stage 16.2 case-definition request gate. |
| `STAGE16_2_ONE_FIBRE_FSI_ENABLE` | `1` | One-fibre case request for definition/check only. |
| `STAGE16_2_NPTS` | `8` | Number of fibre control points. |
| `STAGE16_2_MIN_WALL_CLEARANCE` | `1.0e-2` | Minimum normalized wall/domain clearance. |
| `STAGE16_2_MIN_POINT_SPACING` | `1.0e-6` | Minimum control-point spacing. |
| `STAGE16_2_MAX_INITIAL_VELOCITY` | `1.0e-8` | Maximum allowed initial velocity magnitude. |
| `STAGE16_2_MAX_INITIAL_ACCELERATION` | `1.0e-8` | Maximum allowed initial acceleration magnitude. |
| `STAGE16_2_MAX_STRUCTURE_UPDATE` | `1.0e-12` | Future structure-update safety bound. |
| `STAGE16_2_MAX_RHS_INCREMENT` | `1.0e-8` | Future RHS-increment safety bound. |
| `STAGE16_2_DIAGNOSTIC_ONLY` | `1` | Stage 16.2 must remain diagnostic-only. |

## Added case-definition API

`src/fibre_stage16_one_fibre_case.f90` provides:

- `stage16_one_fibre_case_reset()`
- `stage16_one_fibre_case_define_default()`
- `stage16_one_fibre_case_load_from_environment()`
- `stage16_one_fibre_case_validate()`
- `stage16_one_fibre_case_get_npts()`
- `stage16_one_fibre_case_get_position(index, position)`
- `stage16_one_fibre_case_get_velocity(index, velocity)`
- `stage16_one_fibre_case_get_acceleration(index, acceleration)`
- `stage16_one_fibre_case_write_diagnostics(unit)`

The API exposes case data for future stages but does not connect it to the
production time loop.

## Audit coverage

Stage 16.2 verifies:

1. Stage 16.2 wrapper/doc/source/helper/check target and minimal build
   registration exist.
2. Stage 14 and Stage 15 closure files exist when required.
3. Stage 16.1 config PASS evidence exists when required.
4. The default one-fibre case is valid, finite, exactly one fibre, and has valid
   `N_f`, spacing, wall/domain clearance, velocity bounds, and acceleration
   bounds.
5. Invalid `N_f`, geometry, velocity, acceleration, wall/contact request, and
   multi-fibre request are rejected.
6. Stage 16.2 does not insert production hooks into `src/xcompact3d.f90`.
7. Stage 16.2 does not activate structure advance, RHS modification,
   bending/tension solves, wall/contact handling beyond rejection flags,
   multi-fibre logic beyond rejection flags, pressure/projection/Poisson/RK3 or
   channel-forcing changes, or legacy IBM forcing outside the approved chain.
8. Stage 11.5, Stage 13.6, Stage 14.5, Stage 15.1-15.11, Stage 16.0, and Stage
   16.1 diagnostics remain present.
9. Rank0-safe diagnostic writing is preserved for Stage 11, Stage 13, Stage 14,
   Stage 15, and Stage 16 diagnostics.
10. No Stage 16.2 script requires ripgrep without a grep fallback.

## Required summary fields

`stage16_outputs/fibre_stage16_2_one_fibre_case_definition.dat` contains at
least:

```text
stage16_2_requested_status
one_fibre_case_definition_status
one_fibre_count_status
npts
npts_valid_status
position_finite_status
velocity_finite_status
acceleration_finite_status
initial_velocity_bound_status
initial_acceleration_bound_status
min_point_spacing
point_spacing_status
min_wall_clearance
wall_clearance_status
domain_containment_status
invalid_npts_rejection_status
invalid_geometry_rejection_status
invalid_velocity_rejection_status
invalid_acceleration_rejection_status
wall_contact_rejection_status
multifibre_rejection_status
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
final_status
```

`final_status 1` is required for PASS.

## Stage 16.2 boundary

Stage 16.2 intentionally adds only case-definition source/check code, a wrapper,
a parser/audit helper, documentation, and minimal standalone build registration.
It does not add one-fibre FSI physics, compute `F_fs_cand`, advance the
structure, solve bending/tension, handle wall/contact interactions, enable
multi-fibre operation, modify fluid RHS, modify pressure/projection/Poisson/RK3
or channel forcing, or activate legacy IBM forcing outside the approved
Stage 11-14 chain.
