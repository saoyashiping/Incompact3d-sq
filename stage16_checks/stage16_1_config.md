# Stage 16.1 configuration / global switches

Stage 16.1 is a configuration-only stage for the first controlled one-flexible-
fibre channel-DNS FSI path. It defines guarded, safe default-off controls that
future Stage 16 steps can use, but it does not compute feedback force, advance
the structure, alter the production time loop, inject RHS forcing, or modify any
fluid solver behavior.

The future coupling context remains:

```text
Stage 15 controlled structure-state update
-> Stage 12 feedback-force candidate
-> Stage 13 Eulerian force-density candidate
-> Stage 14 controlled RHS diagnostic/injection path
```

Stage 16.1 only configures switches and validates that this boundary is still
closed.

## Manual command

Run exactly:

```bash
bash stage16_checks/run_stage16_1_config.sh
```

The wrapper builds and runs only the standalone Stage 16.1 configuration check
target, then writes:

```text
stage16_outputs/fibre_stage16_1_config.dat
```

A successful run prints:

```text
STAGE 16.1 CONFIG VERDICT: PASS
STAGE 16.1 FINAL VERDICT: PASS
```

A failed run prints `FAIL` and explicit failure reasons.

## Safe defaults and environment controls

| Variable | Default | Purpose |
| --- | ---: | --- |
| `BUILD_DIR` | `build_stage9` | CMake build directory. |
| `DECOMP2D_ROOT` | empty | Optional CMake prefix path when configuring a missing build directory. |
| `MPIEXEC` | `mpirun` | Preserved for consistency with closed-stage wrappers. |
| `MPIEXEC_FLAGS` | empty | Preserved for consistency with closed-stage wrappers. |
| `STAGE16_1_RUN_STAGE16_0` | `0` | Optional prerequisite Stage 16.0 rerun. |
| `STAGE16_1_REQUIRE_STAGE14_CLOSED` | `1` | Require `stage14_checks/STAGE14_CLOSED.md`. |
| `STAGE16_1_REQUIRE_STAGE15_CLOSED` | `1` | Require `stage15_checks/STAGE15_CLOSED.md`. |
| `STAGE16_1_REQUIRE_STAGE16_0` | `1` | Require Stage 16.0 preflight PASS evidence. |
| `STAGE16_1_ENABLE` | `0` | Stage 16 master request flag; default off. |
| `STAGE16_1_ONE_FIBRE_FSI_ENABLE` | `0` | One-fibre FSI request flag; default off. |
| `STAGE16_1_STRUCTURE_ADVANCE_ENABLE` | `0` | Structure advance request flag; Stage 16.1 rejects activation. |
| `STAGE16_1_TWO_WAY_RHS_ENABLE` | `0` | Two-way RHS request flag; Stage 16.1 rejects activation. |
| `STAGE16_1_DIAGNOSTIC_ONLY` | `1` | Diagnostic-only guard; must remain enabled in Stage 16.1. |
| `STAGE16_1_FEEDBACK_ALPHA` | `1.0` | Future feedback gain; parsed and bounded only. |
| `STAGE16_1_LAMBDA` | `0.0` | Future RHS gain; parsed and bounded only. |
| `STAGE16_1_MAX_STRUCTURE_UPDATE` | `1.0e-12` | Future structure-update safety bound. |
| `STAGE16_1_MAX_RHS_INCREMENT` | `1.0e-8` | Future RHS-increment safety bound. |

All physics-impacting controls default to disabled. Invalid numeric values,
negative safety bounds, disabled diagnostic-only mode, structure-advance
activation, and two-way RHS activation fail closed in Stage 16.1. A one-fibre FSI
request is accepted only while diagnostic-only mode remains enabled.

## Added configuration API

`src/fibre_stage16_config.f90` provides a guarded configuration layer with these
public routines:

- `stage16_config_reset()`
- `stage16_config_load_from_environment()`
- `stage16_config_validate()`
- `stage16_is_enabled()`
- `stage16_one_fibre_fsi_requested()`
- `stage16_structure_advance_requested()`
- `stage16_two_way_rhs_requested()`
- `stage16_diagnostic_only_requested()`
- `stage16_get_lambda()`
- `stage16_get_feedback_alpha()`
- `stage16_get_max_structure_update()`
- `stage16_get_max_rhs_increment()`
- `stage16_config_write_diagnostics(unit)`

The module is not inserted into `xcompact3d.f90` and does not alter production
DNS numerics.

## Audit coverage

The Stage 16.1 wrapper/helper verifies:

1. The Stage 16.1 wrapper, documentation, source module, standalone check target,
   and build registration exist.
2. Stage 16.0 PASS evidence is present when `STAGE16_1_REQUIRE_STAGE16_0=1`.
3. Stage 14 and Stage 15 closure files are present when required.
4. Default-off behavior is preserved.
5. Environment override parsing is safe and diagnostic-only.
6. Invalid numeric values and invalid flag combinations are rejected.
7. Stage 16.1 does not insert a production hook into `src/xcompact3d.f90`.
8. Stage 16.1 does not activate structure advance, RHS modification,
   bending/tension solves, wall/contact handling, multi-fibre logic,
   pressure/projection/Poisson/RK3/channel-forcing changes, or legacy IBM
   forcing outside the approved chain.
9. Stage 14 small-lambda, Stage 11.5, Stage 13.6, Stage 14.5, rank0-safe, and
   Stage 13 sampling-repair protections remain intact.
10. Stage 15.1 through Stage 15.11 scripts/docs remain present and Stage 15.11
    explicit-failure protection remains active.
11. No Stage 16.1 wrapper path requires ripgrep without a grep fallback.

## Required summary fields

`stage16_outputs/fibre_stage16_1_config.dat` contains at least:

```text
stage16_1_requested_status
default_off_status
env_override_status
master_enable_status
one_fibre_fsi_enable_status
structure_advance_enable_status
two_way_rhs_enable_status
diagnostic_only_status
feedback_alpha
lambda_value
max_structure_update
max_rhs_increment
invalid_numeric_rejection_status
invalid_flag_combination_rejection_status
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
final_status
```

`final_status 1` is required for Stage 16.1 PASS.

## Stage 16.1 boundary

Stage 16.1 intentionally adds only the global-switch configuration layer,
standalone configuration diagnostics, build registration for that standalone
check, and Stage 16.1 wrapper/docs/helper files. It does not add one-fibre FSI
physics, compute `F_fs_cand`, advance the structure, solve bending/tension,
handle wall/contact interactions, enable multi-fibre logic, modify fluid RHS,
modify pressure/projection/Poisson/RK3/channel forcing, or activate legacy IBM
forcing outside the approved Stage 11-14 chain.

## Stage 16.1 audit repair note

Stage 16.1 accepts the already-passed Stage 16.0 preflight evidence in a fresh unpacked tree when `STAGE16_1_ACCEPT_STAGE16_0_CLOSED_EVIDENCE=1`.  In that mode the Stage 16.0 wrapper, helper, documentation, `STAGE14_CLOSED.md`, and `STAGE15_CLOSED.md` are treated as the closed preflight evidence if the runtime `stage16_outputs/fibre_stage16_0_preflight_closure_integrity.dat` file is not present.

The Stage 16.1 helper must distinguish real executable commands from negative audit strings.  A literal `rg` inside a quoted regular expression is not an `rg` dependency.  A literal `stage14_get_injection_gain() == 0.0` inside closure documentation is not a Stage 14 hook-registration regression.  The legitimate Stage 13.5 conservation-sign audit is also not the old production force-density candidate name; only an actual reappearance of `stage13_5_production_force_density_candidate` in production force-density evidence is a failure.
