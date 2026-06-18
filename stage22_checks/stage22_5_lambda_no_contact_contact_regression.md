# Stage 22.5 lambda/no-contact/contact regression

Stage 22.5 is helper-only. It may compute wall contact force and fibre-fibre collision force candidates only inside the new Stage 22.5 helper, compare `lambda_contact` and `lambda_fsi` gated helper responses, and must not commit candidate X/V/A to production structure state, apply contact/collision force to production structure advance, spread contact/collision force to production RHS, call Stage 14 RHS injection, run production DNS, run MPI, build, or modify production source or closed-stage files.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, Stage 22.0 PASS, Stage 22.1 PASS, Stage 22.2 PASS, Stage 22.3 PASS, and Stage 22.4 PASS are accepted through available evidence or source-only closure acceptance. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 through Stage 22.4 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/I/O files remain untouched.

## Helper-only case settings

* Nx = 16
* Ny = 16
* Nz = 16
* dx = 1.0 / 16.0
* dy = 1.0 / 16.0
* dz = 1.0 / 16.0
* no production DNS
* dt = 1.0e-5
* n_steps = 20
* helper_only = true
* actual DNS stepping = false
* production structure commit = false
* n_fibre_values = 1,2
* n_point_per_fibre = 64
* component_dim = 3
* fibre_radius = 0.01
* fibre_length = 0.40
* rho_tilde = 1.0
* gamma = 1.0e-5
* c_fs = 1.0
* k_wall = 1.0e2
* k_collision = 1.0e2
* damping_ratio = 0.2
* m_eff = rho_tilde * ds
* c_wall = 2 * damping_ratio * sqrt(k_wall * m_eff)
* c_collision = 2 * damping_ratio * sqrt(k_collision * m_eff)
* lambda_fsi_values = 0.0, 1.0e-6
* lambda_contact_values = 0.0, 1.0e-6, 1.0

## Safety limits

* max_penetration_allowed = 1.0e-4
* max_contact_force_norm = 1.0e3
* max_collision_force_norm = 1.0e3
* max_total_force_norm = 1.0e3
* max_contact_acceleration = 1.0e3
* max_structure_step_displacement_fraction = 0.1
* segment_length_residual_limit = 1.0e-6
* action_reaction_tol = 1.0e-12
* audit_tol = 1.0e-12
* zero_tol = 1.0e-14
* scaling_relative_tol = 1.0e-10

## Stage 22.5 formulas

```text
u_relative_candidate = u_interp_candidate - V_fibre
F_fs_candidate = c_fs * u_relative_candidate
F_on_structure_from_fluid_candidate = -F_fs_candidate
F_on_fluid_from_structure_candidate = +F_fs_candidate
rhs_delta_candidate = lambda_fsi * f_eulerian_candidate
```

If `lambda_fsi = 0.0`, `rhs_delta_candidate` must be exactly zero within `zero_tol`.

```text
delta_wall = max(0, -g_wall)
F_wall_raw = k_wall * delta_wall * n_wall - c_wall * v_n_minus_wall * n_wall
F_wall_candidate = lambda_contact * F_wall_raw if delta_wall > 0 else 0
delta_ff = max(0, -g_ff)
F_collision_raw_i = k_collision * delta_ff * n_ij - c_collision * v_n_minus_collision * n_ij
F_collision_candidate_i = lambda_contact * F_collision_raw_i if delta_ff > 0 else 0
F_collision_candidate_j = -F_collision_candidate_i
F_total_with_contact_candidate = F_bending_candidate + F_tension_candidate - F_fs_candidate + F_wall_candidate + F_collision_candidate
A_candidate = F_total_with_contact_candidate / rho_tilde
V_next_candidate = V_current + dt * A_candidate
X_next_candidate = X_current + dt * V_current + 0.5 * dt^2 * A_candidate
```

If `lambda_contact = 0.0`, wall/collision candidates must be exactly zero within `zero_tol` even under penetration/overlap. If `delta_wall = 0` or `delta_ff = 0`, contact/collision candidates must be zero for all `lambda_contact` values. Scaling rule: `norm(F_contact at lambda_contact=1.0e-6) / norm(F_contact at lambda_contact=1.0)` must equal `1.0e-6` within `scaling_relative_tol`, unless both norms are zero. No-contact invariance requires zero wall/collision force, invariant `F_total_with_contact_candidate`, and unchanged candidate update for all `lambda_contact` values. No candidate update may be committed to production structure state.

## Controlled regression cases

* Case A: one_fibre_no_contact_lambda_contact_0
* Case B: one_fibre_no_contact_lambda_contact_1
* Case C: one_fibre_wall_penetration_lambda_contact_0
* Case D: one_fibre_wall_penetration_lambda_contact_small
* Case E: one_fibre_wall_penetration_lambda_contact_1
* Case F: two_fibre_overlap_lambda_contact_0
* Case G: two_fibre_overlap_lambda_contact_small
* Case H: two_fibre_overlap_lambda_contact_1
* Case I: lambda_fsi_0_rhs_noop
* Case J: lambda_fsi_small_rhs_bounded

## Required validation groups

Gate and evidence group: Stage 22.5 requested; Stage 22.5 lambda/no-contact/contact regression enabled; Stage 22.4 evidence accepted; Stage 22.3 evidence accepted; Stage 22.2 evidence accepted; Stage 22.1 evidence accepted; Stage 22.0 evidence accepted; Stage 20 closure accepted; Stage 21 closure accepted; source-only closure acceptance preserved; missing old outputs allowed; no previous stage rerun; helper-only mode enabled.

Lambda gate group: lambda_fsi ladder valid; lambda_contact ladder valid; lambda_fsi zero gives RHS no-op; lambda_fsi small gives bounded helper RHS candidate; lambda_contact zero gives wall no-op under penetration; lambda_contact zero gives collision no-op under overlap; lambda_contact small response scales from lambda_contact=1; lambda_contact does not create force without penetration/overlap.

No-contact invariance group: no-contact geometry has zero wall force for all lambda_contact; no-overlap geometry has zero collision force for all lambda_contact; no-contact/no-overlap F_total equals contact-disabled F_total; no-contact/no-overlap candidate update unchanged.

Contact/collision response group: wall penetration force inward; fibre-fibre overlap force repulsive; fibre-fibre action-reaction residual bounded; contact/collision force finite and bounded; contact/collision acceleration bounded; contact/collision energy nonnegative; damping power nonpositive; no attractive wall force; no duplicate pair force; no self-pair force.

Candidate update group: F_total_with_contact_candidate finite; A_candidate finite; V_next_candidate finite; X_next_candidate finite; structure step displacement fraction bounded; segment length residual bounded; candidate update not committed to production state.

Isolation group: production structure advance disabled; production structure state unchanged; RHS coupling disabled except helper-local rhs_delta_candidate; contact/collision not spread to production RHS; Stage 14 RHS injection disabled; production RHS update disabled; IBM forcing modification disabled; production DNS disabled; MPI disabled; build disabled; production restart/statistics/visualization I/O disabled; production multifibre disabled; no source modification; no closed-stage modification.

## Safe defaults

* STAGE22_5_ENABLE=1
* STAGE22_5_LAMBDA_NO_CONTACT_CONTACT_REGRESSION_ENABLE=1
* STAGE22_5_REQUIRE_STAGE22_4_PASS=1
* STAGE22_5_REQUIRE_STAGE22_3_PASS=1
* STAGE22_5_REQUIRE_STAGE22_2_PASS=1
* STAGE22_5_REQUIRE_STAGE22_1_PASS=1
* STAGE22_5_REQUIRE_STAGE22_0_PASS=1
* STAGE22_5_REQUIRE_STAGE20_CLOSED=1
* STAGE22_5_REQUIRE_STAGE21_CLOSED=1
* STAGE22_5_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_5_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_5_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_5_DIAGNOSTIC_ONLY=1
* STAGE22_5_FAIL_CLOSED=1
* STAGE22_5_HELPER_ONLY=1
* STAGE22_5_BUILD_ALLOWED=0
* STAGE22_5_PRODUCTION_DNS_ALLOWED=0
* STAGE22_5_ACTUAL_MPI_ALLOWED=0
* STAGE22_5_PRODUCTION_STRUCTURE_COMMIT_ALLOWED=0
* STAGE22_5_PRODUCTION_STRUCTURE_ADVANCE_ALLOWED=0
* STAGE22_5_RHS_COUPLING_ENABLE=0
* STAGE22_5_CONTACT_TO_RHS_ENABLE=0
* STAGE22_5_COLLISION_TO_RHS_ENABLE=0
* STAGE22_5_STAGE14_RHS_INJECTION_ALLOWED=0
* STAGE22_5_PRODUCTION_RHS_UPDATE_ALLOWED=0
* STAGE22_5_PRODUCTION_RESTART_IO_ALLOWED=0
* STAGE22_5_PRODUCTION_STATISTICS_IO_ALLOWED=0
* STAGE22_5_PRODUCTION_VISUALIZATION_IO_ALLOWED=0
* STAGE22_5_PRODUCTION_MULTIFIBRE_ENABLE=0

## Next stage

Stage 22.6: single-fibre channel FSI micro-case.
