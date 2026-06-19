# Stage 22.4 contact force into structure candidate

Stage 22.4 is helper-only. It may compute `F_wall_candidate` and `F_collision_candidate` only inside the new Stage 22.4 helper and may add these helper-local candidate forces into `F_total_candidate` only inside the helper. It must not commit candidate X/V/A to production structure state, apply contact/collision force to production structure advance, spread contact/collision force to RHS, call Stage 14 RHS injection, run production DNS, run MPI, build, or modify production source or closed-stage files.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, Stage 22.0 PASS, Stage 22.1 PASS, Stage 22.2 PASS, and Stage 22.3 PASS are accepted through available evidence or source-only closure acceptance. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 through Stage 22.3 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/I/O files remain untouched.

## Helper-only case settings

* Nx = 16
* Ny = 16
* Nz = 16
* dx = 1.0 / 16.0
* dy = 1.0 / 16.0
* dz = 1.0 / 16.0
* no production DNS
* dt = 1.0e-5
* n_steps = 10
* helper_only = true
* actual DNS stepping = false
* production structure commit = false
* n_fibre = 2
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
* lambda_fsi = 1.0e-6
* lambda_contact = 1.0

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

## Stage 22.4 helper-local force balance

```text
F_total_without_contact_candidate(q) = F_bending_candidate(q) + F_tension_candidate(q) - F_fs_candidate(q)
F_total_contact_candidate(q) = F_wall_candidate(q) + F_collision_candidate(q)
F_total_with_contact_candidate(q) = F_total_without_contact_candidate(q) + F_total_contact_candidate(q)
A_candidate(q) = F_total_with_contact_candidate(q) / rho_tilde
V_next_candidate(q) = V_current(q) + dt * A_candidate(q)
X_next_candidate(q) = X_current(q) + dt * V_current(q) + 0.5 * dt^2 * A_candidate(q)
```

No production structure force arrays may be modified. `X_next_candidate`, `V_next_candidate`, and `A_candidate` are diagnostic helper fields only. They must not be written to production structure state, used by production DNS, spread to RHS, or coupled into IBM.

## Contact/collision candidate formulas

Use the Stage 22.2 formula: `delta_wall = max(0, -g_wall)` and `F_wall_candidate = lambda_contact * (k_wall * delta_wall * n_wall - c_wall * v_n_minus * n_wall) if delta_wall > 0 else 0`. `F_wall` must be inward only.

Use the Stage 22.3 formula: `delta_ff = max(0, -g_ff)`, `F_i_collision_candidate = lambda_contact * (k_collision * delta_ff * n_ij - c_collision * v_n_minus * n_ij) if delta_ff > 0 else 0`, and `F_j_collision_candidate = -F_i_collision_candidate`. Pair action-reaction must hold.

## Controlled helper cases

* Case A: no_contact_no_collision; no wall penetration; no fibre-fibre overlap; `F_wall_candidate = 0`; `F_collision_candidate = 0`; `F_total_with_contact_candidate = F_total_without_contact_candidate`.
* Case B: wall_small_penetration_only; lower wall small penetration with `g_wall_min = -1.0e-5`; no fibre-fibre overlap; active inward wall force; zero collision force; helper candidate gap remains bounded or improves.
* Case C: fibre_collision_small_overlap_only; fibre-fibre small overlap with `g_ff_min = -1.0e-5`; no wall penetration; active repulsive collision force; zero wall force; `F_i + F_j = 0`; helper candidate fibre-fibre gap remains bounded or improves.
* Case D: combined_wall_and_collision; one wall small penetration and one fibre-fibre small overlap; both candidate forces active; `F_total_with_contact_candidate` finite and bounded; candidate X/V/A finite; no production commit.

## Validation groups

Gate and evidence group: Stage 22.4 requested; Stage 22.4 contact force into structure candidate enabled; Stage 22.3 evidence accepted; Stage 22.2 evidence accepted; Stage 22.1 evidence accepted; Stage 22.0 evidence accepted; Stage 20 closure accepted; Stage 21 closure accepted; source-only closure acceptance preserved; missing old outputs allowed; no previous stage rerun; helper-only mode enabled.

Structure force candidate group: F_bending_candidate finite; F_tension_candidate finite; F_fs_candidate finite; F_total_without_contact_candidate finite; F_wall_candidate finite; F_collision_candidate finite; F_total_contact_candidate finite; F_total_with_contact_candidate finite; F_total formula valid; no contact/collision double counting.

Candidate update group: A_candidate finite; V_next_candidate finite; X_next_candidate finite; acceleration bounded; velocity increment bounded; position increment bounded; structure step displacement fraction bounded; segment length residual bounded; candidate update is not committed to production state.

Controlled case group: no-contact/no-collision case has no contact contribution; wall-only case gives inward wall contribution; collision-only case gives repulsive collision contribution; combined case remains finite and bounded; wall gap remains bounded or improves; fibre-fibre gap remains bounded or improves; contact/collision energy nonnegative; damping power nonpositive; collision action-reaction residual bounded; no attractive wall force; no duplicate pair force; no self-pair force.

Isolation group: production structure advance disabled; production structure state unchanged; RHS coupling disabled; contact/collision not spread to RHS; Stage 14 RHS injection disabled; production RHS update disabled; IBM forcing modification disabled; production DNS disabled; MPI disabled; build disabled; production restart/statistics/visualization I/O disabled; production multifibre disabled; no source modification; no closed-stage modification.

## Safe defaults

* STAGE22_4_ENABLE=1
* STAGE22_4_CONTACT_FORCE_INTO_STRUCTURE_CANDIDATE_ENABLE=1
* STAGE22_4_REQUIRE_STAGE22_3_PASS=1
* STAGE22_4_REQUIRE_STAGE22_2_PASS=1
* STAGE22_4_REQUIRE_STAGE22_1_PASS=1
* STAGE22_4_REQUIRE_STAGE22_0_PASS=1
* STAGE22_4_REQUIRE_STAGE20_CLOSED=1
* STAGE22_4_REQUIRE_STAGE21_CLOSED=1
* STAGE22_4_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_4_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_4_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_4_DIAGNOSTIC_ONLY=1
* STAGE22_4_FAIL_CLOSED=1
* STAGE22_4_HELPER_ONLY=1
* STAGE22_4_BUILD_ALLOWED=0
* STAGE22_4_PRODUCTION_DNS_ALLOWED=0
* STAGE22_4_ACTUAL_MPI_ALLOWED=0
* STAGE22_4_PRODUCTION_STRUCTURE_COMMIT_ALLOWED=0
* STAGE22_4_PRODUCTION_STRUCTURE_ADVANCE_ALLOWED=0
* STAGE22_4_RHS_COUPLING_ENABLE=0
* STAGE22_4_CONTACT_TO_RHS_ENABLE=0
* STAGE22_4_COLLISION_TO_RHS_ENABLE=0
* STAGE22_4_STAGE14_RHS_INJECTION_ALLOWED=0
* STAGE22_4_PRODUCTION_RHS_UPDATE_ALLOWED=0
* STAGE22_4_PRODUCTION_RESTART_IO_ALLOWED=0
* STAGE22_4_PRODUCTION_STATISTICS_IO_ALLOWED=0
* STAGE22_4_PRODUCTION_VISUALIZATION_IO_ALLOWED=0
* STAGE22_4_PRODUCTION_MULTIFIBRE_ENABLE=0

## Next stage

Stage 22.5: lambda/no-contact/contact regression.
