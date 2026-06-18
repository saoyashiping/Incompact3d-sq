# Stage 22.3 fibre-fibre collision force candidate helper test

Stage 22 title: final integrated validation and production-readiness closure.
Stage 22.3 title: fibre-fibre collision force candidate helper test.

Stage 22.3 is helper-only. It may compute a fibre-fibre collision force candidate only inside the new Stage 22.3 helper. It does not compute wall contact force, apply collision/contact force to production structure advance, spread collision/contact force to RHS, call Stage 14 RHS injection, run production DNS, run MPI, build, or modify production source or closed-stage files.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, Stage 22.0 PASS, Stage 22.1 PASS, and Stage 22.2 PASS are accepted through available evidence or source-only closure acceptance. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 through Stage 22.2 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/I/O files remain untouched.

## Helper-only case settings

Synthetic helper grid metadata:

* Nx = 16
* Ny = 16
* Nz = 16
* dx = 1.0 / 16.0
* dy = 1.0 / 16.0
* dz = 1.0 / 16.0
* no production DNS

Time metadata:

* dt = 1.0e-5
* n_steps = 1
* helper_only = true
* actual DNS stepping = false

Fibre metadata:

* n_fibre = 2
* n_point_per_fibre = 64
* component_dim = 3
* fibre_radius_0 = 0.01
* fibre_radius_1 = 0.01
* fibre_radius_sum = 0.02
* fibre_length = 0.40
* rho_tilde = 1.0

Collision candidate parameters:

* k_collision = 1.0e2
* damping_ratio = 0.2
* m_eff = rho_tilde * ds
* c_collision = 2 * damping_ratio * sqrt(k_collision * m_eff)
* lambda_contact = 1.0
* lambda_fsi = 0.0

Safety limits:

* max_penetration_allowed = 1.0e-4
* max_collision_force_norm = 1.0e3
* max_collision_acceleration = 1.0e3
* action_reaction_tol = 1.0e-12
* audit_tol = 1.0e-12
* zero_tol = 1.0e-14

## Fibre-fibre distance and collision formulas

```text
d_ff = minimum point/segment distance between fibre i and fibre j
g_ff = d_ff - (fibre_radius_i + fibre_radius_j)
delta_ff = max(0, -g_ff)
X_i_star
X_j_star
n_ij = (X_i_star - X_j_star) / |X_i_star - X_j_star|
V_rel = V_i_star - V_j_star
v_n = V_rel · n_ij
v_n_minus = min(v_n, 0)
F_elastic_i = k_collision * delta_ff * n_ij
F_damping_i = - c_collision * v_n_minus * n_ij
F_i_candidate = lambda_contact * (F_elastic_i + F_damping_i) if delta_ff > 0 else 0
F_j_candidate = -F_i_candidate
E_collision = 0.5 * k_collision * delta_ff^2
P_damping = F_damping_i · V_rel
```

If `|X_i_star - X_j_star|` is too small to define a normal, Stage 22.3 must fail closed rather than inventing an arbitrary normal. If `delta_ff = 0`, both candidates must be exactly zero within tolerance. Near-contact but non-overlapping geometry must not create force. Collision force must be repulsive, pair action-reaction must satisfy `F_i_candidate + F_j_candidate = 0`, collision force must not be duplicated, and self-pairs must be excluded. For active damping under approach, `P_damping <= 0`.

## Controlled helper cases

* Case A: separated_safe; `g_ff > warning_gap`; `delta_ff = 0`; expected `F_i_candidate_norm = 0` and `F_j_candidate_norm = 0`.
* Case B: near_contact_no_overlap; `g_ff = 1.0e-3`; `delta_ff = 0`; relative velocity directed toward each other; expected zero force on both fibres.
* Case C: small_overlap_order_ij; `g_ff = -1.0e-5`; `delta_ff = 1.0e-5`; relative velocity directed toward each other; expected repulsive force, action-reaction bounded, nonnegative collision energy, nonpositive damping power, and bounded force.
* Case D: small_overlap_reversed_order_ji; same physical geometry as Case C with candidate pair order reversed; expected canonical ordering gives physically equivalent force and no duplicate pair.

Stage 22.3 computes only fibre-fibre collision force candidate inside this helper. Stage 22.3 does not compute wall contact force, does not compute RHS spreading, and does not update production structure state.

## Gate and evidence group

* Stage 22.3 requested
* Stage 22.3 fibre-fibre collision force candidate helper enabled
* Stage 22.2 evidence accepted
* Stage 22.1 evidence accepted
* Stage 22.0 evidence accepted
* Stage 20 closure accepted
* Stage 21 closure accepted
* source-only closure acceptance preserved
* missing old outputs allowed
* no previous stage rerun
* helper-only mode enabled

## Geometry/gap group

* fibre geometry finite
* fibre radii valid
* fibre radius sum valid
* closest point/segment distance finite
* signed gap formula valid
* penetration depth formula valid
* collision normal valid for overlap cases
* relative velocity finite
* no self-pair
* no duplicate pair

## Fibre-fibre collision candidate formula group

* k_collision valid
* damping ratio valid
* m_eff valid
* c_collision formula valid
* elastic collision force formula valid
* damping collision force formula valid
* total collision force formula valid
* lambda_contact applied correctly
* pair action-reaction formula valid

## Controlled case group

* separated case gives zero force
* near-contact no-overlap case gives zero force
* small-overlap case gives repulsive force
* reversed-order case gives equivalent canonical result
* action-reaction residual bounded
* no duplicate pair force
* no self-pair force
* collision force finite
* collision force bounded
* collision acceleration bounded
* collision energy nonnegative
* damping power nonpositive

## Isolation group

* wall contact force disabled
* wall contact force candidate disabled
* contact/collision force application disabled
* collision not added to structure total force
* structure advance disabled
* RHS coupling disabled
* collision not spread to RHS
* Stage 14 RHS injection disabled
* production RHS update disabled
* production DNS disabled
* MPI disabled
* build disabled
* production restart/statistics/visualization I/O disabled
* production multifibre disabled
* no production collision physics activation

## Safe defaults

* STAGE22_3_ENABLE=1
* STAGE22_3_FIBRE_COLLISION_FORCE_CANDIDATE_HELPER_ENABLE=1
* STAGE22_3_REQUIRE_STAGE22_2_PASS=1
* STAGE22_3_REQUIRE_STAGE22_1_PASS=1
* STAGE22_3_REQUIRE_STAGE22_0_PASS=1
* STAGE22_3_REQUIRE_STAGE20_CLOSED=1
* STAGE22_3_REQUIRE_STAGE21_CLOSED=1
* STAGE22_3_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_3_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_3_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_3_DIAGNOSTIC_ONLY=1
* STAGE22_3_FAIL_CLOSED=1
* STAGE22_3_HELPER_ONLY=1
* STAGE22_3_BUILD_ALLOWED=0
* STAGE22_3_PRODUCTION_DNS_ALLOWED=0
* STAGE22_3_ACTUAL_MPI_ALLOWED=0
* STAGE22_3_NX=16
* STAGE22_3_NY=16
* STAGE22_3_NZ=16
* STAGE22_3_DT=1.0e-5
* STAGE22_3_N_STEPS=1
* STAGE22_3_N_FIBRE=2
* STAGE22_3_N_POINT_PER_FIBRE=64
* STAGE22_3_COMPONENT_DIM=3
* STAGE22_3_FIBRE_RADIUS_0=0.01
* STAGE22_3_FIBRE_RADIUS_1=0.01
* STAGE22_3_FIBRE_LENGTH=0.40
* STAGE22_3_RHO_TILDE=1.0
* STAGE22_3_K_COLLISION=1.0e2
* STAGE22_3_DAMPING_RATIO=0.2
* STAGE22_3_LAMBDA_CONTACT=1.0
* STAGE22_3_LAMBDA_FSI=0.0
* STAGE22_3_WALL_CONTACT_FORCE_ENABLE=0
* STAGE22_3_WALL_CONTACT_FORCE_CANDIDATE_ENABLE=0
* STAGE22_3_CONTACT_FORCE_APPLY_ENABLE=0
* STAGE22_3_STRUCTURE_ADVANCE_ENABLE=0
* STAGE22_3_RHS_COUPLING_ENABLE=0
* STAGE22_3_COLLISION_TO_RHS_ENABLE=0
* STAGE22_3_STAGE14_RHS_INJECTION_ALLOWED=0
* STAGE22_3_PRODUCTION_RHS_UPDATE_ALLOWED=0
* STAGE22_3_PRODUCTION_RESTART_IO_ALLOWED=0
* STAGE22_3_PRODUCTION_STATISTICS_IO_ALLOWED=0
* STAGE22_3_PRODUCTION_VISUALIZATION_IO_ALLOWED=0
* STAGE22_3_PRODUCTION_MULTIFIBRE_ENABLE=0
* STAGE22_3_MAX_PENETRATION_ALLOWED=1.0e-4
* STAGE22_3_MAX_COLLISION_FORCE_NORM=1.0e3
* STAGE22_3_MAX_COLLISION_ACCELERATION=1.0e3
* STAGE22_3_ACTION_REACTION_TOL=1.0e-12
* STAGE22_3_AUDIT_TOL=1.0e-12
* STAGE22_3_ZERO_TOL=1.0e-14

## Next stage

Stage 22.4: contact force into structure candidate.
