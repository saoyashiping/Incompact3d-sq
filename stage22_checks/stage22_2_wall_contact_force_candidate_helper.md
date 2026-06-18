# Stage 22.2 wall contact force candidate helper test

Stage 22 title: final integrated validation and production-readiness closure.
Stage 22.2 title: wall contact force candidate helper test.

Stage 22.2 is helper-only. It may compute a wall contact force candidate only inside the new Stage 22.2 helper. It does not compute fibre-fibre collision force, apply contact force to production structure advance, spread contact force to RHS, call Stage 14 RHS injection, run production DNS, run MPI, build, or modify production source or closed-stage files.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, Stage 22.0 PASS, and Stage 22.1 PASS are accepted through available evidence or source-only closure acceptance. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 through Stage 22.1 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/I/O files remain untouched.

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

Channel and fibre metadata:

* y_min = 0.0
* y_max = 1.0
* channel_height = 1.0
* n_fibre = 1
* n_point_per_fibre = 64
* component_dim = 3
* fibre_radius = 0.01
* fibre_length = 0.40
* rho_tilde = 1.0

Wall contact candidate parameters:

* k_wall = 1.0e2
* damping_ratio = 0.2
* m_eff = rho_tilde * ds
* c_wall = 2 * damping_ratio * sqrt(k_wall * m_eff)
* lambda_contact = 1.0
* lambda_fsi = 0.0

Safety limits:

* max_penetration_allowed = 1.0e-4
* max_contact_force_norm = 1.0e3
* max_contact_acceleration = 1.0e3
* audit_tol = 1.0e-12
* zero_tol = 1.0e-14

## Wall contact formulas

For each fibre point q:

```text
d_lower(q) = y_q - y_min
d_upper(q) = y_max - y_q
g_lower(q) = d_lower(q) - fibre_radius
g_upper(q) = d_upper(q) - fibre_radius
g_wall(q) = min(g_lower(q), g_upper(q))
nearest_wall(q) = lower if g_lower(q) <= g_upper(q), upper otherwise
n_lower = (0, 1, 0)
n_upper = (0, -1, 0)
delta_wall(q) = max(0, -g_wall(q))
v_n(q) = V_q · n_wall(q)
v_n_minus(q) = min(v_n(q), 0)
F_elastic(q) = k_wall * delta_wall(q) * n_wall(q)
F_damping(q) = - c_wall * v_n_minus(q) * n_wall(q)
F_wall_candidate(q) = lambda_contact * [F_elastic(q) + F_damping(q)] if delta_wall(q) > 0 else 0
E_wall(q) = 0.5 * k_wall * delta_wall(q)^2
P_damping(q) = F_damping(q) · V_q
```

If `delta_wall = 0`, `F_wall_candidate` must be exactly zero within tolerance. Near-wall but non-penetrating gaps must not create force. Lower-wall force y-component must be >= 0, upper-wall force y-component must be <= 0, and the wall force must never attract fibre into the wall. For active damping under approach, `P_damping(q) <= 0`.

## Controlled helper cases

* Case A: no_contact_lower_safe; nearest wall lower; `g_wall_min > warning_gap`; `delta_wall = 0`; velocity may be zero; expected `F_wall_candidate_norm = 0`.
* Case B: near_wall_no_penetration_lower; nearest wall lower; `g_wall_min = 1.0e-3`; `delta_wall = 0`; velocity directed toward lower wall; expected `F_wall_candidate_norm = 0`.
* Case C: small_penetration_lower; nearest wall lower; `g_wall_min = -1.0e-5`; `delta_wall = 1.0e-5`; velocity directed toward lower wall; expected inward y-component >= 0, bounded force, nonnegative energy, and nonpositive damping power.
* Case D: small_penetration_upper; nearest wall upper; `g_wall_min = -1.0e-5`; `delta_wall = 1.0e-5`; velocity directed toward upper wall; expected inward y-component <= 0, bounded force, nonnegative energy, and nonpositive damping power.

Stage 22.2 computes only wall contact force candidate inside this helper. Stage 22.2 does not compute fibre-fibre collision force, does not compute RHS spreading, and does not update production structure state.

## Gate and evidence group

* Stage 22.2 requested
* Stage 22.2 wall contact force candidate helper enabled
* Stage 22.1 evidence accepted
* Stage 22.0 evidence accepted
* Stage 20 closure accepted
* Stage 21 closure accepted
* source-only closure acceptance preserved
* missing old outputs allowed
* no previous stage rerun
* helper-only mode enabled

## Geometry/gap group

* channel bounds valid
* fibre radius valid
* fibre geometry finite
* wall distance finite
* wall signed gap formula valid
* nearest wall side valid
* penetration depth formula valid

## Wall contact candidate formula group

* k_wall valid
* damping ratio valid
* m_eff valid
* c_wall formula valid
* elastic force formula valid
* damping force formula valid
* total wall contact force formula valid
* lambda_contact applied correctly

## Controlled case group

* no-contact case gives zero force
* near-wall no-penetration case gives zero force
* lower small-penetration case gives inward force
* upper small-penetration case gives inward force
* no attractive wall force
* contact force finite
* contact force bounded
* contact acceleration bounded
* contact energy nonnegative
* damping power nonpositive

## Isolation group

* fibre-fibre collision force disabled
* fibre-fibre collision force candidate disabled
* contact force application disabled
* contact not added to structure total force
* structure advance disabled
* RHS coupling disabled
* contact not spread to RHS
* Stage 14 RHS injection disabled
* production RHS update disabled
* production DNS disabled
* MPI disabled
* build disabled
* production restart/statistics/visualization I/O disabled
* production multifibre disabled
* no production contact physics activation

## Safe defaults

* STAGE22_2_ENABLE=1
* STAGE22_2_WALL_CONTACT_FORCE_CANDIDATE_HELPER_ENABLE=1
* STAGE22_2_REQUIRE_STAGE22_1_PASS=1
* STAGE22_2_REQUIRE_STAGE22_0_PASS=1
* STAGE22_2_REQUIRE_STAGE20_CLOSED=1
* STAGE22_2_REQUIRE_STAGE21_CLOSED=1
* STAGE22_2_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_2_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_2_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_2_DIAGNOSTIC_ONLY=1
* STAGE22_2_FAIL_CLOSED=1
* STAGE22_2_HELPER_ONLY=1
* STAGE22_2_BUILD_ALLOWED=0
* STAGE22_2_PRODUCTION_DNS_ALLOWED=0
* STAGE22_2_ACTUAL_MPI_ALLOWED=0
* STAGE22_2_NX=16
* STAGE22_2_NY=16
* STAGE22_2_NZ=16
* STAGE22_2_DT=1.0e-5
* STAGE22_2_N_STEPS=1
* STAGE22_2_N_FIBRE=1
* STAGE22_2_N_POINT_PER_FIBRE=64
* STAGE22_2_COMPONENT_DIM=3
* STAGE22_2_Y_MIN=0.0
* STAGE22_2_Y_MAX=1.0
* STAGE22_2_FIBRE_RADIUS=0.01
* STAGE22_2_FIBRE_LENGTH=0.40
* STAGE22_2_RHO_TILDE=1.0
* STAGE22_2_K_WALL=1.0e2
* STAGE22_2_DAMPING_RATIO=0.2
* STAGE22_2_LAMBDA_CONTACT=1.0
* STAGE22_2_LAMBDA_FSI=0.0
* STAGE22_2_COLLISION_FORCE_ENABLE=0
* STAGE22_2_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE=0
* STAGE22_2_CONTACT_FORCE_APPLY_ENABLE=0
* STAGE22_2_STRUCTURE_ADVANCE_ENABLE=0
* STAGE22_2_RHS_COUPLING_ENABLE=0
* STAGE22_2_CONTACT_TO_RHS_ENABLE=0
* STAGE22_2_STAGE14_RHS_INJECTION_ALLOWED=0
* STAGE22_2_PRODUCTION_RHS_UPDATE_ALLOWED=0
* STAGE22_2_PRODUCTION_RESTART_IO_ALLOWED=0
* STAGE22_2_PRODUCTION_STATISTICS_IO_ALLOWED=0
* STAGE22_2_PRODUCTION_VISUALIZATION_IO_ALLOWED=0
* STAGE22_2_PRODUCTION_MULTIFIBRE_ENABLE=0
* STAGE22_2_MAX_PENETRATION_ALLOWED=1.0e-4
* STAGE22_2_MAX_CONTACT_FORCE_NORM=1.0e3
* STAGE22_2_MAX_CONTACT_ACCELERATION=1.0e3
* STAGE22_2_AUDIT_TOL=1.0e-12
* STAGE22_2_ZERO_TOL=1.0e-14

## Next stage

Stage 22.3: fibre-fibre collision force candidate helper test.
