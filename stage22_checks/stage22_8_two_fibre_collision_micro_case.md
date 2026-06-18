# Stage 22.8 two-fibre collision micro-case

Stage 22.8 is the first real DNS micro-case with fibre-fibre collision force enabled.  It is intentionally cautious: np=1 only, no restart testing, no mesh/time-step refinement, no np=2/4 consistency check, and no production source or schema modification.

## Source-only and prior-stage acceptance

Stage 22.8 accepts Stage 20 closure, Stage 21 closure, and Stage 22.0 through Stage 22.7 PASS evidence through available committed artifacts or source-only closure acceptance.  Missing old outputs are allowed and previous stages are not rerun.

## Controlled G1 channel target

* grid label: G1
* Nx = 32
* Ny = 33
* Nz = 32
* channel case
* channel_height = 1.0
* y_min = 0.0
* y_max = 1.0
* dt = 5.0e-5
* n_steps = 300
* final_time = 0.015
* CFL_max_limit = 0.3
* contact_CFL_limit = 0.2
* controlled real DNS micro-case = true

## Parallel and case gates

* np = 1 only
* np=2 forbidden
* np=4 forbidden
* n_fibre = 2
* n_point_per_fibre = 64
* component_dim = 3
* fibre_radius = 0.01
* fibre_radius_sum = 0.02
* fibre_length = 0.40
* rho_tilde = 1.0
* gamma = 1.0e-5
* c_fs = 1.0
* initial fibre-fibre gap target = 1.0e-3
* initial severe overlap forbidden
* initial wall penetration absent or bounded by max_penetration_allowed

## FSI and contact/collision gates

* Stage 20 two-way FSI pathway enabled
* lambda_fsi = 1.0e-4
* lambda_contact = 1.0
* fibre-fibre collision force enabled
* fibre-fibre collision candidate enabled
* collision force application enabled
* wall contact force enabled only as safety guard
* uncontrolled production multifibre disabled
* no restart/statistics/visualization schema modification

## Fibre-fibre collision formula

For fibres i and j:

```text
d_ff = minimum point/segment distance between fibre i and fibre j
g_ff = d_ff - (fibre_radius_i + fibre_radius_j)
delta_ff = max(0, -g_ff)
n_ij = (X_i_star - X_j_star) / |X_i_star - X_j_star|
v_n = (V_i_star - V_j_star) dot n_ij
v_n_minus = min(v_n, 0)
F_i_collision = lambda_contact * (k_collision * delta_ff * n_ij - c_collision * v_n_minus * n_ij) if delta_ff > 0 else 0
F_j_collision = -F_i_collision
```

The normal must be well-defined.  Collision force is repulsive, action-reaction residual is bounded, duplicate pair force is forbidden, and self-pair force is forbidden.

## Wall-contact safety guard

Wall contact is enabled only as a safety guard.  It must remain finite and bounded if active, must not be attractive, and must not dominate the two-fibre collision test unless wall penetration occurs.

## Safe defaults

* STAGE22_8_ENABLE=1
* STAGE22_8_TWO_FIBRE_COLLISION_MICRO_CASE_ENABLE=1
* STAGE22_8_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_8_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_8_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_8_CONTROLLED_MICRO_CASE_ALLOWED=1
* STAGE22_8_BUILD_ALLOWED=1
* STAGE22_8_PRODUCTION_DNS_ALLOWED=1
* STAGE22_8_ACTUAL_MPI_ALLOWED=1
* STAGE22_8_NP=1
* STAGE22_8_NP2_ALLOWED=0
* STAGE22_8_NP4_ALLOWED=0
* STAGE22_8_NX=32
* STAGE22_8_NY=33
* STAGE22_8_NZ=32
* STAGE22_8_DT=5.0e-5
* STAGE22_8_N_STEPS=300
* STAGE22_8_FINAL_TIME=0.015
* STAGE22_8_CFL_MAX_LIMIT=0.3
* STAGE22_8_CONTACT_CFL_LIMIT=0.2
* STAGE22_8_N_FIBRE=2
* STAGE22_8_N_POINT_PER_FIBRE=64
* STAGE22_8_COMPONENT_DIM=3
* STAGE22_8_FIBRE_RADIUS=0.01
* STAGE22_8_FIBRE_LENGTH=0.40
* STAGE22_8_LAMBDA_FSI=1.0e-4
* STAGE22_8_LAMBDA_CONTACT=1.0
* STAGE22_8_K_COLLISION=1.0e2
* STAGE22_8_K_WALL=1.0e2
* STAGE22_8_DAMPING_RATIO=0.2
* STAGE22_8_INITIAL_FIBRE_FIBRE_GAP_TARGET=1.0e-3
* STAGE22_8_COLLISION_FORCE_ENABLE=1
* STAGE22_8_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE=1
* STAGE22_8_COLLISION_FORCE_APPLY_ENABLE=1
* STAGE22_8_WALL_CONTACT_FORCE_ENABLE=1
* STAGE22_8_WALL_CONTACT_SAFETY_ENABLE=1
* STAGE22_8_PRODUCTION_MULTIFIBRE_ENABLE=0
* STAGE22_8_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED=0
* STAGE22_8_PRODUCTION_RESTART_TEST_ALLOWED=0
* STAGE22_8_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_8_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_8_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED=0

## Required PASS evidence

The audit writes `stage22_outputs/fibre_stage22_8_two_fibre_collision_micro_case.dat` and must contain:

```text
STAGE 22.8 TWO-FIBRE COLLISION MICRO-CASE VERDICT: PASS
STAGE 22.8 FINAL VERDICT: PASS
```

## Next stage

Stage 22.9: mesh/time-step sensitivity check.
