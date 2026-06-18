# Stage 22.7 single-fibre near-wall contact micro-case

Stage 22.7 is the first real DNS micro-case with wall-contact force enabled: one near-wall flexible fibre with Stage 20 FSI plus Stage 22 wall-contact force, fibre-fibre collision disabled, production multifibre disabled, and `np=1` only. Controlled build and controlled np=1 DNS execution are allowed only through the new Stage 22.7 wrapper. Stage 22.7 must not run np=2/np=4, perform restart testing, perform mesh/time-step refinement, perform np consistency checks, modify production source, or modify DNS-core, IBM, pressure projection, Poisson, RK3/channel forcing, restart/statistics/visualization schemas, or closed-stage files.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, and Stage 22.0 through Stage 22.6 PASS evidence are accepted through available evidence or source-only closure acceptance. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 through Stage 22.6 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/restart/statistics/visualization schema files remain untouched.

## Controlled DNS micro-case settings

Grid target G1: `Nx = 32`, `Ny = 33`, `Nz = 32`, channel case, `channel_height = 1.0`, `y_min = 0.0`, `y_max = 1.0`.

Time settings: `dt = 5.0e-5`, `n_steps = 300`, `final_time = 0.015`, `CFL_max_limit = 0.3`, `contact_CFL_limit = 0.2`, controlled real DNS micro-case true.

Parallel settings: `np = 1 only`; actual MPI allowed only for np=1 if required; `np=2` forbidden; `np=4` forbidden; no parallel consistency check in Stage 22.7. The presence of `MPIEXEC` is not permission to run np>1.

Fibre settings: `n_fibre = 1`, `n_point_per_fibre = 64`, `component_dim = 3`, `fibre_radius = 0.01`, `fibre_length = 0.40`, `rho_tilde = 1.0`, `gamma = 1.0e-5`, `c_fs = 1.0`.

Initial near-wall condition: one fibre near lower wall; initial minimum wall gap target = `1.0e-3`; no severe initial penetration; fail closed if initial penetration exceeds `max_penetration_allowed` or if near-wall placement cannot be guaranteed without source modification.

FSI settings: Stage 20 two-way FSI pathway enabled; `lambda_fsi = 1.0e-4`; `F_on_structure_from_fluid_candidate = -F_fs_candidate`; `F_on_fluid_from_structure_candidate = +F_fs_candidate`; action-reaction residual bounded.

Wall-contact settings: `lambda_contact = 1.0`; wall contact force enabled; `k_wall = 1.0e2`; `damping_ratio = 0.2`; `m_eff = rho_tilde * ds`; `c_wall = 2 * damping_ratio * sqrt(k_wall * m_eff)`; `max_penetration_allowed = 1.0e-4`; `max_contact_force_norm = 1.0e3`; `max_contact_acceleration = 1.0e3`.

Wall contact formula: `d_lower = y - y_min`; `d_upper = y_max - y`; `g_lower = d_lower - fibre_radius`; `g_upper = d_upper - fibre_radius`; `g_wall = min(g_lower, g_upper)`; `delta_wall = max(0, -g_wall)`; lower normal `(0,1,0)`; upper normal `(0,-1,0)`; `v_n = V_fibre · n_wall`; `v_n_minus = min(v_n,0)`; `F_wall = lambda_contact * (k_wall * delta_wall * n_wall - c_wall * v_n_minus * n_wall) if delta_wall > 0 else 0`.

Wall force rules: no force if `delta_wall = 0`; lower-wall y-force >= 0; upper-wall y-force <= 0; no attractive wall force; contact energy nonnegative; damping power nonpositive; penetration <= max_penetration_allowed.

Collision settings: `n_fibre = 1`; fibre-fibre collision force disabled; fibre collision candidate disabled; no self-pair collision; no production multifibre activation.

Restart/statistics/visualization settings: do not perform restart test; do not modify production restart/statistics/visualization schemas; normal outputs are allowed only if produced without schema changes; Stage 22.11 audits production readiness.

Build/execution policy: build only if `STAGE22_7_BUILD_ALLOWED=1`; DNS only if `STAGE22_7_PRODUCTION_DNS_ALLOWED=1` and `STAGE22_7_CONTROLLED_MICRO_CASE_ALLOWED=1`; MPI only for np=1 if `STAGE22_7_ACTUAL_MPI_ALLOWED=1`; fail closed if np is not exactly 1, n_fibre not exactly 1, lambda_contact not exactly 1.0, wall contact gate not enabled, any fibre-fibre collision gate enabled, production multifibre enabled, or grid differs from G1 without documented compatibility substitution.

## Validation groups

Gate/evidence group: Stage 22.7 requested; Stage 22.7 single-fibre near-wall contact micro-case enabled; Stage 22.6 evidence accepted; Stage 22.5 evidence accepted; Stage 22.4 evidence accepted; Stage 22.3 evidence accepted; Stage 22.2 evidence accepted; Stage 22.1 evidence accepted; Stage 22.0 evidence accepted; Stage 20 closure accepted; Stage 21 closure accepted; source-only closure acceptance preserved; missing old outputs allowed; no previous stage rerun; first real wall-contact DNS micro-case declared; cautious mode enabled.

Case configuration group: G1 grid documented; Nx/Ny/Nz valid; dt valid; n_steps valid; final_time valid; np exactly 1; np2 disabled; np4 disabled; n_fibre exactly 1; n_point_per_fibre valid; fibre parameters valid; initial near-wall placement documented; initial wall gap near target; lambda_fsi valid; lambda_contact exactly 1; wall contact gate enabled; fibre-fibre collision gate disabled.

Build/run group: build allowed only through Stage 22.7 wrapper; build directory isolated; no source modification during build; controlled DNS micro-case allowed only through wrapper; np=1 run only; np=2/4 not run; DNS run completed without crash; no NaN/Inf in log; no unknown runtime failure.

Fluid/structure/FSI groups: velocity finite; pressure finite; RHS finite; divergence bounded; CFL bounded; projection stable; Poisson stable; X/V/A finite; segment length residual bounded; structure displacement bounded; bending/tension bounded; FSI force bounded; wall contact force bounded; total structure force bounded; F_fs finite; F_on_structure finite; F_on_fluid finite; action-reaction residual bounded; Lagrangian/Eulerian FSI force bounded; FSI force conservation residual bounded; RHS delta bounded; lambda_fsi response bounded.

Wall-contact group: wall gap metadata finite; wall_gap_min reported; wall penetration max <= max_penetration_allowed; wall contact force finite and bounded; wall contact acceleration bounded; inward wall force direction valid; no attractive wall force; contact energy nonnegative; damping power nonpositive; contact_CFL bounded; wall contact only active under penetration; zero wall contact under zero penetration; wall gap bounded or improved.

Collision-disabled and isolation groups: fibre-fibre collision force disabled; fibre collision candidate disabled; n_fibre=1 collision inactive; no self-pair collision; collision force norm zero; collision energy zero; no restart test; no statistics/visualization schema modification; no Stage 10–21 files modified; no Stage 22.0–22.6 files modified; no src/CMake modification; no production DNS/RHS/IBM source modification; no projection/Poisson/RK3/channel forcing modification; no production restart/statistics/visualization schema modification; no production multifibre activation.

## Safe defaults

* STAGE22_7_ENABLE=1
* STAGE22_7_SINGLE_FIBRE_NEAR_WALL_CONTACT_MICRO_CASE_ENABLE=1
* STAGE22_7_REQUIRE_STAGE22_6_PASS=1
* STAGE22_7_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_7_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_7_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_7_DIAGNOSTIC_ONLY=0
* STAGE22_7_FAIL_CLOSED=1
* STAGE22_7_CONTROLLED_MICRO_CASE_ALLOWED=1
* STAGE22_7_BUILD_ALLOWED=1
* STAGE22_7_PRODUCTION_DNS_ALLOWED=1
* STAGE22_7_ACTUAL_MPI_ALLOWED=1
* STAGE22_7_NP=1
* STAGE22_7_NP2_ALLOWED=0
* STAGE22_7_NP4_ALLOWED=0
* STAGE22_7_NX=32
* STAGE22_7_NY=33
* STAGE22_7_NZ=32
* STAGE22_7_DT=5.0e-5
* STAGE22_7_N_STEPS=300
* STAGE22_7_FINAL_TIME=0.015
* STAGE22_7_CFL_MAX_LIMIT=0.3
* STAGE22_7_CONTACT_CFL_LIMIT=0.2
* STAGE22_7_N_FIBRE=1
* STAGE22_7_N_POINT_PER_FIBRE=64
* STAGE22_7_FIBRE_RADIUS=0.01
* STAGE22_7_FIBRE_LENGTH=0.40
* STAGE22_7_LAMBDA_FSI=1.0e-4
* STAGE22_7_LAMBDA_CONTACT=1.0
* STAGE22_7_K_WALL=1.0e2
* STAGE22_7_DAMPING_RATIO=0.2
* STAGE22_7_INITIAL_WALL_GAP_TARGET=1.0e-3
* STAGE22_7_WALL_CONTACT_FORCE_ENABLE=1
* STAGE22_7_COLLISION_FORCE_ENABLE=0
* STAGE22_7_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE=0
* STAGE22_7_COLLISION_FORCE_APPLY_ENABLE=0
* STAGE22_7_COLLISION_TO_RHS_ENABLE=0
* STAGE22_7_PRODUCTION_MULTIFIBRE_ENABLE=0
* BUILD_DIR=build_stage22_7

## Next stage

Stage 22.8: two-fibre collision micro-case.
