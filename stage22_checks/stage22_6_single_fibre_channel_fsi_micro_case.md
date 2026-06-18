# Stage 22.6 single-fibre channel FSI micro-case

Stage 22.6 is the first controlled real DNS micro-case in Stage 22: one flexible fibre in channel flow with Stage 20 two-way FSI enabled, contact/collision disabled, `lambda_contact=0`, and `np=1` only. Controlled build and controlled np=1 DNS execution are allowed only through the new Stage 22.6 wrapper. Stage 22.6 must not run np=2 or np=4, enable contact or collision force, apply contact/collision to structure advance, spread contact/collision to RHS, introduce new RHS hooks, modify production source, modify DNS-core/IBM/projection/Poisson/RK3/channel forcing/restart/statistics/visualization schemas, activate production multifibre logic, perform restart testing, perform refinement, or perform np=1/2/4 consistency checks.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, and Stage 22.0 through Stage 22.5 PASS evidence are accepted through available evidence or source-only closure acceptance. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 through Stage 22.5 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/restart/statistics/visualization schema files remain untouched.

## Controlled DNS micro-case settings

Grid target G1: `Nx = 32`, `Ny = 33`, `Nz = 32`, channel case, `channel_height = 1.0`, `y_min = 0.0`, `y_max = 1.0`.

Time settings: `dt = 1.0e-4`, `n_steps = 200`, `final_time = 0.02`, `CFL_max_limit = 0.3`, `helper-only = false`, `controlled real DNS micro-case = true`.

Parallel settings: `np = 1 only`; actual MPI is allowed only for np=1 if required by the executable; `np=2` forbidden; `np=4` forbidden; no parallel consistency check in Stage 22.6. The presence of `MPIEXEC` is not permission to run np>1.

Fibre settings: `n_fibre = 1`, `n_point_per_fibre = 64`, `component_dim = 3`, `fibre_radius = 0.01`, `fibre_length = 0.40`, `rho_tilde = 1.0`, `gamma = 1.0e-5`, `c_fs = 1.0`.

FSI settings: Stage 20 two-way FSI pathway enabled; `lambda_fsi = 1.0e-4`; `F_on_structure_from_fluid_candidate = -F_fs_candidate`; `F_on_fluid_from_structure_candidate = +F_fs_candidate`; action-reaction residual must be bounded.

Contact/collision settings: `lambda_contact = 0.0`; wall contact force disabled; fibre-fibre collision force disabled; wall contact candidate disabled; fibre collision candidate disabled; contact/collision force application disabled; contact/collision not added to F_total; contact/collision not spread to RHS; contact/collision energy zero; contact/collision force norm zero.

Restart/statistics/visualization settings: do not perform restart test; do not modify production restart schema; do not modify production statistics schema; do not modify production visualization schema; normal existing outputs are allowed only if produced by the controlled micro-case without schema changes; Stage 22.11 audits production readiness.

Build/execution policy: build is allowed only if `STAGE22_6_BUILD_ALLOWED=1`; DNS run is allowed only if `STAGE22_6_PRODUCTION_DNS_ALLOWED=1` and `STAGE22_6_CONTROLLED_MICRO_CASE_ALLOWED=1`; MPI execution is allowed only for np=1 and only if `STAGE22_6_ACTUAL_MPI_ALLOWED=1`; wrapper must fail closed if np is not exactly 1, `lambda_contact` is not exactly 0.0, `n_fibre` is not exactly 1, any contact/collision gate is nonzero, or grid differs from G1 except documented compatibility substitution.

Recommended environment compatibility: `DECOMP2D_ROOT=/home/sq/opt/2decomp-fft-xcompact3d-v2.0.4`, `BUILD_DIR=build_stage22_6`, `MPIEXEC=mpirun`, `MPIEXEC_FLAGS=--mca btl self,vader,tcp`.

## Validation groups

Gate/evidence group: Stage 22.6 requested; Stage 22.6 single-fibre channel FSI micro-case enabled; Stage 22.5 evidence accepted; Stage 22.4 evidence accepted; Stage 22.3 evidence accepted; Stage 22.2 evidence accepted; Stage 22.1 evidence accepted; Stage 22.0 evidence accepted; Stage 20 closure accepted; Stage 21 closure accepted; source-only closure acceptance preserved; missing old outputs allowed; no previous stage rerun; first real DNS micro-case declared; cautious mode enabled.

Case configuration group: G1 grid documented; Nx/Ny/Nz valid; dt valid; n_steps valid; final_time valid; np exactly 1; np2 disabled; np4 disabled; n_fibre exactly 1; n_point_per_fibre valid; fibre parameters valid; lambda_fsi valid; lambda_contact zero; contact/collision gates disabled.

Build/run group: build allowed only through Stage 22.6 wrapper; build directory isolated; no source modification during build; controlled DNS micro-case allowed only through Stage 22.6 wrapper; np=1 run only; np=2/4 not run; DNS run completed without crash; no NaN/Inf in log; no unknown runtime failure.

Fluid/structure/FSI groups: velocity finite; pressure finite; RHS finite; divergence bounded; CFL bounded; projection stable; Poisson stable; X/V/A finite; segment length residual bounded; structure displacement bounded; bending/tension bounded; F_fs finite; F_on_structure finite; F_on_fluid finite; action-reaction residual bounded; Lagrangian total force bounded; Eulerian force integral bounded; force conservation residual bounded; RHS delta bounded; lambda_fsi response bounded.

Contact/collision disabled and gap groups: lambda_contact zero; wall contact force disabled; fibre-fibre collision force disabled; contact/collision force norm zero; contact/collision energy zero; contact/collision not added to structure total force; contact/collision not spread to RHS; Stage 14 contact/collision RHS injection disabled; no wall contact activation; no fibre collision activation; wall gap metadata finite; wall_gap_min reported; wall penetration max <= max_penetration_allowed; n_fibre=1 so fibre-fibre diagnostics are not active as collision; no self-pair contamination.

Production isolation group: no Stage 10–21 files modified; no Stage 22.0–22.5 files modified; no src modification; no CMake modification; no production DNS/RHS/IBM source modification; no pressure projection modification; no Poisson modification; no RK3/channel forcing modification; no restart/statistics/visualization schema modification; no production multifibre activation; no production hook activation beyond already closed Stage 20 path.

## Safe defaults

* STAGE22_6_ENABLE=1
* STAGE22_6_SINGLE_FIBRE_CHANNEL_FSI_MICRO_CASE_ENABLE=1
* STAGE22_6_REQUIRE_STAGE22_5_PASS=1
* STAGE22_6_REQUIRE_STAGE22_4_PASS=1
* STAGE22_6_REQUIRE_STAGE22_3_PASS=1
* STAGE22_6_REQUIRE_STAGE22_2_PASS=1
* STAGE22_6_REQUIRE_STAGE22_1_PASS=1
* STAGE22_6_REQUIRE_STAGE22_0_PASS=1
* STAGE22_6_REQUIRE_STAGE20_CLOSED=1
* STAGE22_6_REQUIRE_STAGE21_CLOSED=1
* STAGE22_6_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_6_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_6_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_6_DIAGNOSTIC_ONLY=0
* STAGE22_6_FAIL_CLOSED=1
* STAGE22_6_CONTROLLED_MICRO_CASE_ALLOWED=1
* STAGE22_6_BUILD_ALLOWED=1
* STAGE22_6_PRODUCTION_DNS_ALLOWED=1
* STAGE22_6_ACTUAL_MPI_ALLOWED=1
* STAGE22_6_NP=1
* STAGE22_6_NP2_ALLOWED=0
* STAGE22_6_NP4_ALLOWED=0
* STAGE22_6_NX=32
* STAGE22_6_NY=33
* STAGE22_6_NZ=32
* STAGE22_6_DT=1.0e-4
* STAGE22_6_N_STEPS=200
* STAGE22_6_FINAL_TIME=0.02
* STAGE22_6_CFL_MAX_LIMIT=0.3
* STAGE22_6_N_FIBRE=1
* STAGE22_6_N_POINT_PER_FIBRE=64
* STAGE22_6_COMPONENT_DIM=3
* STAGE22_6_FIBRE_RADIUS=0.01
* STAGE22_6_FIBRE_LENGTH=0.40
* STAGE22_6_RHO_TILDE=1.0
* STAGE22_6_GAMMA=1.0e-5
* STAGE22_6_C_FS=1.0
* STAGE22_6_LAMBDA_FSI=1.0e-4
* STAGE22_6_LAMBDA_CONTACT=0.0
* STAGE22_6_WALL_CONTACT_FORCE_ENABLE=0
* STAGE22_6_COLLISION_FORCE_ENABLE=0
* STAGE22_6_WALL_CONTACT_FORCE_CANDIDATE_ENABLE=0
* STAGE22_6_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE=0
* STAGE22_6_CONTACT_FORCE_APPLY_ENABLE=0
* STAGE22_6_COLLISION_FORCE_APPLY_ENABLE=0
* STAGE22_6_CONTACT_TO_RHS_ENABLE=0
* STAGE22_6_COLLISION_TO_RHS_ENABLE=0
* STAGE22_6_STAGE14_CONTACT_COLLISION_RHS_INJECTION_ALLOWED=0
* STAGE22_6_PRODUCTION_RESTART_TEST_ALLOWED=0
* STAGE22_6_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_6_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_6_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_6_PRODUCTION_MULTIFIBRE_ENABLE=0

## Next stage

Stage 22.7: single-fibre near-wall contact micro-case.
