# Stage 22.0 preflight boundary

Stage 22 title: final integrated validation and production-readiness closure.
Stage 22.0 title: final integrated validation and production-readiness preflight.

Stage 22 is the final large stage. No Stage 23 is planned. Stage 22.12 will produce final closure evidence if all previous Stage 22 sub-stages pass.

## Source-only and closed-stage policy

Stage 22.0 is diagnostic-only and fail-closed. It does not run DNS, MPI, builds, previous stages, new physics, wall contact force, fibre-fibre collision force, contact/collision application, structure advance, RHS coupling, Stage 14 RHS injection, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, production restart/statistics/visualization I/O, or production multifibre logic.

Stage 20 closure and Stage 21 closure are accepted from available closure evidence when present, and from the user-reported closed state for source-only archives when closure files or old outputs are absent. Missing old outputs are allowed. Previous stages are not rerun. Closed-stage immutability is required: Stage 10 through Stage 21 files, `src`, CMake, cmake, and production DNS/RHS/IBM/I/O paths must not be modified by Stage 22.0.

## Stage 22 sub-stage plan

* Stage 22.0: final integrated validation and production-readiness preflight
* Stage 22.1: full helper-chain reconstruction
* Stage 22.2: wall contact force candidate helper test
* Stage 22.3: fibre-fibre collision force candidate helper test
* Stage 22.4: contact force into structure candidate
* Stage 22.5: lambda/no-contact/contact regression
* Stage 22.6: single-fibre channel FSI micro-case
* Stage 22.7: single-fibre near-wall contact micro-case
* Stage 22.8: two-fibre collision micro-case
* Stage 22.9: mesh/time-step sensitivity check
* Stage 22.10: np=1/2/4 parallel consistency
* Stage 22.11: restart/statistics/visualization production-readiness audit
* Stage 22.12: final total closure

## Conservative global mesh ladder

* G0 helper-only: synthetic Eulerian helper grid 16 x 16 x 16; no production DNS.
* G1 coarse DNS micro-case: Nx = 32, Ny = 33, Nz = 32, default dt = 1.0e-4.
* G2 medium DNS micro-case: Nx = 48, Ny = 49, Nz = 48, default dt = 5.0e-5.
* G3 optional refinement check: Nx = 64, Ny = 65, Nz = 64, default dt = 2.5e-5.

G1/G2/G3 are later-stage targets only. Stage 22.0 does not run any of these cases.

## Global time-step safety limits

* CFL_max_limit = 0.3
* contact_CFL_limit = 0.2
* max_structure_step_displacement_fraction = 0.1
* helper_dt = 1.0e-5
* G1_dt = 1.0e-4
* G2_dt = 5.0e-5
* G3_dt = 2.5e-5

## Fibre parameter ladder

* n_fibre_default = 1
* n_fibre_collision_test = 2
* n_point_default = 64
* n_point_refined = 128
* fibre_radius = 0.01
* fibre_length = 0.40
* rho_tilde = 1.0
* gamma = 1.0e-5
* c_fs = 1.0

## Contact parameter ladder

* k_wall_default = 1.0e2
* k_collision_default = 1.0e2
* damping_ratio = 0.2
* max_penetration_allowed = 1.0e-4
* max_contact_force_norm = 1.0e3
* max_contact_acceleration = 1.0e3

## Lambda ladder

* lambda_fsi_zero = 0.0
* lambda_fsi_small = 1.0e-6
* lambda_fsi_prod_test = 1.0e-4
* lambda_contact_zero = 0.0
* lambda_contact_small = 1.0e-6
* lambda_contact_prod_test = 1.0

## Final integrated validation metrics

Fluid metrics: velocity finite; pressure finite; divergence bounded; CFL bounded; RHS delta bounded.

Structure metrics: X finite; V finite; A finite; segment length residual bounded; bending/tension residual bounded; structure displacement per step bounded.

FSI metrics: action-reaction residual bounded; Lagrangian total force bounded; Eulerian force integral bounded; force conservation residual bounded; lambda=0 no-op; small-lambda bounded response.

Contact/collision metrics: wall_gap_min bounded; fibre_fibre_gap_min bounded; wall_penetration_max <= max_penetration_allowed; fibre_fibre_penetration_max <= max_penetration_allowed; contact force bounded; contact energy nonnegative; damping power nonpositive; no attractive wall force; fibre-fibre action-reaction residual bounded.

Parallel/I-O metrics: np=1/2/4 consistency; restart consistency; statistics finite; visualization finite; no production I/O schema contamination.

## Conservative Stage 22.0 default gates

Disabled: production_dns_execution; actual_mpi_execution; build_execution; contact_force_computation; collision_force_computation; contact_force_apply; collision_force_apply; contact_in_structure_advance; collision_in_structure_advance; contact_to_rhs; collision_to_rhs; stage14_rhs_injection; production_rhs_update; production_restart_io; production_statistics_io; production_visualization_io; production_multifibre_activation.

Enabled: diagnostic_only; fail_closed; do_not_rerun_previous_stages; source_only_closure_acceptance; closed_stage_immutability_check; production_contamination_check; no_stage23_declaration_check.

## Safe defaults

* STAGE22_0_ENABLE=1
* STAGE22_0_PREFLIGHT_ENABLE=1
* STAGE22_0_REQUIRE_STAGE20_CLOSED=1
* STAGE22_0_REQUIRE_STAGE21_CLOSED=1
* STAGE22_0_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_0_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_0_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_0_DIAGNOSTIC_ONLY=1
* STAGE22_0_FAIL_CLOSED=1
* STAGE22_0_STAGE22_IS_FINAL_STAGE=1
* STAGE22_0_NO_STAGE23_PLANNED=1
* STAGE22_0_BUILD_ALLOWED=0
* STAGE22_0_PRODUCTION_DNS_ALLOWED=0
* STAGE22_0_ACTUAL_MPI_ALLOWED=0
* STAGE22_0_CONTACT_FORCE_ENABLE=0
* STAGE22_0_COLLISION_FORCE_ENABLE=0
* STAGE22_0_CONTACT_FORCE_APPLY_ENABLE=0
* STAGE22_0_STRUCTURE_ADVANCE_ENABLE=0
* STAGE22_0_RHS_COUPLING_ENABLE=0
* STAGE22_0_STAGE14_RHS_INJECTION_ALLOWED=0
* STAGE22_0_PRODUCTION_RHS_UPDATE_ALLOWED=0
* STAGE22_0_PRODUCTION_RESTART_IO_ALLOWED=0
* STAGE22_0_PRODUCTION_STATISTICS_IO_ALLOWED=0
* STAGE22_0_PRODUCTION_VISUALIZATION_IO_ALLOWED=0
* STAGE22_0_PRODUCTION_MULTIFIBRE_ENABLE=0

## Next stage

Stage 22.1: full helper-chain reconstruction.
