# Stage 22.10 np=1/2/4 parallel consistency

Stage 22.10 is a bounded np=1/2/4 parallel consistency screen for the controlled G1 micro-case family.  It compares np=1, np=2, and np=4 summaries only; does not run G2 or G3; does not perform restart testing; does not perform mesh refinement; introduces no new physics; and does not modify production source, CMake, DNS/RHS/IBM, pressure projection, Poisson, RK3/channel forcing, restart, statistics, visualization, or closed-stage files.

## Source-only and prior-stage acceptance

Stage 22.10 accepts Stage 20 closure, Stage 21 closure, and Stage 22.0 through Stage 22.9 PASS evidence through available committed artifacts or source-only closure acceptance.  Missing old outputs are allowed and previous stages are not rerun.

## Default G1 parallel-consistency plan

* grid label: G1
* Nx = 32
* Ny = 33
* Nz = 32
* channel case
* channel_height = 1.0
* y_min = 0.0
* y_max = 1.0
* dt = 1.0e-4
* n_steps = 100
* final_time = 0.01
* CFL_max_limit = 0.3
* controlled parallel consistency check = true
* np values = 1, 2, 4
* np=1 reference
* np=2 comparison
* np=4 comparison
* all runs use identical physical case parameters
* isolated outputs per np
* no restart handoff between np values

## Case groups

### single_fibre_fsi_no_contact

* n_fibre = 1
* n_point_per_fibre = 64
* fibre_radius = 0.01
* fibre_length = 0.40
* rho_tilde = 1.0
* gamma = 1.0e-5
* c_fs = 1.0
* lambda_fsi = 1.0e-4
* lambda_contact = 0.0
* wall contact disabled
* fibre-fibre collision disabled
* purpose: verify baseline FSI parallel consistency

### single_fibre_near_wall_contact

* n_fibre = 1
* n_point_per_fibre = 64
* lambda_fsi = 1.0e-4
* lambda_contact = 1.0
* wall contact enabled
* fibre-fibre collision disabled
* purpose: verify wall-contact metadata/force parallel consistency
* helper-level contact ownership/order consistency fallback is documented if a full wall-contact DNS comparison is not safely supported

### two_fibre_collision_metadata_np124

* n_fibre = 2
* n_point_per_fibre = 64
* lambda_contact = 1.0
* fibre-fibre collision metadata enabled
* full two-fibre DNS collision rerun is not required by default
* candidate pair ownership/order/determinism is checked across np=1,2,4

## Acceptance rules

* np=1, np=2, and np=4 summaries must complete without crash.
* unsupported np values are forbidden.
* G2 disabled.
* G3 disabled.
* restart test disabled.
* mesh refinement disabled.
* no new physics introduced.
* physical case parameters identical across np values.
* all metrics finite and bounded.
* CFL_max <= 0.3 for all np values.
* force conservation residual <= 1.0e-8.
* action-reaction residual <= 1.0e-10.
* decomposition differences must be within explicit tolerances.
* candidate order deterministic across np.
* pair ownership deterministic across np.
* no duplicate pair across np.
* no self-pair across np.
* no worse instability at np=2 or np=4 relative to np=1.

## Safe defaults

* STAGE22_10_ENABLE=1
* STAGE22_10_NP124_PARALLEL_CONSISTENCY_ENABLE=1
* STAGE22_10_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_10_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_10_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_10_FAIL_CLOSED=1
* STAGE22_10_CONTROLLED_PARALLEL_CHECK_ALLOWED=1
* STAGE22_10_BUILD_ALLOWED=1
* STAGE22_10_PRODUCTION_DNS_ALLOWED=1
* STAGE22_10_ACTUAL_MPI_ALLOWED=1
* STAGE22_10_NP_VALUES=1,2,4
* STAGE22_10_ALLOWED_NP_VALUES=1,2,4
* STAGE22_10_G1_ENABLE=1
* STAGE22_10_G2_ENABLE=0
* STAGE22_10_G3_ENABLE=0
* STAGE22_10_NX=32
* STAGE22_10_NY=33
* STAGE22_10_NZ=32
* STAGE22_10_DT=1.0e-4
* STAGE22_10_N_STEPS=100
* STAGE22_10_FINAL_TIME=0.01
* STAGE22_10_CFL_MAX_LIMIT=0.3
* STAGE22_10_CASE_A_SINGLE_FIBRE_FSI_ENABLE=1
* STAGE22_10_CASE_B_SINGLE_FIBRE_WALL_CONTACT_ENABLE=1
* STAGE22_10_CASE_C_TWO_FIBRE_METADATA_ENABLE=1
* STAGE22_10_FULL_TWO_FIBRE_COLLISION_DNS_REQUIRED=0
* STAGE22_10_LAMBDA_FSI=1.0e-4
* STAGE22_10_LAMBDA_CONTACT_OFF=0.0
* STAGE22_10_LAMBDA_CONTACT_ON=1.0
* STAGE22_10_PRODUCTION_RESTART_TEST_ALLOWED=0
* STAGE22_10_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_10_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_10_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_10_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED=0
* STAGE22_10_PARALLEL_SIGNATURE_ABS_TOL=1.0e-8
* STAGE22_10_PARALLEL_SIGNATURE_REL_TOL=1.0e-6

## Required PASS evidence

The audit writes `stage22_outputs/fibre_stage22_10_np124_parallel_consistency.dat` and must contain:

```text
STAGE 22.10 NP124 PARALLEL CONSISTENCY VERDICT: PASS
STAGE 22.10 FINAL VERDICT: PASS
```

## Next stage

Stage 22.11: restart/statistics/visualization production-readiness audit.
