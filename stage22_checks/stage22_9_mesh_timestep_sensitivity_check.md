# Stage 22.9 mesh/time-step sensitivity check

Stage 22.9 is a cautious sensitivity screen, not a full convergence campaign.  It compares G1 and G2 controlled single-rank summaries only, keeps np=1, does not run np=2 or np=4, does not run G3 by default, does not perform restart testing, and does not modify production source, CMake, DNS/RHS/IBM, pressure projection, Poisson, RK3/channel forcing, restart, statistics, visualization, or closed-stage files.

## Source-only and prior-stage acceptance

Stage 22.9 accepts Stage 20 closure, Stage 21 closure, and Stage 22.0 through Stage 22.8 PASS evidence through available committed artifacts or source-only closure acceptance.  Missing old outputs are allowed and previous stages are not rerun.

## Default sensitivity cases

### Case G1_BASE

* grid label: G1
* Nx = 32
* Ny = 33
* Nz = 32
* dt = 1.0e-4
* n_steps = 200
* final_time = 0.02
* np = 1
* n_fibre = 1
* n_point_per_fibre = 64
* lambda_fsi = 1.0e-4
* lambda_contact = 0.0
* contact/collision disabled for baseline sensitivity

### Case G2_MEDIUM

* grid label: G2
* Nx = 48
* Ny = 49
* Nz = 48
* dt = 5.0e-5
* n_steps = 300
* final_time = 0.015
* np = 1
* n_fibre = 1
* n_point_per_fibre = 64
* lambda_fsi = 1.0e-4
* lambda_contact = 0.0
* contact/collision disabled for baseline sensitivity

### Case G3_OPTIONAL

* grid label: G3
* Nx = 64
* Ny = 65
* Nz = 64
* dt = 2.5e-5
* n_steps = 400
* n_point_per_fibre = 128
* run_allowed_by_default = false
* requires explicit future enable flag: STAGE22_9_G3_OPTIONAL_ENABLE=1

G3 optional documented but not run by default.  Stage 22.9 must not run G3 unless a future instruction explicitly permits it and the optional enable flag is set.

## Metrics and acceptance envelope

The audit compares finite fluid, structure, FSI, gap, force, energy, damping, and boundedness summaries.  Required sensitivity metrics include velocity norm ratio, structure displacement ratio, FSI force norm ratio, RHS delta norm ratio, wall gap difference, segment length residual ratio, and force conservation residual ratio.

Broad safety envelope:

* CFL_max <= 0.3
* structure_step_displacement_fraction <= 0.1
* segment_length_residual_max <= 1.0e-6
* force_conservation_residual <= 1.0e-8
* action_reaction_residual <= 1.0e-10
* all G1-to-G2 ratios finite
* G1-to-G2 structure displacement ratio <= 10
* G1-to-G2 force norm ratio <= 10
* G1-to-G2 RHS delta norm ratio <= 10
* wall penetration max <= 1.0e-4
* no strict convergence claim made

## Execution policy

* Build is allowed only through the Stage 22.9 wrapper.
* DNS run is allowed only through the controlled Stage 22.9 sensitivity wrapper.
* MPI execution is np=1 only.
* np=2 and np=4 are forbidden.
* restart testing is forbidden.
* source/schema modification attempts must fail closed.
* uncontrolled production multifibre is forbidden.

## Safe defaults

* STAGE22_9_ENABLE=1
* STAGE22_9_MESH_TIMESTEP_SENSITIVITY_CHECK_ENABLE=1
* STAGE22_9_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_9_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_9_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_9_FAIL_CLOSED=1
* STAGE22_9_CONTROLLED_SENSITIVITY_CHECK_ALLOWED=1
* STAGE22_9_BUILD_ALLOWED=1
* STAGE22_9_PRODUCTION_DNS_ALLOWED=1
* STAGE22_9_ACTUAL_MPI_ALLOWED=1
* STAGE22_9_NP=1
* STAGE22_9_NP2_ALLOWED=0
* STAGE22_9_NP4_ALLOWED=0
* STAGE22_9_G1_ENABLE=1
* STAGE22_9_G2_ENABLE=1
* STAGE22_9_G3_OPTIONAL_ENABLE=0
* STAGE22_9_N_FIBRE=1
* STAGE22_9_N_POINT_PER_FIBRE=64
* STAGE22_9_LAMBDA_FSI=1.0e-4
* STAGE22_9_LAMBDA_CONTACT=0.0
* STAGE22_9_CFL_MAX_LIMIT=0.3
* STAGE22_9_RATIO_SAFETY_LIMIT=10.0
* STAGE22_9_PRODUCTION_RESTART_TEST_ALLOWED=0
* STAGE22_9_PRODUCTION_RESTART_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_9_PRODUCTION_STATISTICS_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_9_PRODUCTION_VISUALIZATION_SCHEMA_MODIFICATION_ALLOWED=0
* STAGE22_9_UNCONTROLLED_PRODUCTION_MULTIFIBRE_ALLOWED=0

## Required PASS evidence

The audit writes `stage22_outputs/fibre_stage22_9_mesh_timestep_sensitivity_check.dat` and must contain:

```text
STAGE 22.9 MESH TIMESTEP SENSITIVITY CHECK VERDICT: PASS
STAGE 22.9 FINAL VERDICT: PASS
```

## Next stage

Stage 22.10: np=1/2/4 parallel consistency.
