# Stage 22.1 full helper-chain reconstruction

Stage 22 title: final integrated validation and production-readiness closure.
Stage 22.1 title: full helper-chain reconstruction.

Stage 22.1 is helper-only. It reconstructs, in a synthetic helper environment only, the Stage 20 FSI candidate chain and the Stage 21 wall/contact/collision diagnostic metadata chain. It does not run production DNS, MPI, builds, previous stages, real contact/collision force computation, contact/collision force application, structure advance, RHS coupling, Stage 14 RHS injection, IBM, DNS-core, pressure projection, Poisson, RK3/channel forcing, production restart/statistics/visualization I/O, or production multifibre logic.

## Source-only and closed-stage policy

Stage 20 closure, Stage 21 closure, and Stage 22.0 PASS evidence are accepted from available evidence when present. Source-only archives are accepted when old closure files or outputs are absent. Missing old outputs are allowed. No previous stage is rerun. Stage 10 through Stage 21 files and Stage 22.0 files remain immutable. Existing `src`, CMake, cmake, and production DNS/RHS/IBM/I/O files remain untouched.

## Helper-only case settings

Synthetic Eulerian helper grid:

* Nx = 16
* Ny = 16
* Nz = 16
* dx = 1.0 / 16.0
* dy = 1.0 / 16.0
* dz = 1.0 / 16.0

Time integration metadata:

* dt = 1.0e-5
* n_steps = 5
* actual DNS stepping = false
* helper reconstruction only = true

Fibre settings:

* n_fibre = 2
* n_point_per_fibre = 64
* component_dim = 3
* fibre_radius = 0.01
* fibre_length = 0.40
* rho_tilde = 1.0
* gamma = 1.0e-5
* c_fs = 1.0

Lambda settings:

* lambda_fsi = 1.0e-6
* lambda_contact = 0.0

Contact/collision settings:

* contact_force_enable = false
* collision_force_enable = false
* wall_contact_force_candidate_enable = false
* fibre_collision_force_candidate_enable = false
* contact_force_apply_enable = false
* contact_in_structure_advance_enable = false
* contact_to_rhs_enable = false

## Stage 20 FSI helper reconstruction group

1. Fluid-to-structure candidate:

```text
u_relative_candidate = u_interp_candidate - V_fibre
F_fs_candidate = C_fs * u_relative_candidate
F_on_structure_from_fluid_candidate = -F_fs_candidate
```

2. Structure force candidate:

```text
F_total_without_contact_candidate = F_bending_candidate + F_tension_candidate - F_fs_candidate
```

No contact/collision force is included in Stage 22.1.

3. Structure response candidate:

```text
A_candidate = F_total_without_contact_candidate / rho_tilde
V_next_candidate = V_current + dt * A_candidate
X_next_candidate = X_current + dt * V_current + 0.5 * dt^2 * A_candidate
```

This remains helper-local only. No production structure commit is made.

4. Structure-to-fluid reaction candidate:

```text
F_on_fluid_from_structure_candidate = +F_fs_candidate
F_on_structure_from_fluid_candidate + F_on_fluid_from_structure_candidate = 0
```

5. Lagrangian-to-Eulerian force-density candidate:

```text
f_eulerian_candidate = sum_q F_on_fluid_from_structure_candidate(q) * delta_h(x_i - X_q) * ds / dV
force_conservation_residual = sum_i f_eulerian_candidate(i) * dV - sum_q F_on_fluid_from_structure_candidate(q) * ds
```

6. Lambda-gated RHS candidate:

```text
rhs_delta_candidate = lambda_fsi * f_eulerian_candidate
```

This is helper-local only. No production RHS update is made.

Validation items: u_interp_candidate finite; V_fibre finite; u_relative_candidate finite; F_fs_candidate finite; F_on_structure_from_fluid_candidate finite; F_on_fluid_from_structure_candidate finite; action_reaction_residual bounded; f_eulerian_candidate finite; force_conservation_residual bounded; rhs_delta_candidate finite and bounded; lambda_fsi applied correctly.

## Stage 21 contact/collision metadata reconstruction group

1. Wall gap metadata:

```text
d_lower = y - y_min
d_upper = y_max - y
g_lower = d_lower - fibre_radius
g_upper = d_upper - fibre_radius
g_wall = min(g_lower, g_upper)
wall_penetration_depth = max(0, -g_wall)
```

2. Fibre-fibre gap metadata:

```text
d_ff = minimum helper point/segment distance
g_ff = d_ff - 2 * fibre_radius
fibre_fibre_penetration_depth = max(0, -g_ff)
```

3. Warning/fail metadata: safe / near-contact / overlap / fail-closed classification; risk_level; risk_label; warning_trigger; fail_closed_trigger.

4. Candidate registry metadata: candidate_id; candidate_type; candidate_key; canonical_pair_key; canonical_sort_key; gap_value; penetration_depth; risk_level; risk_label; diagnostic_only; force_computation_allowed = false; force_application_allowed = false.

5. Ownership metadata: owner_rank_np1; owner_rank_np2; owner_rank_np4; local_candidate_id_np1; local_candidate_id_np2; local_candidate_id_np4; rank_candidate_count_np1; rank_candidate_count_np2; rank_candidate_count_np4.

6. Deterministic ordering metadata: global_order_index; canonical_sort_key; sorted_order_reference; repeated_eval_order_hash; ordering_deterministic.

7. Persistence metadata: schema_name; schema_version; serialization_hash; reload_hash; reconstruction_hash; roundtrip_hash; roundtrip_equal.

8. Integrated diagnostic summary: wall_gap_min; fibre_fibre_gap_min; max_risk_level; candidate_count; force_disabled_summary; production_isolation_summary.

Stage 22.1 validates that contact metadata exists but contact/collision forces remain disabled.

## Force-disabled group

* lambda_contact = 0.0
* contact force disabled
* collision force disabled
* wall contact force candidate disabled
* fibre collision force candidate disabled
* contact/collision force application disabled
* contact/collision force not added to structure total force
* contact/collision force not spread to RHS
* Stage 14 RHS injection disabled
* production RHS update disabled

## Production-isolation group

* no production DNS
* no MPI
* no build
* no src modification
* no CMake modification
* no production IBM/RHS/DNS-core/projection/Poisson/RK3/I/O modification
* no production hook activation

## Safe defaults

* STAGE22_1_ENABLE=1
* STAGE22_1_FULL_HELPER_CHAIN_RECONSTRUCTION_ENABLE=1
* STAGE22_1_REQUIRE_STAGE22_0_PASS=1
* STAGE22_1_REQUIRE_STAGE20_CLOSED=1
* STAGE22_1_REQUIRE_STAGE21_CLOSED=1
* STAGE22_1_ALLOW_SOURCE_ONLY_ARCHIVE=1
* STAGE22_1_ALLOW_MISSING_OLD_OUTPUTS=1
* STAGE22_1_DO_NOT_RERUN_PREVIOUS_STAGES=1
* STAGE22_1_DIAGNOSTIC_ONLY=1
* STAGE22_1_FAIL_CLOSED=1
* STAGE22_1_HELPER_ONLY=1
* STAGE22_1_BUILD_ALLOWED=0
* STAGE22_1_PRODUCTION_DNS_ALLOWED=0
* STAGE22_1_ACTUAL_MPI_ALLOWED=0
* STAGE22_1_NX=16
* STAGE22_1_NY=16
* STAGE22_1_NZ=16
* STAGE22_1_DT=1.0e-5
* STAGE22_1_N_STEPS=5
* STAGE22_1_N_FIBRE=2
* STAGE22_1_N_POINT_PER_FIBRE=64
* STAGE22_1_COMPONENT_DIM=3
* STAGE22_1_FIBRE_RADIUS=0.01
* STAGE22_1_FIBRE_LENGTH=0.40
* STAGE22_1_RHO_TILDE=1.0
* STAGE22_1_GAMMA=1.0e-5
* STAGE22_1_C_FS=1.0
* STAGE22_1_LAMBDA_FSI=1.0e-6
* STAGE22_1_LAMBDA_CONTACT=0.0
* STAGE22_1_CONTACT_FORCE_ENABLE=0
* STAGE22_1_COLLISION_FORCE_ENABLE=0
* STAGE22_1_WALL_CONTACT_FORCE_CANDIDATE_ENABLE=0
* STAGE22_1_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE=0
* STAGE22_1_CONTACT_FORCE_APPLY_ENABLE=0
* STAGE22_1_STRUCTURE_ADVANCE_ENABLE=0
* STAGE22_1_RHS_COUPLING_ENABLE=0
* STAGE22_1_STAGE14_RHS_INJECTION_ALLOWED=0
* STAGE22_1_PRODUCTION_RHS_UPDATE_ALLOWED=0
* STAGE22_1_PRODUCTION_RESTART_IO_ALLOWED=0
* STAGE22_1_PRODUCTION_STATISTICS_IO_ALLOWED=0
* STAGE22_1_PRODUCTION_VISUALIZATION_IO_ALLOWED=0
* STAGE22_1_PRODUCTION_MULTIFIBRE_ENABLE=0
* STAGE22_1_AUDIT_TOL=1.0e-12
* STAGE22_1_ZERO_TOL=1.0e-14

## Next stage

Stage 22.2: wall contact force candidate helper test.
