# Stage 19.1 production physical structure config gate

Stage 19.1 is configuration-only. It defines fail-closed feature gates for later
production-side physical-structure integration, but it does not create
production X/V/A state, structure buffers, hooks, advance APIs, commit paths,
force input paths, RHS spreading, Stage 14 RHS injection calls, restart /
statistics / visualization I/O, contact/collision forces, or production
multifibre logic.

## Required logical gates and safe defaults

| Gate | Safe default |
| --- | --- |
| `stage19_physical_structure_enable` | `false` |
| `stage19_physical_structure_diagnostic_only` | `true` |
| `stage19_physical_structure_single_fibre_only` | `true` |
| `stage19_physical_structure_fail_closed` | `true` |
| `stage19_physical_structure_state_enable` | `false` |
| `stage19_physical_structure_init_enable` | `false` |
| `stage19_physical_structure_force_candidate_enable` | `false` |
| `stage19_physical_structure_advance_candidate_enable` | `false` |
| `stage19_physical_structure_commit_enable` | `false` |
| `stage19_physical_structure_hook_enable` | `false` |
| `stage19_fluid_force_input_enable` | `false` |
| `stage19_rhs_spreading_enable` | `false` |
| `stage19_stage14_rhs_injection_enable` | `false` |
| `stage19_restart_io_enable` | `false` |
| `stage19_statistics_io_enable` | `false` |
| `stage19_visualization_io_enable` | `false` |
| `stage19_contact_model_enable` | `false` |
| `stage19_fibre_fibre_collision_enable` | `false` |
| `stage19_multifibre_production_enable` | `false` |

## Consistency rules

1. If `stage19_physical_structure_enable=false`, all runtime activation gates
   remain false: state, init, force candidate, advance candidate, commit, hook,
   fluid-force input, RHS spreading, and Stage 14 RHS injection.
2. `stage19_physical_structure_diagnostic_only=true` implies no commit, RHS
   spreading, Stage 14 RHS injection, restart I/O, statistics I/O, or
   visualization I/O activation.
3. `stage19_physical_structure_single_fibre_only=true` implies
   `stage19_multifibre_production_enable=false` and
   `stage19_fibre_fibre_collision_enable=false`.
4. `stage19_contact_model_enable=false` means no wall contact, penalty,
   repulsive, lubrication, friction, adhesion, or contact damping force.
5. `stage19_rhs_spreading_enable=false` implies no Stage 14 RHS injection, no
   fluid RHS modification, and no IBM forcing activation.
6. `stage19_physical_structure_fail_closed=true` means invalid or unknown gate
   combinations fail the helper, as does any accidental production activation.

## Evidence and false-positive policy

Stage 19.1 accepts Stage 19.0 PASS evidence when present. In source-only archives,
it safely accepts preserved Stage 19.0 source-only Stage 18 closure acceptance so
that Stage 18.0 through Stage 18.11 do not need to be rerun and individual old
Stage 18 outputs are not required. Helper-local `stage19_outputs` files are not
production restart/statistics/visualization I/O.

The gate uses pure Python and shell only; it does not require `rg`, does not run
MPI, does not run production DNS, does not build, and does not treat
documentation strings or negative-check names as production activation.

## Manual command

```bash
stage19_checks/run_stage19_1_physical_structure_config_gate.sh
```

Expected PASS evidence:

```text
STAGE 19.1 PHYSICAL STRUCTURE CONFIG GATE VERDICT: PASS
STAGE 19.1 FINAL VERDICT: PASS
final_status PASS
```
