# Stage 19.10: controlled restart/statistics/visualization state audit

Stage 19.10 is a helper-local audit of candidate restart, statistics, and visualization
state for the Stage 19 controlled one-fibre structure path. It rebuilds the deterministic
Stage 19.8-style single-fibre response while preserving Stage 19.9 lambda=0/no-fluid-coupling
semantics, then checks candidate schema completeness, shapes, finite values, ownership
metadata, and marker invariance.

## Scope

This stage is audit-only. It may construct helper-local candidate schemas for future restart,
statistics, and visualization compatibility, but it must not write production restart,
statistics, or visualization output and must not modify production I/O source. It must not
insert a runtime hook, activate a production advance or commit path, call Stage 14 RHS
injection, spread force to the Eulerian RHS, modify IBM/DNS-core/projection/Poisson/RK3, run
MPI, or run production DNS.

## Required candidate schemas

The restart candidate carries final-step structure state and metadata: `n_fibre`, `n_point`,
`component_dim`, `fibre_length`, `ds`, `dt`, `step_index`, `X_current`, `V_current`,
`A_current`, force candidates, `lambda_coupling`, effective coupling norm, ownership IDs,
`state_valid`, and `container_initialized`.

The statistics candidate carries per-step displacement, velocity, acceleration, force, and
effective-coupling norms plus finite, bounded, and no-fluid-coupling flags. The visualization
candidate carries vector fields (`X_current`, `V_current`, `A_current`, `F_total_candidate`),
scalar magnitudes, and point IDs.

## Safe defaults

The wrapper defaults to a single fibre with `n_point=64`, `component_dim=3`,
`fibre_length=1.0`, `dt=1.0e-5`, `n_steps=5`, `lambda_coupling=0.0`, `rho_l=1.0`,
`rho_tilde=1.0`, `bending_stiffness=1.0e-5`, `gamma=1.0e-5`, sine initialization, constant
zero tension, and zero controlled helper force. Audit-only gates are enabled, while production
restart/statistics/visualization I/O, fluid-force input, RHS spreading, Stage 14 injection,
and runtime hooks are disabled.

## Validation policy

The helper accepts preserved Stage 19.0-19.9 evidence by existing PASS artifacts or source-only
helper/wrapper/documentation files, and accepts Stage 18 closure source-only evidence without
requiring Stage 18.0-18.11 reruns. It checks only helper-local candidates and marker arrays;
helper-local diagnostic output under `stage19_outputs` is not production I/O contamination.

## Manual command

```bash
stage19_checks/run_stage19_10_restart_stats_visu_state_audit.sh
```

Expected PASS evidence is written to:

```text
stage19_outputs/fibre_stage19_10_restart_stats_visu_state_audit.dat
```

The wrapper prints:

```text
STAGE 19.10 RESTART STATS VISU STATE AUDIT VERDICT: PASS
STAGE 19.10 FINAL VERDICT: PASS
```
