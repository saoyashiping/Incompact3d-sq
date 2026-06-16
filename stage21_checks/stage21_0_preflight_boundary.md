# Stage 21.0: wall/contact/collision preflight boundary

Stage 21.0 establishes the Stage 21 wall/contact/collision safety boundary for flexible-fibre FSI. It is a configuration, documentation, and gate audit only; it does not compute real wall distances, fibre-fibre distances, wall gaps, fibre-fibre gaps, contact forces, collision forces, contact-augmented structure advances, or RHS contact forcing.

Stage 21.0 accepts Stage 20 closure evidence when Stage 20.11 PASS evidence or `stage20_checks/STAGE20_CLOSED.md` exists. If source-only closure acceptance is enabled, missing old output folders and old closure files are allowed and no previous stage is rerun.

## Mathematical boundary for future stages

These formulas are documented for later Stage 21 substages and are not active physics in Stage 21.0:

* Wall signed gap: `g_wall = d_wall - a_f`.
* Fibre-fibre signed gap: `g_ff = d_ff - (a_i + a_j)`.
* Non-penetration requirement: `g_wall >= 0` and `g_ff >= 0`.
* Future wall contact candidate force: `F_wall_candidate = k_w * delta_wall * n_wall - c_w * v_n_minus * n_wall`.
* Future fibre-fibre collision candidate force: `F_i_candidate = k_c * delta_ij * n_ij - c_c * v_n_minus * n_ij`, with `F_j_candidate = -F_i_candidate`.
* Fibre-fibre action-reaction requirement: `F_i_candidate + F_j_candidate = 0`.

## Safe default gates

Stage 21.0 keeps wall contact, fibre collision, distance/gap audits, near-contact warnings, contact-force candidates, contact application, structure-advance contact, RHS contact, production multifibre, production DNS, and actual MPI disabled by default. Diagnostic-only and fail-closed modes are enabled.

## No production activation

The Stage 21.0 wrapper and helper do not build targets, run MPI, run production DNS, call Stage 14 RHS injection, modify IBM/DNS-core/projection/Poisson/RK3/channel-forcing paths, modify restart/statistics/visualization production I/O, or modify closed Stage 10-20 files.

## Output

```text
stage21_outputs/fibre_stage21_0_preflight_boundary.dat
```

## Manual command

```bash
stage21_checks/run_stage21_0_preflight_boundary.sh
```

## Expected PASS evidence

```text
STAGE 21.0 PREFLIGHT BOUNDARY VERDICT: PASS
STAGE 21.0 FINAL VERDICT: PASS
```

## Next stage

Stage 21.1: wall distance and signed-gap audit
