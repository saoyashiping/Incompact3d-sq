# P0_12 Plan: end-to-end synthetic closed-loop + np=1/2/4 consistency

P0_12 wires the default-off synthetic production path from velocity sampling through fibre-state attachment, hydrodynamic input, structure handoff, dry-step, commit gate, reaction-force spreading, Eulerian force buffer, and lambda-gated RHS coupling.

Scope:
- Single deterministic synthetic step only.
- No long-running DNS and no production DNS-FSI readiness claim.
- No pressure/projection/RK3/channel-forcing logic changes.
- P0_12 path remains disabled unless `FIBRE_PROD_SYNTHETIC_CLOSED_LOOP_ENABLE` is explicitly enabled.

Exit criteria:
- Lambda=0 full path leaves RHS unchanged.
- Small lambda produces finite, bounded, non-uniform RHS increments from the force buffer only.
- Deterministic signatures are consistent across np=1/2/4 runs.
- P0_2 through P0_11 checks remain available and non-regressed.
