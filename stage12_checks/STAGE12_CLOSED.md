# Stage 12 Closed

Stage 12 is closed after the Stage 12.11 total smoke gate passed.

## Purpose

Stage 12 established the production Lagrangian feedback force candidate diagnostic path. It computes and audits the Lagrangian force candidate without enabling Eulerian force density, IBM spreading, RHS injection, two-way force density, or fibre structure advance.

## Closed sub-stages

- Stage 12.0 config and global switches.
- Stage 12.1 Lagrangian force candidate buffer skeleton.
- Stage 12.2 prescribed fibre/control-point velocity skeleton.
- Stage 12.3 controlled feedback force formula checks.
- Stage 12.4 feedback sign convention audit.
- Stage 12.5 force work / power diagnostic candidate.
- Stage 12.6 production feedback candidate hook skeleton.
- Stage 12.7 np=1 force-candidate no-contamination invariance.
- Stage 12.8 np=1/2/4 parallel force-candidate consistency.
- Stage 12.9 restart / stats / visu / coarse I/O compatibility.
- Stage 12.10 RHS / spreading / structure contamination audit.
- Stage 12.11 total smoke closure.

## Governing diagnostic model

```text
U_f = I_h[u](X_f)
slip = U_f - V_f
F_fluid_to_fibre_cand = alpha * slip
F_fibre_to_fluid_cand = -F_fluid_to_fibre_cand
P_pair + P_slip = 0
f_fsi = 0
RHS_stage12 = RHS_stage11 = RHS_stage10 = RHS_stage9
```

## Closure statement

Stage 12 closes the Lagrangian feedback force candidate diagnostic path only. No Eulerian force density, no IBM spreading, no RHS injection, no two-way force density, and no fibre structure advance are activated.

Real Eulerian force-density construction begins in Stage 13.

## Next recommended stage

Stage 13 production Eulerian force-density candidate / spreading skeleton.
