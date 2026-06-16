# Stage 21.3: near-contact warning and fail-closed gate

Stage 21.3 is a diagnostic-only near-contact warning and fail-closed gate for helper-local fibre-fibre signed gaps. It classifies controlled signed-gap cases as `SAFE`, `NEAR_CONTACT`, `OVERLAP`, or `FAIL_CLOSED` and records risk levels without computing or applying contact/collision forces.

Stage 21.3 accepts Stage 21.2 PASS evidence when present and preserves Stage 20/Stage 21 source-only closure acceptance. Missing old outputs and closure files are allowed when source-only acceptance is enabled, and no previous stage is rerun.

## Classification mathematics

Given signed fibre-fibre gap `g_ff`, warning threshold `g_warn`, and fail threshold `delta_fail`:

* `SAFE` when `g_ff > g_warn`, with risk level `0`.
* `NEAR_CONTACT` when `0 <= g_ff <= g_warn`, with risk level `1`.
* `OVERLAP` when `g_ff < 0` and `max(0, -g_ff) <= delta_fail`, with risk level `2`.
* `FAIL_CLOSED` when `max(0, -g_ff) > delta_fail`, with risk level `3`.
* Penetration depth is `delta_p = max(0, -g_ff)`.
* Warning trigger is active for near-contact, overlap, and fail-closed states.
* Fail-closed trigger is active only when `delta_p > delta_fail`.

## Controlled diagnostic cases

* Case A: `g_ff = 0.08`, expected `SAFE`.
* Case B: `g_ff = 0.02`, expected `NEAR_CONTACT`.
* Case C: `g_ff = -0.00002`, expected `OVERLAP`.
* Case D: `g_ff = -0.01`, expected `FAIL_CLOSED`.

The default warning gap is `0.05`, and the default penetration fail limit is `1.0e-4`.

## No production activation

Stage 21.3 does not compute collision force, contact force, wall force, structure advance updates, RHS updates, Stage 14 RHS injection, MPI, production DNS, production multifibre logic, or production source/CMake changes.

## Output

```text
stage21_outputs/fibre_stage21_3_near_contact_warning_gate.dat
```

## Manual command

```bash
stage21_checks/run_stage21_3_near_contact_warning_gate.sh
```

## Expected PASS evidence

```text
STAGE 21.3 NEAR-CONTACT WARNING GATE VERDICT: PASS
STAGE 21.3 FINAL VERDICT: PASS
```

## Next stage

Stage 21.4: contact candidate registry
