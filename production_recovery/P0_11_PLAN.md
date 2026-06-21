# P0_11 Plan: Force-Buffer to RHS Lambda-Gated Synthetic Coupling

Goal: pass a nonzero production force buffer through the existing lambda-gated main-hook/RHS-adapter path and verify strict lambda=0 no-contamination plus bounded, linearly scaled small-lambda increments.

Scope:
- Reuse the P0_10 force-buffer path.
- Default-off xcompact3d gate via `FIBRE_PROD_FORCE_BUFFER_RHS_GATE_ENABLE`.
- No pressure/projection/RK3/channel-forcing changes.
- No long production DNS and no production-ready claim.
