# P0.9 Plan: Controlled Structure Dry-Step Commit Gate

Goal: add a default-off commit gate that can accept finite/bounded P0.8 dry-step trial states and optionally commit only structure-side coordinates and velocity storage.

Scope:
- No pressure/projection/RK3/channel-forcing changes.
- No Eulerian force-buffer writes, spreading, RHS adapter calls, or RHS feedback.
- Commit gate remains disabled unless `FIBRE_PROD_STRUCTURE_DRY_COMMIT_ENABLE=1` is set.
- Production-run status remains STILL BLOCKED after this stage.
