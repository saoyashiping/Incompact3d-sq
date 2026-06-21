# P0_10 Plan: Reaction-Force Candidate + Eulerian Spreading Buffer

Goal: construct a structure-to-fluid reaction-force candidate from structure-side input storage and spread it into a production Eulerian force buffer while preserving strict no-RHS-feedback behavior.

Scope:
- Use sign convention `reaction_force_candidate = -structure_input_force`.
- Allow writes only to `fibre_prod_force_buffer_type%fx/fy/fz` in this stage.
- Do not call RHS adapter or main-hook force-buffer-to-RHS APIs.
- Keep the xcompact3d path default-off via explicit reaction-spreading environment gates.
- Production-run status remains STILL BLOCKED after this stage.
