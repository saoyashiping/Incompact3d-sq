# Stage 17.9 closed-loop wall/contact compatibility

Stage 17.9 verifies that the closed Stage 16 one-fibre closed-loop diagnostic chain can coexist with the closed Stage 17.3--17.8 wall/contact/collision-safety diagnostics as read-only, no-force side diagnostics.

The approved Stage 16 chain remains:

* Stage 11 fluid-to-fibre sampling
* Stage 12 feedback-force candidate
* Stage 16.4 structure-side fluid-on-fibre force input
* controlled structure-state diagnostics
* Stage 13 Eulerian force-density candidate
* Stage 14 RHS diagnostic / controlled path

The approved Stage 17 side-diagnostic chain remains:

* Stage 17.3 effective-radius wall gaps
* Stage 17.4 containment / penetration fail-closed evidence
* Stage 17.5 near-wall/contact-state classification
* Stage 17.6 segment wall-clearance safety
* Stage 17.7 no-force contact placeholder records
* Stage 17.8 standalone mock fibre-fibre placeholder geometry

Stage 17.9 is non-invasive.  It does not insert production hooks, does not modify Stage 16, Stage 13, or Stage 14 code, does not apply wall/contact/fibre-fibre force, does not alter fluid RHS, and does not activate production multi-fibre FSI.

All contact/collision response channels are required to remain zero:

```text
contact_force_norm = 0
contact_rhs_increment_norm = 0
contact_structure_update_norm = 0
fibre_fibre_force_norm = 0
fibre_fibre_rhs_increment_norm = 0
fibre_fibre_structure_update_norm = 0
```

The wrapper writes `stage17_outputs/fibre_stage17_9_closed_loop_wall_contact_compatibility.dat`, creates `stage17_outputs` if needed, and by default runs no production physics, no MPI, and no build.

The helper preserves the corrected false-positive-safe audit pattern from Stage 16 and Stage 17.0--17.8: targeted evidence checks only, no Markdown-as-code scanning, no broad repository-wide scan, no negative-check string regressions, no brittle old failure-label matching, no `rg`-only dependency, source-only archives without `.git` metadata are not treated as DNS-core contamination, and the Stage 14 small-lambda hook check uses `src/fibre_stage14_production_rhs_injection.f90` plus `src/xcompact3d.f90`.
