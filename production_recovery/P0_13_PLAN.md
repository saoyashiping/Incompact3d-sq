# P0_13 Plan: no-fibre regression + P0 closure/P1 entry gate

P0_13 is the final P0 closure stage. It verifies that all production fibre paths remain default-off, no-fibre xcompact3d behavior is protected, restart/statistics/visualization sources remain compatible, and P0_2-P0_12 regression checks remain wired.

Scope:
- Source/static compatibility audits plus existing P0 micro-checks.
- No long-time DNS, no real P1 single-fibre case, and no paper-scale DNS claim.
- P1 entry can be allowed only after P0_13 PASS; paper-scale long-time DNS remains blocked.
