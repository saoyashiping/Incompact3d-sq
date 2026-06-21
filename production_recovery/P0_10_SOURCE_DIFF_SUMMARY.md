# P0_10 Source Diff Summary

- Added `fibre_prod_reaction_force_candidate` for sign-convention and diagnostics on Lagrangian reaction-force candidates.
- Added `fibre_prod_reaction_spreading_buffer` for IBM spreading into a production force buffer with finite/bounded/conservation diagnostics and no RHS path.
- Added `fibre_prod_reaction_spreading_buffer_check` covering reaction sign, buffer nonzero/finite/bounded checks, conservation scale, no-contamination, and invalid input guards.
- Added read-only structure-input access in `fibre_prod_state`.
- Added a default-off `xcompact3d` P0_10 diagnostic path and CMake target coverage.
