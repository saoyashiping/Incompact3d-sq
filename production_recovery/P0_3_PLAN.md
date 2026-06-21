# P0.3 Runtime Grid/State Bridge Plan

## Scope

P0.3 establishes a runtime bridge that maps xcompact3d RHS array shapes into fibre production containers and validates lambda=0 no-contamination through the bridge/main-hook/force-buffer call path.

## Non-goals

- No long production DNS runs.
- No full FSI closure.
- No real fibre structure advance or motion.
- No wall-contact or fibre-fibre collision forces.
- No pressure/projection/RK3/channel-forcing algorithm changes.
- No production-ready claim.

## Objectives

1. Add `fibre_prod_runtime_bridge_type` and bridge lifecycle APIs.
2. Allocate a zeroed `fibre_prod_force_buffer_type` matching RHS shape.
3. Ensure bridge initialization and reset never modify DNS RHS arrays.
4. Route lambda=0 through bridge + force buffer + main hook with strict RHS no-contamination.
5. Keep lambda>0 without physical force generation from producing pseudo RHS response in xcompact3d.
6. Preserve the P0.2 force-buffer-to-RHS validation path.

## Next stage

P0.4 must add runtime velocity interpolation bridge and one-way fibre sampling. Production-run status remains blocked.
