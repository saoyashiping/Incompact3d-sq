# Fiber rigid coupling ramp input stabilization

This note documents how `coupling_ramp_steps` is read and how to verify it at runtime.

## Where to set `coupling_ramp_steps`

Set `coupling_ramp_steps` inside `&FiberParam ... /End`, for example:

```fortran
&FiberParam
  rigid_coupling_test_active = T
  ibm_beta = -20.0
  coupling_ramp_steps = 5
/End
```

## Multiple `&FiberParam` blocks

The reader now consumes all `&FiberParam` blocks in the file and applies them in order.
If multiple blocks are present, the last value of each field is the effective runtime value.

For clarity and reviewability, prefer a single `&FiberParam` block in production inputs.

## Runtime check

On the first coupling step, check:

```text
DEBUG_RAMP coupling_ramp_steps=...
```

Expected examples:

- input `coupling_ramp_steps = 5`  -> `DEBUG_RAMP coupling_ramp_steps=         5`
- input `coupling_ramp_steps = 10` -> `DEBUG_RAMP coupling_ramp_steps=        10`

If the first-step debug line still shows `0`, then the input value was not read successfully.
