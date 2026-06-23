# P1_3 RHS formula grep fix

## Root cause

P1_3 real xcompact3d segment1 completed successfully and emitted the valid diagnostic line:

```text
P1_3 RHS increment diagnostic: nonzero finite bounded PASS formula=lambda*penalty_beta*force_buffer values=...
```

However, `P1_3_VALIDATION_COMMAND.sh` required the nonexistent literal token:

```text
RHS increment formula: PASS
```

Therefore the validation failed with:

```text
Reason: segment 1 missing RHS formula PASS
```

This was a validation-script token mismatch, not a physical DNS-FSI failure.

## Fix

The P1_3 validation script now accepts the actual emitted PASS diagnostic:

```text
P1_3 RHS increment diagnostic: ... nonzero finite bounded PASS ... formula=lambda*penalty_beta*force_buffer
```

The legacy token `RHS increment formula: PASS` is still accepted for compatibility.

## Scope

Only validation/evidence logic is changed. No xcompact3d physical solver, pressure/projection/RK3/channel forcing, fibre reaction force, force buffer, or RHS adapter physics is modified.
