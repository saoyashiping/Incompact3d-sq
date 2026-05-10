#!/usr/bin/env python3
from pathlib import Path

vals = {}
for line in Path('stage6_outputs/fibre_stage6_layout_guard_check.dat').read_text().splitlines():
    parts = line.split()
    if len(parts) >= 2:
        vals[parts[0]] = float(parts[1])

def need(k):
    if k not in vals:
        raise AssertionError(f'missing key: {k}')
    return vals[k]

def is_zero(k, tol=1e-14):
    return abs(need(k)) <= tol

assert need('stage6_layout_uniform_collocated_supported_flag') == 1
assert need('stage6_layout_uniform_collocated_blocked_flag') == 0
assert need('stage6_layout_uniform_ordinary_path_allowed_flag') == 1
assert need('stage6_layout_uniform_guard_status') == 1

assert need('stage6_layout_nonuniform_y_detected_flag') == 1
assert need('stage6_layout_nonuniform_y_blocked_flag') == 1
assert need('stage6_layout_nonuniform_y_ordinary_path_allowed_flag') == 0
assert need('stage6_layout_nonuniform_y_block_reason_code') > 0

assert need('stage6_layout_staggered_detected_flag') == 1
assert need('stage6_layout_staggered_blocked_flag') == 1
assert need('stage6_layout_staggered_ordinary_path_allowed_flag') == 0
assert need('stage6_layout_staggered_block_reason_code') > 0

assert need('stage6_layout_unknown_detected_flag') == 1
assert need('stage6_layout_unknown_blocked_flag') == 1
assert need('stage6_layout_unknown_ordinary_path_allowed_flag') == 0
assert need('stage6_layout_unknown_block_reason_code') > 0

assert need('stage6_layout_blocked_interpolation_called_count') == 0
assert need('stage6_layout_blocked_spreading_called_count') == 0
assert need('stage6_layout_blocked_rhs_injection_called_count') == 0
assert need('stage6_layout_blocked_fluid_update_called_count') == 0

assert need('stage6_layout_blocked_nonzero_buffer_max_abs') > 0
assert is_zero('stage6_layout_blocked_rhs_change_max')
assert need('stage6_layout_blocked_injected_flag') == 0
assert need('stage6_layout_blocked_modified_flag') == 0

assert need('stage6_layout_uniform_controlled_guard_pass_flag') == 1
assert need('stage6_layout_uniform_controlled_rhs_allowed_flag') == 1
assert need('stage6_layout_uniform_controlled_layout_status') == 1

assert is_zero('stage6_layout_default_supported_rhs_change_max')
assert need('stage6_layout_default_supported_injected_flag') == 0
assert need('stage6_layout_default_supported_modified_flag') == 0

assert need('stage6_layout_invalid_rejected_flag') == 1
assert is_zero('stage6_layout_invalid_rhs_change_max')
assert need('stage6_layout_production_enabled_rejected_flag') == 1
assert is_zero('stage6_layout_production_enabled_rhs_change_max')

assert is_zero('stage6_layout_velocity_change_max')
assert need('stage6_layout_fluid_update_called_flag') == 0
assert need('stage6_layout_pressure_poisson_modified_flag') == 0
assert need('stage6_layout_projection_modified_flag') == 0
assert need('stage6_layout_real_projection_called_flag') == 0
assert need('stage6_layout_guard_check_status') == 1

print('STAGE 6.7 LAYOUT GUARD CHECK PASSED')
