#!/usr/bin/env python3
from pathlib import Path
vals = {}
for line in Path('stage6_outputs/fibre_stage6_micro_smoke_check.dat').read_text().splitlines():
    p = line.split()
    if len(p) >= 2:
        vals[p[0]] = float(p[1])

def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]

def z(k, t=1e-14):
    return abs(need(k)) <= t

assert need('stage6_micro_controlled_config_valid_flag') == 1
assert need('stage6_micro_layout_guard_pass_flag') == 1
assert need('stage6_micro_hook_called_flag') == 1
assert need('stage6_micro_rhs_injection_called_flag') == 1
assert need('stage6_micro_controlled_path_status') == 1
assert need('stage6_micro_buffer_max_abs') > 0
assert need('stage6_micro_rhs_change_max') > 0
assert z('stage6_micro_rhs_expected_error')
assert need('stage6_micro_injected_flag') == 1
assert need('stage6_micro_modified_flag') == 1
assert need('stage6_micro_rejected_flag') == 0
assert z('stage6_micro_component_x_error')
assert z('stage6_micro_component_y_error')
assert z('stage6_micro_component_z_error')
assert z('stage6_micro_ustar_expected_error')
assert need('stage6_micro_ustar_change_norm') > 0
assert need('stage6_micro_dt_positive_flag') == 1
assert z('stage6_micro_default_rhs_change_max')
assert z('stage6_micro_default_ustar_diff_max')
assert need('stage6_micro_default_injected_flag') == 0
assert need('stage6_micro_default_modified_flag') == 0
assert need('stage6_micro_blocked_layout_blocked_flag') == 1
assert z('stage6_micro_blocked_layout_rhs_change_max')
assert need('stage6_micro_blocked_layout_injected_flag') == 0
assert need('stage6_micro_blocked_layout_modified_flag') == 0
assert need('stage6_micro_blocked_layout_rhs_injection_called_flag') == 0
assert need('stage6_micro_invalid_rejected_flag') == 1
assert z('stage6_micro_invalid_rhs_change_max')
assert need('stage6_micro_invalid_injected_flag') == 0
assert need('stage6_micro_invalid_modified_flag') == 0
assert need('stage6_micro_production_enabled_rejected_flag') == 1
assert z('stage6_micro_production_enabled_rhs_change_max')
assert need('stage6_micro_production_enabled_injected_flag') == 0
assert need('stage6_micro_production_enabled_modified_flag') == 0
assert z('stage6_micro_zero_buffer_rhs_change_max')
assert need('stage6_micro_zero_buffer_injected_flag') == 1
assert need('stage6_micro_zero_buffer_modified_flag') == 0
assert need('stage6_micro_pressure_poisson_modified_flag') == 0
assert need('stage6_micro_pressure_rhs_modified_flag') == 0
assert need('stage6_micro_projection_modified_flag') == 0
assert need('stage6_micro_real_projection_called_flag') == 0
assert need('stage6_micro_post_projection_velocity_modified_flag') == 0
assert need('stage6_micro_production_dns_called_flag') == 0
assert need('stage6_micro_production_enabled_by_default_flag') == 0
assert need('stage6_micro_smoke_check_status') == 1
print('STAGE 6.8 MICRO SMOKE CHECK PASSED')
