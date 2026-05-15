#!/usr/bin/env python3
from pathlib import Path

vals = {}
for line in Path('stage7_outputs/fibre_stage7_rho_convention_audit.dat').read_text().splitlines():
    parts = line.split()
    if len(parts) >= 2:
        vals[parts[0]] = float(parts[1])

expected_one = [
    'stage7_rho_stage6_rhs_hook_called_flag',
    'stage7_rho_rho2_hook_called_flag',
    'stage7_rho_rho4_hook_called_flag',
    'stage7_rho_rho2_injected_flag',
    'stage7_rho_rho4_injected_flag',
    'stage7_rho_force_buffer_independent_of_rho_flag',
    'stage7_rho_rhs_divides_once_flag',
    'stage7_rho_stage6_config_validator_called_flag',
    'stage7_rho_invalid_rho_rejected_by_stage6_flag',
    'stage7_rho_convention_audit_status',
]
for key in expected_one:
    assert vals[key] == 1.0, f'{key} expected 1, got {vals[key]}'

expected_zero = [
    'stage7_rho_double_division_detected_flag',
    'stage7_rho_invalid_rho_rhs_allowed_flag',
]
for key in expected_zero:
    assert vals[key] == 0.0, f'{key} expected 0, got {vals[key]}'

assert vals['stage7_rho_force_buffer_change_with_rho_max'] <= 1e-14
assert vals['stage7_rho_rho2_expected_error_max'] <= 1e-14
assert vals['stage7_rho_rho4_expected_error_max'] <= 1e-14
assert vals['stage7_rho_scaling_error'] <= 1e-14

print('STAGE 7 RHO CONVENTION AUDIT PASSED')
