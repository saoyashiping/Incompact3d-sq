#!/usr/bin/env python3
from pathlib import Path

vals = {}
for line in Path('stage5_outputs/fibre_stage5_closed_loop_diagnostics_audit.dat').read_text().splitlines():
    parts = line.split()
    if len(parts) >= 2:
        vals[parts[0]] = float(parts[1])

expected_one = [
    'stage5_closed_loop_eulerian_force_computed_flag',
    'stage5_closed_loop_lagrangian_force_computed_flag',
    'stage5_closed_loop_eulerian_power_computed_flag',
    'stage5_closed_loop_lagrangian_power_computed_flag',
    'stage5_closed_loop_used_spreading_buffer_flag',
    'stage5_closed_loop_used_lagrangian_force_flag',
    'stage5_closed_loop_used_interpolated_velocity_flag',
    'stage5_closed_loop_real_spreading_module_called_flag',
    'stage5_closed_loop_real_interpolation_module_called_flag',
    'stage5_closed_loop_no_tautological_force_flag',
    'stage5_closed_loop_no_tautological_power_flag',
    'stage5_closed_loop_diagnostics_audit_status',
]
for key in expected_one:
    assert vals[key] == 1.0, f'{key} expected 1, got {vals[key]}'

assert vals['stage5_closed_loop_actual_force_conservation_abs_error'] <= 1e-12
assert vals['stage5_closed_loop_actual_force_conservation_relative_error'] <= 1e-12
assert vals['stage5_closed_loop_actual_power_abs_error'] <= 1e-12
assert vals['stage5_closed_loop_actual_power_relative_error'] <= 1e-12

print('STAGE 5 CLOSED-LOOP DIAGNOSTICS AUDIT PASSED')
