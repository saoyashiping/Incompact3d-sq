#!/usr/bin/env python3
from pathlib import Path
p=Path('stage8_outputs/fibre_stage8_lagrangian_state_check.dat')
vals={}
for ln in p.read_text().splitlines():
    s=ln.split()
    if len(s)>=2:
        vals[s[0]]=float(s[1])
ones='''stage8_state_stage7_closed_marker_exists stage8_state_stage7_total_smoke_output_exists stage8_state_stage7_total_smoke_status stage8_state_stage7_closed_marker_status stage8_state_stage8_0_output_exists stage8_state_stage8_0_status stage8_state_stage8_1_output_exists stage8_state_stage8_1_status stage8_state_dependency_status stage8_state_allocate_valid_flag stage8_state_arrays_allocated_flag stage8_state_initialization_status stage8_state_safe_geometry_status stage8_state_safe_classification_status stage8_state_nearwall_blocked_flag stage8_state_nearwall_status stage8_state_outside_y_blocked_flag stage8_state_outside_y_status stage8_state_invalid_nlag_rejected_flag stage8_state_invalid_length_rejected_flag stage8_state_invalid_direction_rejected_flag stage8_state_invalid_nan_rejected_flag stage8_state_invalid_input_status stage8_state_clear_status stage8_state_destroy_deallocated_flag stage8_state_destroy_status stage8_state_placeholders_zero_status stage8_state_noop_safety_status stage8_lagrangian_state_check_status'''.split()
zeros='''stage8_state_allocate_rejected_flag stage8_state_safe_blocked_count stage8_state_safe_unsafe_count stage8_state_noop_rhs_modified_flag stage8_state_pressure_poisson_modified_flag stage8_state_projection_modified_flag stage8_state_real_projection_called_flag stage8_state_production_dns_called_flag stage8_state_fluid_update_called_flag stage8_state_fibre_advance_called_flag'''.split()
for k in ones: assert vals.get(k,0)==1, k
for k in zeros: assert vals.get(k,1)==0, k
assert vals['stage8_state_initial_zero_error_max']<=1e-14
assert vals['stage8_state_safe_total_ds_error']<=1e-12
assert vals['stage8_state_safe_segment_error_max']<=1e-12
assert vals['stage8_state_safe_point_count']>=2
assert vals['stage8_state_nearwall_blocked_count']>0
assert vals['stage8_state_nearwall_unsafe_count']>0
assert vals['stage8_state_outside_y_blocked_count']>0
assert vals['stage8_state_clear_zero_error_max']<=1e-14
assert vals['stage8_state_velocity_placeholder_zero_error']<=1e-14
assert vals['stage8_state_force_placeholder_zero_error']<=1e-14
assert vals['stage8_state_noop_rhs_change_max']<=1e-14
print('STAGE 8.2 LAGRANGIAN STATE CHECK PASSED')
