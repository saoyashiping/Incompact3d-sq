#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2: vals[s[0]]=float(s[1])
ones='''stage8_v2f_stage7_closed_marker_exists stage8_v2f_stage7_total_smoke_output_exists stage8_v2f_stage7_total_smoke_status stage8_v2f_stage7_closed_marker_status stage8_v2f_stage8_0_output_exists stage8_v2f_stage8_0_status stage8_v2f_stage8_1_output_exists stage8_v2f_stage8_1_status stage8_v2f_stage8_2_output_exists stage8_v2f_stage8_2_status stage8_v2f_dependency_status stage8_v2f_constant_status stage8_v2f_linear_status stage8_v2f_poiseuille_status stage8_v2f_periodic_x_shift_status stage8_v2f_periodic_z_shift_status stage8_v2f_periodic_status stage8_v2f_nearwall_blocked_flag stage8_v2f_nearwall_status stage8_v2f_outside_y_blocked_flag stage8_v2f_outside_y_status stage8_v2f_invalid_layout_blocked_flag stage8_v2f_invalid_layout_status stage8_v2f_other_placeholders_zero_status stage8_v2f_clear_geometry_preserved_flag stage8_v2f_clear_status stage8_v2f_noop_safety_status stage8_velocity_to_fibre_check_status'''.split()
zeros='''stage8_v2f_constant_blocked_count stage8_v2f_constant_unsafe_count stage8_v2f_noop_rhs_modified_flag stage8_v2f_pressure_poisson_modified_flag stage8_v2f_projection_modified_flag stage8_v2f_real_projection_called_flag stage8_v2f_production_dns_called_flag stage8_v2f_fluid_update_called_flag stage8_v2f_fibre_advance_called_flag'''.split()
for k in ones: assert vals.get(k,0)==1,k
for k in zeros: assert vals.get(k,1)==0,k
assert vals['stage8_v2f_constant_valid_count']>=2
assert vals['stage8_v2f_constant_error_max']<=1e-12
assert vals['stage8_v2f_linear_error_max']<=1e-11
assert vals['stage8_v2f_poiseuille_error_max']<=1e-11
assert vals['stage8_v2f_periodic_x_shift_error_max']<=1e-12
assert vals['stage8_v2f_periodic_z_shift_error_max']<=1e-12
assert vals['stage8_v2f_nearwall_blocked_count']>0
assert vals['stage8_v2f_nearwall_unsafe_count']>0
assert vals['stage8_v2f_nearwall_velocity_write_error_max']<=1e-14
assert vals['stage8_v2f_outside_y_blocked_count']>0
assert vals['stage8_v2f_outside_y_velocity_write_error_max']<=1e-14
assert vals['stage8_v2f_invalid_layout_blocked_count']>0
assert vals['stage8_v2f_invalid_layout_velocity_write_error_max']<=1e-14
assert vals['stage8_v2f_v_fibre_zero_error_max']<=1e-14
assert vals['stage8_v2f_slip_zero_error_max']<=1e-14
assert vals['stage8_v2f_force_zero_error_max']<=1e-14
assert vals['stage8_v2f_clear_velocity_error_max']<=1e-14
assert vals['stage8_v2f_noop_rhs_change_max']<=1e-14
print('STAGE 8.3 VELOCITY-TO-FIBRE CHECK PASSED')
