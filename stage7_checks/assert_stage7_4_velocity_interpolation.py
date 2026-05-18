#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_velocity_interpolation_check.dat').read_text().splitlines():
 p=l.split()
 if len(p)>=2: vals[p[0]]=float(p[1])
expected_one=['stage7_vel_stage6_dependency_status','stage7_vel_stage7_0_dependency_status','stage7_vel_stage7_1_dependency_status','stage7_vel_stage7_2_dependency_status','stage7_vel_stage7_3_dependency_status','stage7_vel_collocated_layout_valid_flag','stage7_vel_component_layout_valid_flag','stage7_vel_invalid_layout_rejected_flag','stage7_vel_layout_validation_status','stage7_vel_constant_status','stage7_vel_linear_status','stage7_vel_poiseuille_status','stage7_vel_component_consistency_status','stage7_vel_periodic_x_shift_status','stage7_vel_periodic_z_shift_status','stage7_vel_nearwall_blocked_flag','stage7_vel_y_outside_low_blocked_flag','stage7_vel_y_outside_high_blocked_flag','stage7_vel_y_outside_status','stage7_velocity_interpolation_check_status']
for k in expected_one: assert vals[k]==1.0,f'{k} expected1 got {vals[k]}'
expected_zero=['stage7_vel_noop_rhs_modified_flag','stage7_vel_pressure_poisson_modified_flag','stage7_vel_projection_modified_flag','stage7_vel_production_dns_called_flag','stage7_vel_fluid_update_called_flag','stage7_vel_fibre_advance_called_flag']
for k in expected_zero: assert vals[k]==0.0,f'{k} expected0 got {vals[k]}'
assert vals['stage7_vel_constant_error_max']<=1e-12
assert vals['stage7_vel_linear_error_max']<=1e-11
assert vals['stage7_vel_poiseuille_error_max']<=1e-11
assert vals['stage7_vel_u_component_error_max']<=1e-11
assert vals['stage7_vel_v_component_error_max']<=1e-11
assert vals['stage7_vel_w_component_error_max']<=1e-11
assert vals['stage7_vel_periodic_x_shift_error_max']<=1e-12
assert vals['stage7_vel_periodic_z_shift_error_max']<=1e-12
assert vals['stage7_vel_noop_rhs_change_max']<=1e-14
assert vals['stage7_vel_nearwall_unsafe_count']>0
print('STAGE 7.4 VELOCITY INTERPOLATION CHECK PASSED')
