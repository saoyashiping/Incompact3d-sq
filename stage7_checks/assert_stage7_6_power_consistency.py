#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_power_consistency_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2: vals[p[0]]=float(p[1])
for k in ['stage7_power_stage6_dependency_status','stage7_power_stage7_0_dependency_status','stage7_power_stage7_1_dependency_status','stage7_power_stage7_2_dependency_status','stage7_power_stage7_3_dependency_status','stage7_power_stage7_4_dependency_status','stage7_power_stage7_5_dependency_status','stage7_power_single_consistency_status','stage7_power_multipoint_consistency_status','stage7_power_periodic_x_wrap_status','stage7_power_periodic_z_wrap_status','stage7_power_periodic_consistency_status','stage7_power_nonuniform_volume_used_flag','stage7_power_volume_weighting_status','stage7_power_nearwall_blocked_flag','stage7_power_y_outside_low_blocked_flag','stage7_power_y_outside_high_blocked_flag','stage7_power_y_outside_status','stage7_power_force_density_convention_flag','stage7_power_no_rho_division_flag','stage7_power_consistency_check_status']:
 assert vals[k]==1.0,k
for k in ['stage7_power_multipoint_unsafe_count','stage7_power_stage6_rhs_hook_called_flag','stage7_power_noop_rhs_modified_flag','stage7_power_pressure_poisson_modified_flag','stage7_power_projection_modified_flag','stage7_power_production_dns_called_flag','stage7_power_fluid_update_called_flag','stage7_power_fibre_advance_called_flag']:
 assert vals[k]==0.0,k
assert vals['stage7_power_single_abs_error']<=1e-12
assert vals['stage7_power_single_relative_error']<=1e-12
assert vals['stage7_power_multipoint_valid_count']>=5
assert vals['stage7_power_multipoint_abs_error']<=1e-12
assert vals['stage7_power_multipoint_relative_error']<=1e-12
assert vals['stage7_power_periodic_abs_error']<=1e-12
assert vals['stage7_power_periodic_relative_error']<=1e-12
assert vals['stage7_power_volume_weighting_error_max']<=1e-12
assert vals['stage7_power_nearwall_unsafe_count']>0
assert vals['stage7_power_nearwall_buffer_change_max']<=1e-14
assert vals['stage7_power_noop_rhs_change_max']<=1e-14
print('STAGE 7.6 POWER CONSISTENCY CHECK PASSED')
