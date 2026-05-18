#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage8_outputs/fibre_stage8_twoway_force_density_check.dat').read_text().splitlines():
 s=ln.split();
 if len(s)>=2: vals[s[0]]=float(s[1])
ones='''stage8_twoway_stage7_closed_marker_exists stage8_twoway_stage7_total_smoke_output_exists stage8_twoway_stage7_total_smoke_status stage8_twoway_stage7_closed_marker_status stage8_twoway_stage8_0_output_exists stage8_twoway_stage8_0_status stage8_twoway_stage8_1_output_exists stage8_twoway_stage8_1_status stage8_twoway_stage8_2_output_exists stage8_twoway_stage8_2_status stage8_twoway_stage8_3_output_exists stage8_twoway_stage8_3_status stage8_twoway_stage8_4_output_exists stage8_twoway_stage8_4_status stage8_twoway_stage8_5_output_exists stage8_twoway_stage8_5_status stage8_twoway_dependency_status stage8_twoway_force_conservation_status stage8_twoway_component_conservation_status stage8_twoway_force_density_convention_flag stage8_twoway_no_rho_division_flag stage8_twoway_no_rho_status stage8_twoway_nonuniform_volume_used_flag stage8_twoway_volume_scaling_status stage8_twoway_periodic_x_status stage8_twoway_periodic_z_status stage8_twoway_periodic_status stage8_twoway_blocked_status stage8_twoway_zero_force_status stage8_twoway_no_rhs_no_projection_status stage8_twoway_noop_safety_status stage8_twoway_force_density_check_status'''.split()
zeros='''stage8_twoway_safe_blocked_count stage8_twoway_safe_unsafe_count stage8_twoway_rhs_hook_called_flag stage8_twoway_rhs_modified_flag stage8_twoway_pressure_poisson_modified_flag stage8_twoway_projection_modified_flag stage8_twoway_real_projection_called_flag stage8_twoway_production_dns_called_flag stage8_twoway_fluid_update_called_flag stage8_twoway_fibre_advance_called_flag stage8_twoway_noop_rhs_modified_flag'''.split()
for k in ones: assert vals.get(k,0)==1,k
for k in zeros: assert vals.get(k,1)==0,k
assert vals['stage8_twoway_safe_valid_count']>=2
assert vals['stage8_twoway_force_density_norm_max']>1e-14
assert vals['stage8_twoway_force_conservation_abs_error']<=1e-12
assert vals['stage8_twoway_force_conservation_relative_error']<=1e-12
assert vals['stage8_twoway_force_conservation_x_error']<=1e-12
assert vals['stage8_twoway_force_conservation_y_error']<=1e-12
assert vals['stage8_twoway_force_conservation_z_error']<=1e-12
assert vals['stage8_twoway_force_buffer_change_with_rho_max']<=1e-14
assert vals['stage8_twoway_periodic_x_force_error']<=1e-12
assert vals['stage8_twoway_periodic_z_force_error']<=1e-12
assert vals['stage8_twoway_blocked_count']>0
assert vals['stage8_twoway_blocked_force_buffer_norm_max']<=1e-14
assert vals['stage8_twoway_blocked_force_buffer_write_error_max']<=1e-14
assert vals['stage8_twoway_zero_force_buffer_norm_max']<=1e-14
assert vals['stage8_twoway_zero_force_conservation_error']<=1e-14
assert vals['stage8_twoway_noop_rhs_change_max']<=1e-14
print('STAGE 8.6 TWOWAY FORCE-DENSITY CHECK PASSED')
