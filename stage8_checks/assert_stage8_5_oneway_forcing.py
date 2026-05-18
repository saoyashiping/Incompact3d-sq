#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage8_outputs/fibre_stage8_oneway_forcing_check.dat').read_text().splitlines():
 s=ln.split();
 if len(s)>=2: vals[s[0]]=float(s[1])
ones='''stage8_oneway_stage7_closed_marker_exists stage8_oneway_stage7_total_smoke_output_exists stage8_oneway_stage7_total_smoke_status stage8_oneway_stage7_closed_marker_status stage8_oneway_stage8_0_output_exists stage8_oneway_stage8_0_status stage8_oneway_stage8_1_output_exists stage8_oneway_stage8_1_status stage8_oneway_stage8_2_output_exists stage8_oneway_stage8_2_status stage8_oneway_stage8_3_output_exists stage8_oneway_stage8_3_status stage8_oneway_stage8_4_output_exists stage8_oneway_stage8_4_status stage8_oneway_dependency_status stage8_oneway_constant_status stage8_oneway_zero_slip_status stage8_oneway_linear_status stage8_oneway_poiseuille_status stage8_oneway_galilean_status stage8_oneway_structure_power_status stage8_oneway_blocked_status stage8_oneway_zero_beta_rejected_flag stage8_oneway_negative_beta_rejected_flag stage8_oneway_nan_beta_rejected_flag stage8_oneway_invalid_beta_status stage8_oneway_no_spreading_no_rhs_status stage8_oneway_noop_safety_status stage8_oneway_forcing_check_status'''.split()
zeros='''stage8_oneway_constant_blocked_count stage8_oneway_constant_unsafe_count stage8_oneway_spreading_called_flag stage8_oneway_rhs_hook_called_flag stage8_oneway_rhs_modified_flag stage8_oneway_noop_rhs_modified_flag stage8_oneway_pressure_poisson_modified_flag stage8_oneway_projection_modified_flag stage8_oneway_real_projection_called_flag stage8_oneway_production_dns_called_flag stage8_oneway_fluid_update_called_flag stage8_oneway_fibre_advance_called_flag'''.split()
for k in ones: assert vals.get(k,0)==1,k
for k in zeros: assert vals.get(k,1)==0,k
assert vals['stage8_oneway_constant_valid_count']>=2
assert vals['stage8_oneway_constant_velocity_error_max']<=1e-12
assert vals['stage8_oneway_constant_slip_error_max']<=1e-12
assert vals['stage8_oneway_constant_force_error_max']<=1e-12
assert vals['stage8_oneway_zero_slip_error_max']<=1e-14
assert vals['stage8_oneway_zero_force_error_max']<=1e-14
assert vals['stage8_oneway_linear_velocity_error_max']<=1e-11
assert vals['stage8_oneway_linear_slip_error_max']<=1e-11
assert vals['stage8_oneway_linear_force_error_max']<=1e-11
assert vals['stage8_oneway_poiseuille_velocity_error_max']<=1e-11
assert vals['stage8_oneway_poiseuille_force_error_max']<=1e-11
assert vals['stage8_oneway_galilean_slip_error_max']<=1e-12
assert vals['stage8_oneway_galilean_force_error_max']<=1e-12
assert vals['stage8_oneway_structure_power_error']<=1e-12
assert vals['stage8_oneway_blocked_count']>0
assert vals['stage8_oneway_blocked_velocity_write_error_max']<=1e-14
assert vals['stage8_oneway_blocked_slip_write_error_max']<=1e-14
assert vals['stage8_oneway_blocked_force_write_error_max']<=1e-14
assert vals['stage8_oneway_noop_rhs_change_max']<=1e-14
print('STAGE 8.5 ONEWAY FORCING CHECK PASSED')
