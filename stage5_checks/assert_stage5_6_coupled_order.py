#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage5_outputs/fibre_stage5_coupled_order_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
for k in ['stage5_order_boundary_before_interpolation_flag','stage5_order_interpolation_before_feedback_flag','stage5_order_feedback_before_spreading_flag','stage5_order_spreading_before_rhs_injection_flag','stage5_order_rhs_injection_before_fluid_update_flag','stage5_order_f_ext_before_structure_advance_flag','stage5_order_two_way_sequence_status']:
 assert v[k]==1
assert v['stage5_order_f_ext_match_error']<=1e-12
assert v['stage5_order_rhs_expected_error']<=1e-14
assert v['stage5_order_velocity_update_error']<=1e-14
assert v['stage5_order_fluid_update_called']==1
assert v['stage5_order_structure_advance_called']==1
assert v['stage5_order_fluid_velocity_change_norm']>0
assert v['stage5_order_structure_motion_norm']>0
assert v['stage5_order_length_error']<=1e-8
assert v['stage5_order_nan_detected']==0
assert v['stage5_order_solver_failure_count']==0
assert v['stage5_order_postadvance_force_staleness_norm']>0
assert v['stage5_order_oneway_structure_advance_called']==1
assert v['stage5_order_oneway_rhs_injected_flag']==0
assert v['stage5_order_oneway_fluid_update_called']==0
assert v['stage5_order_oneway_velocity_change_norm']<=1e-14
for k in ['stage5_order_disabled_interpolation_called','stage5_order_disabled_feedback_called','stage5_order_disabled_spreading_called','stage5_order_disabled_rhs_injected_flag','stage5_order_disabled_fluid_update_called','stage5_order_disabled_structure_advance_called']:
 assert v[k]==0
assert v['stage5_order_invalid_rejected_flag']==1
assert v['stage5_order_invalid_rhs_injected_flag']==0
assert v['stage5_order_invalid_fluid_update_called']==0
assert v['stage5_order_invalid_structure_advance_called']==0
assert v['stage5_order_nearwall_unsafe_count']>0
assert v['stage5_order_nearwall_blocked_flag']==1
for k in ['stage5_order_nearwall_interpolation_called','stage5_order_nearwall_feedback_called','stage5_order_nearwall_spreading_called','stage5_order_nearwall_rhs_injected_flag','stage5_order_nearwall_fluid_update_called','stage5_order_nearwall_structure_advance_called']:
 assert v[k]==0
assert v['stage5_order_nearwall_f_ext_norm']<=1e-14
assert v['stage5_order_pressure_poisson_modified_flag']==0
assert v['stage5_order_main_dns_hooked_flag']==0
assert v['stage5_order_synthetic_only_flag']==1
assert v['stage5_coupled_order_check_status']==1
print('STAGE 5.6 COUPLED ORDER CHECK PASSED')
