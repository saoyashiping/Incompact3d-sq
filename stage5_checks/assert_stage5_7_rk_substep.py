#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage5_outputs/fibre_stage5_rk_substep_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
for k in ['stage5_rk_velocity_12_difference_norm','stage5_rk_velocity_23_difference_norm','stage5_rk_velocity_13_difference_norm','stage5_rk_u_lag_12_difference_norm','stage5_rk_u_lag_23_difference_norm','stage5_rk_u_lag_13_difference_norm','stage5_rk_force_12_difference_norm','stage5_rk_force_23_difference_norm','stage5_rk_force_13_difference_norm','stage5_rk_buffer_12_difference_norm','stage5_rk_buffer_23_difference_norm','stage5_rk_buffer_13_difference_norm']:
 assert v[k]>0
for k in ['stage5_rk_velocity_distinct_flag','stage5_rk_u_lag_distinct_flag','stage5_rk_force_distinct_flag','stage5_rk_buffer_distinct_flag','stage5_rk_stale_force_detected_flag']:
 assert v[k]==1
for k in ['stage5_rk_force_match_error_substep1','stage5_rk_force_match_error_substep2','stage5_rk_force_match_error_substep3','stage5_rk_force_match_error_max']:
 assert v[k]<=1e-12
for k in ['stage5_rk_rhs_match_error_substep1','stage5_rk_rhs_match_error_substep2','stage5_rk_rhs_match_error_substep3','stage5_rk_rhs_match_error_max']:
 assert v[k]<=1e-14
for k in ['stage5_rk_f_ext_refresh_error_substep1','stage5_rk_f_ext_refresh_error_substep2','stage5_rk_f_ext_refresh_error_substep3','stage5_rk_f_ext_refresh_error_max']:
 assert v[k]<=1e-12
assert v['stage5_rk_stale_force_error_substep2']>0
assert v['stage5_rk_stale_force_error_substep3']>0
assert v['stage5_rk_oneway_force_refresh_error_max']<=1e-12
assert v['stage5_rk_oneway_rhs_injected_flag']==0
assert v['stage5_rk_oneway_rhs_change_max']<=1e-14
assert v['stage5_rk_disabled_interpolation_called']==0
assert v['stage5_rk_disabled_rhs_injected_flag']==0
assert v['stage5_rk_invalid_rejected_flag']==1
assert v['stage5_rk_invalid_rhs_injected_flag']==0
assert v['stage5_rk_nearwall_unsafe_count']>0
assert v['stage5_rk_nearwall_blocked_flag']==1
for k in ['stage5_rk_nearwall_interpolation_called','stage5_rk_nearwall_feedback_called','stage5_rk_nearwall_spreading_called','stage5_rk_nearwall_rhs_injected_flag']:
 assert v[k]==0
assert v['stage5_rk_nearwall_f_ext_norm']<=1e-14
assert v['stage5_rk_pressure_poisson_modified_flag']==0
assert v['stage5_rk_main_dns_hooked_flag']==0
assert v['stage5_rk_synthetic_only_flag']==1
assert v['stage5_rk_substep_check_status']==1
print('STAGE 5.7 RK SUBSTEP CHECK PASSED')
