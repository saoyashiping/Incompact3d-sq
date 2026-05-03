#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage4_outputs/fibre_stage4_smoke_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage4_smoke_prior_stage_count']==9
assert v['stage4_smoke_prior_stage_status']==1
assert v['stage4_smoke_safe_unsafe_count']==0
assert v['stage4_smoke_safe_center_displacement_norm']>0
assert v['stage4_smoke_safe_center_velocity_norm']>0
assert v['stage4_smoke_safe_f_ext_norm']>0
assert v['stage4_smoke_safe_length_error']<=1e-8
assert v['stage4_smoke_safe_nan_detected']==0
assert v['stage4_smoke_safe_solver_failure_count']==0
assert v['stage4_smoke_force_conservation_error']<=1e-10
assert v['stage4_smoke_force_conservation_relative_error']<=1e-10
assert v['stage4_smoke_force_buffer_max_abs']>0
assert v['stage4_smoke_power_abs_error']<=1e-10
assert v['stage4_smoke_power_relative_error']<=1e-10
assert abs(v['stage4_smoke_power_abs_error']-abs(v['stage4_smoke_power_eulerian']-v['stage4_smoke_power_lagrangian']))<=1e-12
assert v['stage4_smoke_power_error_consistency_check']<=1e-12
assert v['stage4_smoke_power_nonzero_flag']==1
assert v['stage4_smoke_f_ext_refresh_error_max']<=1e-12
assert v['stage4_smoke_postadvance_force_staleness_norm']>0
assert v['stage4_smoke_velocity_change_max']<=1e-14
assert v['stage4_smoke_nearwall_unsafe_count']>0
assert v['stage4_smoke_nearwall_blocked_flag']==1
assert v['stage4_smoke_nearwall_interpolation_called']==0
assert v['stage4_smoke_nearwall_feedback_called']==0
assert v['stage4_smoke_nearwall_structure_advance_called']==0
assert v['stage4_smoke_nearwall_f_ext_norm']<=1e-14
assert v['stage4_smoke_nonuniform_y_blocked_flag']==1
assert v['stage4_smoke_unknown_layout_blocked_flag']==1
assert v['stage4_smoke_staggered_layout_blocked_flag']==1
assert v['stage4_smoke_nonzero_buffer_rhs_change_max']<=1e-14
assert v['stage4_smoke_fluid_rhs_modified']==0
assert v['stage4_smoke_apply_ibm_to_fluid_rhs']==0
assert v['stage4_smoke_rhs_disabled_flag']==1
assert v['stage4_smoke_status']==1
print('STAGE 4.9 SMOKE CHECK PASSED')
