#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_rk_rhs_sync_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
for k in ['stage6_rk_rhs_substep_policy_status','stage6_rk_rhs_current_substep_required_flag','stage6_rk_rhs_force_recompute_required_flag','stage6_rk_rhs_stale_force_forbidden_flag','stage6_rk_buffer_distinct_flag','stage6_rk_rhs_distinct_flag','stage6_rk_stale_force_detected_flag','stage6_rk_rhs_sync_check_status']:
    assert int(vals[k])==1
for k in ['stage6_rk_buffer_12_difference_norm','stage6_rk_buffer_23_difference_norm','stage6_rk_buffer_13_difference_norm','stage6_rk_rhs_12_difference_norm','stage6_rk_rhs_23_difference_norm','stage6_rk_rhs_13_difference_norm','stage6_rk_stale_force_error_substep2','stage6_rk_stale_force_error_substep3']:
    assert vals[k]>0
for k in ['stage6_rk_rhs_match_error_substep1','stage6_rk_rhs_match_error_substep2','stage6_rk_rhs_match_error_substep3','stage6_rk_rhs_match_error_max','stage6_rk_component_x_error_max','stage6_rk_component_y_error_max','stage6_rk_component_z_error_max','stage6_rk_default_rhs_change_max','stage6_rk_invalid_rhs_change_max','stage6_rk_production_enabled_rhs_change_max','stage6_rk_velocity_change_max']:
    assert vals[k]<=1e-14
assert int(vals['stage6_rk_default_injected_count'])==0 and int(vals['stage6_rk_default_modified_count'])==0
assert int(vals['stage6_rk_invalid_rejected_flag'])==1 and int(vals['stage6_rk_invalid_injected_count'])==0
assert int(vals['stage6_rk_production_enabled_rejected_flag'])==1 and int(vals['stage6_rk_production_enabled_injected_count'])==0
assert int(vals['stage6_rk_fluid_update_called_flag'])==0 and int(vals['stage6_rk_pressure_poisson_modified_flag'])==0 and int(vals['stage6_rk_projection_modified_flag'])==0 and int(vals['stage6_rk_real_projection_called_flag'])==0
print('STAGE 6.6 RK RHS SYNC CHECK PASSED')
