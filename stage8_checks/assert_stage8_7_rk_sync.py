#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage8_outputs/fibre_stage8_rk_sync_check.dat').read_text().splitlines():
 s=ln.split();
 if len(s)>=2: vals[s[0]]=float(s[1])
ones='''stage8_rk_stage7_closed_marker_exists stage8_rk_stage7_total_smoke_output_exists stage8_rk_stage7_total_smoke_status stage8_rk_stage7_closed_marker_status stage8_rk_stage8_0_output_exists stage8_rk_stage8_0_status stage8_rk_stage8_1_output_exists stage8_rk_stage8_1_status stage8_rk_stage8_2_output_exists stage8_rk_stage8_2_status stage8_rk_stage8_3_output_exists stage8_rk_stage8_3_status stage8_rk_stage8_4_output_exists stage8_rk_stage8_4_status stage8_rk_stage8_5_output_exists stage8_rk_stage8_5_status stage8_rk_stage8_6_output_exists stage8_rk_stage8_6_status stage8_rk_dependency_status stage8_rk_three_substep_workflow_status stage8_rk_force_recomputed_each_substep_flag stage8_rk_force_recompute_status stage8_rk_no_stale_velocity_flag stage8_rk_no_stale_slip_flag stage8_rk_no_stale_feedback_flag stage8_rk_no_stale_state_status stage8_rk_buffer_clear_each_substep_flag stage8_rk_buffer_clear_status stage8_rk_clear_before_velocity_flag stage8_rk_velocity_before_slip_flag stage8_rk_slip_before_feedback_flag stage8_rk_feedback_before_force_density_flag stage8_rk_force_density_before_rhs_gate_flag stage8_rk_rhs_gate_before_projection_gate_flag stage8_rk_event_order_status stage8_rk_force_conservation_status stage8_rk_blocked_substep_rejected_flag stage8_rk_blocked_substep_status stage8_rk_no_rhs_no_projection_status stage8_rk_noop_safety_status stage8_rk_sync_check_status'''.split()
zeros='''stage8_rk_substep_rejected_count stage8_rk_rhs_hook_called_flag stage8_rk_rhs_modified_flag stage8_rk_pressure_poisson_modified_flag stage8_rk_projection_modified_flag stage8_rk_real_projection_called_flag stage8_rk_production_dns_called_flag stage8_rk_fluid_update_called_flag stage8_rk_fibre_advance_called_flag stage8_rk_noop_rhs_modified_flag'''.split()
for k in ones: assert vals.get(k,0)==1,k
for k in zeros: assert vals.get(k,1)==0,k
assert vals['stage8_rk_nsubstep']==3
assert vals['stage8_rk_substep_valid_count']==3
assert vals['stage8_rk_velocity_interpolation_count']==3
assert vals['stage8_rk_slip_assembly_count']==3
assert vals['stage8_rk_feedback_candidate_count']==3
assert vals['stage8_rk_force_density_candidate_count']==3
assert vals['stage8_rk_force_buffer_clear_count']==3
assert vals['stage8_rk_force_density_signature_min_separation']>1e-12
assert vals['stage8_rk_velocity_signature_min_separation']>1e-12
assert vals['stage8_rk_slip_signature_min_separation']>1e-12
assert vals['stage8_rk_feedback_signature_min_separation']>1e-12
assert vals['stage8_rk_buffer_preclear_nonzero_count']>=2
assert vals['stage8_rk_buffer_clear_zero_error_max']<=1e-14
assert vals['stage8_rk_force_conservation_error_substep_1']<=1e-12
assert vals['stage8_rk_force_conservation_error_substep_2']<=1e-12
assert vals['stage8_rk_force_conservation_error_substep_3']<=1e-12
assert vals['stage8_rk_force_conservation_error_max']<=1e-12
assert vals['stage8_rk_blocked_substep_blocked_count']>0
assert vals['stage8_rk_blocked_substep_force_buffer_norm_max']<=1e-14
assert vals['stage8_rk_noop_rhs_change_max']<=1e-14
print('STAGE 8.7 RK SYNC CHECK PASSED')
