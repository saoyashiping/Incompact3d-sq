#!/usr/bin/env python3
from pathlib import Path
vals={}
for line in Path('stage8_outputs/fibre_stage8_config_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2:
        vals[p[0]]=float(p[1])

def need(k):
    assert k in vals,f'missing {k}'
    return vals[k]

expected_one=[
 'stage8_config_stage7_closed_marker_exists','stage8_config_stage7_total_smoke_output_exists','stage8_config_stage7_total_smoke_status','stage8_config_stage7_closed_marker_status','stage8_config_stage7_dependency_status',
 'stage8_default_valid_flag','stage8_default_config_status',
 'stage8_controlled_enable_stage8','stage8_controlled_controlled_test_enabled','stage8_controlled_runtime_grid_bridge_enabled','stage8_controlled_lagrangian_state_enabled','stage8_controlled_velocity_to_fibre_enabled','stage8_controlled_valid_flag','stage8_controlled_config_status',
 'stage8_invalid_production_dns_rejected_flag','stage8_invalid_production_twoway_rejected_flag','stage8_invalid_no_controlled_test_rejected_flag','stage8_invalid_rho_rejected_flag','stage8_invalid_config_status',
 'stage8_rhs_candidate_controlled_allowed_flag','stage8_rhs_candidate_production_dns_blocked_flag','stage8_rhs_candidate_production_twoway_blocked_flag','stage8_rhs_candidate_gate_status',
 'stage8_config_noop_safety_status','stage8_config_check_status'
]
expected_zero=[
 'stage8_default_enable_stage8','stage8_default_controlled_test_enabled','stage8_default_production_dns_enabled','stage8_default_production_two_way_enabled','stage8_default_rhs_allowed_flag','stage8_default_rejected_flag',
 'stage8_controlled_production_dns_enabled','stage8_controlled_production_two_way_enabled','stage8_controlled_feedback_candidate_enabled','stage8_controlled_force_density_candidate_enabled','stage8_controlled_rhs_candidate_enabled','stage8_controlled_rhs_allowed_flag','stage8_controlled_rejected_flag',
 'stage8_config_noop_rhs_modified_flag','stage8_config_pressure_poisson_modified_flag','stage8_config_projection_modified_flag','stage8_config_real_projection_called_flag','stage8_config_production_dns_called_flag','stage8_config_fluid_update_called_flag','stage8_config_fibre_advance_called_flag'
]
for k in expected_one:
    assert need(k)==1.0,k
for k in expected_zero:
    assert need(k)==0.0,k
assert abs(need('stage8_config_noop_rhs_change_max'))<=1e-14
print('STAGE 8.0 CONFIG CHECK PASSED')
