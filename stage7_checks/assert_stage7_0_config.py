#!/usr/bin/env python3
from pathlib import Path
vals={}
for line in Path('stage7_outputs/fibre_stage7_config_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2: vals[p[0]]=float(p[1])
def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]
expected_one=['stage7_stage6_closed_marker_status','stage7_stage6_dependency_status','stage7_default_valid_flag','stage7_default_config_status','stage7_controlled_enable_stage7','stage7_controlled_nonuniform_y_ibm_enabled','stage7_controlled_component_specific_layout_enabled','stage7_controlled_real_layout_interpolation_enabled','stage7_controlled_real_layout_spreading_enabled','stage7_controlled_real_layout_rhs_candidate_enabled','stage7_controlled_test_enabled','stage7_controlled_valid_flag','stage7_controlled_rhs_allowed_flag','stage7_controlled_config_status','stage7_invalid_production_dns_rejected_flag','stage7_invalid_production_twoway_rejected_flag','stage7_invalid_interp_without_spread_rejected_flag','stage7_invalid_spread_without_interp_rejected_flag','stage7_invalid_capability_without_controlled_test_rejected_flag','stage7_invalid_rho_rejected_flag','stage7_config_check_status']
expected_zero=['stage7_default_enable_stage7','stage7_default_nonuniform_y_ibm_enabled','stage7_default_component_specific_layout_enabled','stage7_default_real_layout_interpolation_enabled','stage7_default_real_layout_spreading_enabled','stage7_default_real_layout_rhs_candidate_enabled','stage7_default_controlled_test_enabled','stage7_default_production_dns_enabled','stage7_default_production_two_way_enabled','stage7_default_rejected_flag','stage7_default_rhs_allowed_flag','stage7_controlled_production_dns_enabled','stage7_controlled_production_two_way_enabled','stage7_controlled_rejected_flag','stage7_invalid_production_dns_valid_flag','stage7_invalid_production_dns_rhs_allowed_flag','stage7_invalid_production_twoway_valid_flag','stage7_invalid_production_twoway_rhs_allowed_flag','stage7_invalid_interp_without_spread_valid_flag','stage7_invalid_spread_without_interp_valid_flag','stage7_invalid_capability_without_controlled_test_rhs_allowed_flag','stage7_invalid_rho_valid_flag','stage7_invalid_rho_rhs_allowed_flag','stage7_config_noop_rhs_modified_flag','stage7_config_pressure_poisson_modified_flag','stage7_config_projection_modified_flag','stage7_config_real_projection_called_flag','stage7_config_production_dns_called_flag','stage7_config_fluid_update_called_flag','stage7_config_fibre_advance_called_flag']
for k in expected_one: assert need(k)==1, k
for k in expected_zero: assert need(k)==0, k
assert abs(need('stage7_config_noop_rhs_change_max'))<=1e-14
assert Path('stage6_outputs/STAGE6_CLOSED.md').exists()
print('STAGE 7.0 CONFIG CHECK PASSED')
