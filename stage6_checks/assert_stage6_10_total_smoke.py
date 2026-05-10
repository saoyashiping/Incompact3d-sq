#!/usr/bin/env python3
from pathlib import Path
vals={}
for line in Path('stage6_outputs/fibre_stage6_total_smoke_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2:
        vals[p[0]]=float(p[1])

def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]

zeros=['stage6_total_production_two_way_enabled_by_default','stage6_total_default_rhs_allowed_flag','stage6_total_default_main_hook_enabled_flag','stage6_total_rhs_after_projection_flag','stage6_total_pressure_poisson_direct_modify_flag','stage6_total_post_projection_velocity_modified_flag','stage6_total_blocked_rhs_injection_called_count','stage6_total_pressure_poisson_modified_flag','stage6_total_projection_modified_flag','stage6_total_real_projection_called_flag','stage6_total_production_dns_called_flag']
ones=['stage6_total_stage6_0_output_exists','stage6_total_stage6_1_output_exists','stage6_total_stage6_2_output_exists','stage6_total_stage6_3_output_exists','stage6_total_stage6_4_output_exists','stage6_total_stage6_5_output_exists','stage6_total_stage6_6_output_exists','stage6_total_stage6_7_output_exists','stage6_total_stage6_8_output_exists','stage6_total_stage6_9_output_exists','stage6_total_all_prior_outputs_exist','stage6_total_default_config_safe_flag','stage6_total_controlled_injected_flag','stage6_total_controlled_modified_flag','stage6_total_rhs_before_projection_flag','stage6_total_rk_substep_policy_status','stage6_total_rk_current_substep_required_flag','stage6_total_rk_stale_force_forbidden_flag','stage6_total_rk_stale_force_detected_flag','stage6_total_uniform_collocated_supported_flag','stage6_total_nonuniform_y_blocked_flag','stage6_total_staggered_blocked_flag','stage6_total_unknown_layout_blocked_flag','stage6_total_closure_marker_exists','stage6_total_closure_summary_status','stage6_total_smoke_check_status']
for k in zeros: assert need(k)==0, k
for k in ones: assert need(k)==1, k
for k in ['stage6_total_controlled_rhs_expected_error','stage6_total_controlled_component_x_error','stage6_total_controlled_component_y_error','stage6_total_controlled_component_z_error','stage6_total_rk_rhs_match_error_max']:
    assert abs(need(k))<=1e-14, k

m=Path('stage6_outputs/STAGE6_CLOSED.md')
assert m.exists()
text=m.read_text()
for s in ['STAGE 6 CLOSED','default no-op production safety','controlled RHS injection','RHS-before-projection','pressure Poisson not directly modified','RK substep','layout guard','production two-way disabled by default']:
    assert s in text, s
print('STAGE 6.10 TOTAL SMOKE CHECK PASSED')
