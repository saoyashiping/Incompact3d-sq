#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_projection_order_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert int(vals['stage6_projection_rhs_before_projection_flag'])==1
assert int(vals['stage6_projection_rhs_after_projection_flag'])==0
assert int(vals['stage6_projection_after_rhs_policy_flag'])==1
assert int(vals['stage6_projection_order_status'])==1
assert int(vals['stage6_projection_momentum_rhs_modified_flag'])==1
assert int(vals['stage6_projection_pressure_rhs_modified_flag'])==0
assert int(vals['stage6_projection_pressure_poisson_direct_modify_flag'])==0
for k in ['stage6_projection_rhs_expected_error','stage6_projection_component_x_error','stage6_projection_component_y_error','stage6_projection_component_z_error','stage6_projection_ustar_expected_error','stage6_projection_default_rhs_change_max','stage6_projection_default_ustar_change_max','stage6_projection_invalid_rhs_change_max','stage6_projection_production_enabled_rhs_change_max','stage6_projection_pressure_rhs_diff_max']:
    assert vals[k]<=1e-14
assert vals['stage6_projection_ustar_change_norm']>0 and int(vals['stage6_projection_dt_positive_flag'])==1
assert int(vals['stage6_projection_post_projection_velocity_modified_flag'])==0 and int(vals['stage6_projection_post_projection_direct_forcing_forbidden_flag'])==1
assert int(vals['stage6_projection_default_injected_flag'])==0
assert int(vals['stage6_projection_invalid_rejected_flag'])==1 and int(vals['stage6_projection_invalid_injected_flag'])==0
assert int(vals['stage6_projection_production_enabled_rejected_flag'])==1 and int(vals['stage6_projection_production_enabled_injected_flag'])==0
assert int(vals['stage6_projection_pressure_poisson_modified_flag'])==0 and int(vals['stage6_projection_projection_modified_flag'])==0 and int(vals['stage6_projection_real_projection_called_flag'])==0
assert int(vals['stage6_projection_order_check_status'])==1
print('STAGE 6.5 PROJECTION ORDER CHECK PASSED')
