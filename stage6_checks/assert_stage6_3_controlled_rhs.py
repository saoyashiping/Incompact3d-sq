#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_controlled_rhs_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert int(vals['stage6_controlled_hook_called_flag'])==1
assert vals['stage6_controlled_expected_error']<=1e-14
for k in ['stage6_controlled_component_x_error','stage6_controlled_component_y_error','stage6_controlled_component_z_error','stage6_controlled_rho2_expected_error','stage6_controlled_rho4_expected_error','stage6_controlled_rho_scaling_error']:
    assert vals[k]<=1e-14
assert int(vals['stage6_controlled_injected_flag'])==1 and int(vals['stage6_controlled_modified_flag'])==1 and int(vals['stage6_controlled_rejected_flag'])==0
assert vals['stage6_controlled_zero_buffer_rhs_change_max']<=1e-14 and int(vals['stage6_controlled_zero_buffer_injected_flag'])==1 and int(vals['stage6_controlled_zero_buffer_modified_flag'])==0
assert vals['stage6_controlled_default_rhs_change_max']<=1e-14 and int(vals['stage6_controlled_default_injected_flag'])==0
assert int(vals['stage6_controlled_invalid_rejected_flag'])==1
assert int(vals['stage6_controlled_production_enabled_rejected_flag'])==1
assert int(vals['stage6_controlled_invalid_rho_rejected_flag'])==1
assert vals['stage6_controlled_velocity_change_max']<=1e-14 and int(vals['stage6_controlled_fluid_update_called_flag'])==0
assert int(vals['stage6_controlled_pressure_poisson_modified_flag'])==0 and int(vals['stage6_controlled_projection_modified_flag'])==0 and int(vals['stage6_controlled_pressure_rhs_modified_flag'])==0
assert int(vals['stage6_controlled_rhs_hook_check_status'])==1
print('STAGE 6.3 CONTROLLED RHS CHECK PASSED')
