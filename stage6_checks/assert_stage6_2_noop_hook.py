#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_noop_hook_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert int(vals['stage6_noop_default_hook_called_flag'])==1
for k in ['stage6_noop_default_rhs_modified_flag','stage6_noop_default_injected_flag','stage6_noop_oneway_rhs_modified_flag','stage6_noop_oneway_injected_flag','stage6_noop_invalid_injected_flag','stage6_noop_invalid_rhs_modified_flag','stage6_noop_production_enabled_injected_flag','stage6_noop_invalid_rho_injected_flag','stage6_noop_fluid_update_called_flag','stage6_noop_pressure_poisson_modified_flag','stage6_noop_projection_modified_flag','stage6_noop_pressure_rhs_modified_flag']:
    assert int(vals[k])==0
for k in ['stage6_noop_default_rhs_change_max','stage6_noop_oneway_rhs_change_max','stage6_noop_invalid_rhs_change_max','stage6_noop_production_enabled_rhs_change_max','stage6_noop_invalid_rho_rhs_change_max','stage6_noop_velocity_change_max']:
    assert vals[k]<=1e-14
assert int(vals['stage6_noop_invalid_rejected_flag'])==1
assert int(vals['stage6_noop_production_enabled_rejected_flag'])==1
assert int(vals['stage6_noop_invalid_rho_rejected_flag'])==1
assert int(vals['stage6_noop_default_main_dns_safe_flag'])==1
assert int(vals['stage6_noop_production_enabled_by_default_flag'])==0
assert int(vals['stage6_noop_controlled_test_only_flag'])==1
assert int(vals['stage6_noop_hook_check_status'])==1
print('STAGE 6.2 NOOP HOOK CHECK PASSED')
