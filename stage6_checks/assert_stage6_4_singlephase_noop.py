#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_singlephase_noop_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
assert int(vals['stage6_dns_noop_default_enable_stage6'])==0
assert int(vals['stage6_dns_noop_default_main_rhs_hook_enabled'])==0
assert int(vals['stage6_dns_noop_default_controlled_test_enabled'])==0
assert int(vals['stage6_dns_noop_default_production_two_way_enabled'])==0
assert int(vals['stage6_dns_noop_default_rhs_allowed_flag'])==0
for k in ['stage6_dns_noop_rhs_diff_max','stage6_dns_noop_rhsx_diff_max','stage6_dns_noop_rhsy_diff_max','stage6_dns_noop_rhsz_diff_max','stage6_dns_noop_nonzero_buffer_rhs_diff_max','stage6_dns_noop_velocity_diff_max','stage6_dns_noop_ux_diff_max','stage6_dns_noop_uy_diff_max','stage6_dns_noop_uz_diff_max','stage6_dns_noop_pressure_diff_max','stage6_dns_noop_pressure_rhs_diff_max','stage6_dns_noop_divergence_diff_max','stage6_dns_noop_checksum_abs_error']:
    assert vals[k]<=1e-14
assert vals['stage6_dns_noop_nonzero_buffer_max_abs']>0
assert int(vals['stage6_dns_noop_injected_flag'])==0 and int(vals['stage6_dns_noop_modified_flag'])==0
assert int(vals['stage6_dns_noop_fluid_update_called_flag'])==0
assert int(vals['stage6_dns_noop_projection_modified_flag'])==0 and int(vals['stage6_dns_noop_pressure_poisson_modified_flag'])==0
assert int(vals['stage6_dns_noop_divergence_status'])==1
assert int(vals['stage6_dns_noop_production_enabled_by_default_flag'])==0 and int(vals['stage6_dns_noop_controlled_test_only_flag'])==1 and int(vals['stage6_dns_noop_default_main_dns_safe_flag'])==1
assert int(vals['stage6_dns_noop_regression_check_status'])==1
print('STAGE 6.4 SINGLEPHASE NOOP CHECK PASSED')
