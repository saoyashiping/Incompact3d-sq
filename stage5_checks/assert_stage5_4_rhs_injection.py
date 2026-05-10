#!/usr/bin/env python3
from pathlib import Path
v={}
for ln in Path('stage5_outputs/fibre_stage5_rhs_injection_check.dat').read_text().splitlines():
 s=ln.split()
 if len(s)>=2:v[s[0]]=float(s[1])
assert v['stage5_injection_zero_buffer_change_max']<=1e-14
assert v['stage5_injection_zero_buffer_expected_error']<=1e-14
assert v['stage5_injection_zero_buffer_injected_flag']==1
assert v['stage5_injection_zero_buffer_modified_flag']==0
assert v['stage5_injection_nonzero_buffer_max_abs']>0
assert v['stage5_injection_nonzero_rhs_change_max']>0
assert v['stage5_injection_nonzero_expected_error']<=1e-14
assert v['stage5_injection_nonzero_injected_flag']==1
assert v['stage5_injection_nonzero_modified_flag']==1
assert v['stage5_injection_component_x_error']<=1e-14
assert v['stage5_injection_component_y_error']<=1e-14
assert v['stage5_injection_component_z_error']<=1e-14
assert v['stage5_injection_disabled_change_max']<=1e-14
assert v['stage5_injection_disabled_injected_flag']==0
assert v['stage5_injection_disabled_modified_flag']==0
assert v['stage5_injection_oneway_change_max']<=1e-14
assert v['stage5_injection_oneway_injected_flag']==0
assert v['stage5_injection_oneway_modified_flag']==0
assert v['stage5_injection_invalid_rejected_flag']==1
assert v['stage5_injection_invalid_change_max']<=1e-14
assert v['stage5_injection_invalid_injected_flag']==0
assert v['stage5_injection_rho2_expected_error']<=1e-14
assert v['stage5_injection_rho4_expected_error']<=1e-14
assert v['stage5_injection_rho_scaling_error']<=1e-14
assert v['stage5_injection_pressure_poisson_modified_flag']==0
assert v['stage5_injection_main_dns_hooked_flag']==0
assert v['stage5_injection_synthetic_only_flag']==1
assert v['stage5_rhs_injection_check_status']==1
print('STAGE 5.4 RHS INJECTION CHECK PASSED')
