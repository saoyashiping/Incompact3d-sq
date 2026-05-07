#!/usr/bin/env python3
from pathlib import Path
vals={}
for ln in Path('stage6_outputs/fibre_stage6_rhs_audit_check.dat').read_text().splitlines():
    if not ln.strip(): continue
    k,v=ln.split(None,1); vals[k]=float(v)
for k in ['stage6_rhs_location_identified','stage6_u_rhs_array_identified','stage6_v_rhs_array_identified','stage6_w_rhs_array_identified','stage6_momentum_rhs_target_status','stage6_rhs_unit_convention_status','stage6_projection_order_status','stage6_rk_substep_policy_status','stage6_rhs_audit_config_safety_status','stage6_rhs_audit_report_exists','stage6_rhs_audit_report_status','stage6_rhs_audit_check_status']:
    assert int(vals[k])==1
assert vals['stage6_rhs_audit_noop_rhs_change_max']<=1e-14
assert int(vals['stage6_rhs_audit_rhs_modified_flag'])==0
r=Path('stage6_outputs/stage6_rhs_audit_report.md').read_text()
for p in ['Stage 6.1 RHS Audit Report','acceleration form','f_ibm / rho_fluid','momentum RHS','before pressure projection','not directly modified','current substep','audit only','Production two-way enabled by default: false']:
    assert p in r
print('STAGE 6.1 RHS AUDIT CHECK PASSED')
