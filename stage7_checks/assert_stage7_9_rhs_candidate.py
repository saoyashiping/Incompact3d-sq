#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_rhs_candidate_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2: vals[p[0]]=float(p[1])
for k in ['stage7_rhs_stage6_dependency_status','stage7_rhs_stage7_0_dependency_status','stage7_rhs_stage7_1_dependency_status','stage7_rhs_stage7_2_dependency_status','stage7_rhs_stage7_3_dependency_status','stage7_rhs_stage7_4_dependency_status','stage7_rhs_stage7_5_dependency_status','stage7_rhs_stage7_6_dependency_status','stage7_rhs_stage7_7_dependency_status','stage7_rhs_stage7_8_dependency_status','stage7_rhs_default_safety_status','stage7_rhs_controlled_hook_called_flag','stage7_rhs_controlled_injected_flag','stage7_rhs_controlled_modified_flag','stage7_rhs_controlled_integration_status','stage7_rhs_rho_scaling_status','stage7_rhs_no_double_division_status','stage7_rhs_blocked_safety_status','stage7_rhs_production_dns_rejected_flag','stage7_rhs_production_twoway_rejected_flag','stage7_rhs_production_safety_status','stage7_rhs_no_projection_no_dns_status','stage7_rhs_candidate_force_conservation_status','stage7_rhs_candidate_check_status']:
 assert vals[k]==1.0,k
for k in ['stage7_rhs_default_injected_flag','stage7_rhs_default_modified_flag','stage7_rhs_controlled_candidate_blocked_count','stage7_rhs_controlled_candidate_unsafe_count','stage7_rhs_double_division_detected_flag','stage7_rhs_blocked_injected_flag','stage7_rhs_production_injected_flag','stage7_rhs_pressure_poisson_modified_flag','stage7_rhs_projection_modified_flag','stage7_rhs_real_projection_called_flag','stage7_rhs_production_dns_called_flag','stage7_rhs_fluid_update_called_flag','stage7_rhs_fibre_advance_called_flag']:
 assert vals[k]==0.0,k
print('STAGE 7.9 RHS CANDIDATE CHECK PASSED')
