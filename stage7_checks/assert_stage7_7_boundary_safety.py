#!/usr/bin/env python3
from pathlib import Path
vals={}
for l in Path('stage7_outputs/fibre_stage7_boundary_safety_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2: vals[p[0]]=float(p[1])
for k in ['stage7_boundary_stage6_dependency_status','stage7_boundary_stage7_0_dependency_status','stage7_boundary_stage7_1_dependency_status','stage7_boundary_stage7_2_dependency_status','stage7_boundary_stage7_3_dependency_status','stage7_boundary_stage7_4_dependency_status','stage7_boundary_stage7_5_dependency_status','stage7_boundary_stage7_6_dependency_status','stage7_boundary_safe_interior_allowed_flag','stage7_boundary_safe_interior_status','stage7_boundary_nearwall_blocked_flag','stage7_boundary_nearwall_status','stage7_boundary_y_low_blocked_flag','stage7_boundary_y_high_blocked_flag','stage7_boundary_y_outside_status','stage7_boundary_invalid_coord_blocked_flag','stage7_boundary_invalid_coord_status','stage7_boundary_invalid_layout_blocked_flag','stage7_boundary_invalid_layout_status','stage7_boundary_periodic_x_allowed_flag','stage7_boundary_periodic_z_allowed_flag','stage7_boundary_periodic_allowed_status','stage7_boundary_blocked_no_buffer_write_flag','stage7_boundary_safety_check_status']:
 assert vals[k]==1.0,k
for k in ['stage7_boundary_noop_rhs_modified_flag','stage7_boundary_pressure_poisson_modified_flag','stage7_boundary_projection_modified_flag','stage7_boundary_production_dns_called_flag','stage7_boundary_fluid_update_called_flag','stage7_boundary_fibre_advance_called_flag']:
 assert vals[k]==0.0,k
assert vals['stage7_boundary_safe_point_count']>=5
assert vals['stage7_boundary_nearwall_blocked_count']>0
assert vals['stage7_boundary_nearwall_buffer_change_max']<=1e-14
assert vals['stage7_boundary_blocked_buffer_change_max']<=1e-14
assert vals['stage7_boundary_noop_rhs_change_max']<=1e-14
print('STAGE 7.7 BOUNDARY SAFETY CHECK PASSED')
