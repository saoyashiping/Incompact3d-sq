#!/usr/bin/env python3
from pathlib import Path
vals={}
for line in Path('stage7_outputs/fibre_stage7_grid_metadata_check.dat').read_text().splitlines():
    p=line.split()
    if len(p)>=2: vals[p[0]]=float(p[1])
def need(k):
    assert k in vals, f'missing {k}'
    return vals[k]
ones=['stage7_grid_stage6_closed_marker_status','stage7_grid_stage7_0_dependency_status','stage7_grid_x_uniform_flag','stage7_grid_z_uniform_flag','stage7_grid_periodic_x_flag','stage7_grid_periodic_z_flag','stage7_grid_y_nonuniform_flag','stage7_grid_y_monotonic_flag','stage7_grid_dy_min_positive_flag','stage7_grid_volume_positive_flag','stage7_grid_ymin_correct_flag','stage7_grid_ymax_correct_flag','stage7_grid_wall_bounds_status','stage7_grid_uniform_y_detected_flag','stage7_grid_invalid_size_rejected_flag','stage7_grid_invalid_dy_rejected_flag','stage7_grid_invalid_boundary_rejected_flag','stage7_grid_metadata_check_status']
zeros=['stage7_grid_y_uniform_flag','stage7_grid_uniform_y_nonuniform_flag','stage7_grid_noop_rhs_modified_flag','stage7_grid_pressure_poisson_modified_flag','stage7_grid_projection_modified_flag','stage7_grid_real_projection_called_flag','stage7_grid_production_dns_called_flag','stage7_grid_fluid_update_called_flag','stage7_grid_fibre_advance_called_flag']
for k in ones: assert need(k)==1, k
for k in zeros: assert need(k)==0, k
assert need('stage7_grid_dy_min')>0
assert need('stage7_grid_dy_max')>need('stage7_grid_dy_min')
assert need('stage7_grid_volume_min')>0
assert abs(need('stage7_grid_total_volume_error'))<=1e-12
assert abs(need('stage7_grid_uniform_y_volume_error'))<=1e-12
assert abs(need('stage7_grid_noop_rhs_change_max'))<=1e-14
print('STAGE 7.1 GRID METADATA CHECK PASSED')
