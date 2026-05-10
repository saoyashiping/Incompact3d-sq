#!/usr/bin/env python3
from pathlib import Path
v={}
for l in Path('stage7_outputs/fibre_stage7_structure_high_gamma_check.dat').read_text().splitlines():
 p=l.split();
 if len(p)>=2:v[p[0]]=float(p[1])
assert v['stage7_highgamma_nan_detected']==0
assert v['stage7_highgamma_solver_failure_count']==0
assert v['stage7_highgamma_length_error_max']<=1e-8
assert v['stage7_highgamma_stretch_error_max']<=1e-8
assert v['stage7_highgamma_energy_finite_flag']==1
assert v['stage7_highgamma_curvature_finite_flag']==1
assert v['stage7_highgamma_momentum_finite_flag']==1
assert v['stage7_highgamma_structure_check_status']==1
print('STAGE 7 HIGH GAMMA CHECK PASSED')
