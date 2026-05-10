#!/usr/bin/env python3
from pathlib import Path
p=Path('stage5_outputs/fibre_stage5_closed_loop_check.dat')
vals={}
for ln in p.read_text().splitlines():
    if not ln.strip():
        continue
    k,v=ln.split(None,1)
    try: vals[k]=float(v)
    except: vals[k]=v
assert int(vals['stage5_closed_loop_check_status'])==1
assert vals['stage5_closed_loop_completed_steps']==vals['stage5_closed_loop_step_count']
assert int(vals['stage5_closed_loop_two_way_enabled_flag'])==1
assert vals['stage5_closed_loop_fluid_velocity_change_norm']>0
assert vals['stage5_closed_loop_fibre_center_displacement_norm']>0
assert vals['stage5_closed_loop_fibre_center_velocity_norm']>0
assert vals['stage5_closed_loop_f_ext_norm_max']>0
assert vals['stage5_closed_loop_force_conservation_relative_error_max']<=1e-10
assert vals['stage5_closed_loop_force_conservation_abs_error_max']<=1e-10
assert vals['stage5_closed_loop_power_relative_error_max']<=1e-10
assert vals['stage5_closed_loop_power_abs_error_max']<=1e-10
assert vals['stage5_closed_loop_momentum_exchange_relative_error_max']<=1e-10
assert vals['stage5_closed_loop_momentum_exchange_abs_error_max']<=1e-10
assert int(vals['stage5_closed_loop_pressure_poisson_modified_flag'])==0
assert int(vals['stage5_closed_loop_main_dns_hooked_flag'])==0
assert int(vals['stage5_closed_loop_synthetic_only_flag'])==1
print('STAGE 5.8 CLOSED LOOP CHECK PASSED')
