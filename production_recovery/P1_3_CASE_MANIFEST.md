# P1_3 case manifest

- case: `production_recovery/P1_3_case/input.i3d`
- source template: `production_recovery/P1_2_case/input.i3d`
- itype: 3 real channel
- nx: 128
- ny: 129
- nz: 96
- dt: 2.5e-5
- segment1_steps: 150
- restart_segment2_steps: 150
- total_guarded_steps: 300
- fibre_count: 1
- fibre_nnode: 65
- lambda_fsi: 1.0e-4
- penalty_beta: 2.0
- status: guarded restart/statistics/visualization compatibility only; not production DNS.
