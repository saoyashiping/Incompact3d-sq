# P1_2 case manifest

- case: `production_recovery/P1_2_case/input.i3d`
- source template: `production_recovery/P1_1_case/input.i3d`
- itype: 3 real channel
- nx: 96
- ny: 97
- nz: 96
- dt: 2.5e-5
- ilast: 150
- fibre_count: 1
- fibre_nnode: 49
- lambda_fsi: 1.0e-5 and 1.0e-4
- penalty_beta: 2.0
- status: guarded small-lambda two-way validation only; not production DNS.
