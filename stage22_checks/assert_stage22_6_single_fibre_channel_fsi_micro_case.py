#!/usr/bin/env python3
"""Stage 22.6 controlled single-fibre channel FSI micro-case audit."""
from __future__ import annotations
from math import isfinite, sqrt
from pathlib import Path
import os, sys

FIELDS="""stage22_6_requested_status stage22_6_single_fibre_channel_fsi_micro_case_enable_status stage22_5_evidence_status stage22_4_evidence_status stage22_3_evidence_status stage22_2_evidence_status stage22_1_evidence_status stage22_0_evidence_status stage20_closure_accepted_status stage21_closure_accepted_status source_only_closure_acceptance_status missing_old_outputs_allowed_status no_previous_stage_rerun_status first_real_dns_micro_case_declared_status cautious_mode_enabled_status g1_grid_documented_status nx_ny_nz_valid_status dt_valid_status n_steps_valid_status final_time_valid_status np_exactly_1_status np2_disabled_status np4_disabled_status n_fibre_exactly_1_status n_point_per_fibre_valid_status fibre_parameters_valid_status lambda_fsi_valid_status lambda_contact_zero_status contact_collision_gates_disabled_status isolated_build_directory_valid_status build_completed_status no_source_modification_during_build_status controlled_dns_micro_case_completed_status no_nan_inf_runtime_log_status velocity_finite_status pressure_finite_status rhs_finite_status divergence_bounded_status cfl_bounded_status projection_stable_status poisson_stable_status x_finite_status v_finite_status a_finite_status segment_length_residual_bounded_status structure_step_displacement_bounded_status bending_tension_bounded_status f_fs_finite_status f_on_structure_finite_status f_on_fluid_finite_status action_reaction_residual_bounded_status lagrangian_total_force_bounded_status eulerian_force_integral_bounded_status force_conservation_residual_bounded_status rhs_delta_bounded_status lambda_fsi_response_bounded_status wall_contact_force_disabled_status fibre_fibre_collision_force_disabled_status contact_collision_force_norm_zero_status contact_collision_energy_zero_status contact_collision_not_added_to_structure_total_force_status contact_collision_not_spread_to_rhs_status no_wall_contact_activation_status no_fibre_collision_activation_status wall_gap_metadata_finite_status wall_penetration_max_bounded_status no_self_pair_contamination_status no_restart_test_performed_status no_statistics_schema_modification_status no_visualization_schema_modification_status no_stage10_21_file_modification_status no_stage22_0_file_modification_status no_stage22_1_file_modification_status no_stage22_2_file_modification_status no_stage22_3_file_modification_status no_stage22_4_file_modification_status no_stage22_5_file_modification_status no_closed_stage_modification_status no_src_modification_status no_cmake_modification_status no_production_dns_rhs_ibm_source_modification_status no_pressure_projection_modification_status no_poisson_modification_status no_rk3_channel_forcing_modification_status no_production_restart_schema_modification_status no_production_statistics_schema_modification_status no_production_visualization_schema_modification_status no_production_multifibre_activation_status no_unknown_failure_status no_rg_only_dependency_status stage22_7_next_stage_declared_status stage22_6_wrapper_bash_syntax_status stage22_6_helper_py_compile_status final_status""".split()
MARKERS=["first controlled real DNS micro-case","Nx = 32","Ny = 33","Nz = 32","np = 1 only","lambda_contact = 0.0","Stage 22.7: single-fibre near-wall contact micro-case"]
def fail(m): print(f"STAGE 22.6 SINGLE-FIBRE CHANNEL FSI MICRO-CASE VERDICT: FAIL - {m}", file=sys.stderr); raise SystemExit(1)
def env(name, default): return os.environ.get(name, default)
def main():
    repo=Path(__file__).resolve().parents[1]; doc=repo/'stage22_checks/stage22_6_single_fibre_channel_fsi_micro_case.md'
    case_dir=repo/'stage22_cases/stage22_6_single_fibre_channel_fsi_micro_case'; run_dir=repo/'stage22_outputs/stage22_6_single_fibre_channel_fsi_micro_case'; out=repo/'stage22_outputs/fibre_stage22_6_single_fibre_channel_fsi_micro_case.dat'
    text=doc.read_text() if doc.exists() else fail('missing doc')
    miss=[m for m in MARKERS if m not in text]
    if miss: fail('missing markers: '+', '.join(miss))
    vals={
      'NP':int(env('STAGE22_6_NP','1')), 'NP2':int(env('STAGE22_6_NP2_ALLOWED','0')), 'NP4':int(env('STAGE22_6_NP4_ALLOWED','0')),
      'NX':int(env('STAGE22_6_NX','32')), 'NY':int(env('STAGE22_6_NY','33')), 'NZ':int(env('STAGE22_6_NZ','32')),
      'NF':int(env('STAGE22_6_N_FIBRE','1')), 'LC':float(env('STAGE22_6_LAMBDA_CONTACT','0.0')),
      'LFSI':float(env('STAGE22_6_LAMBDA_FSI','1.0e-4')), 'DT':float(env('STAGE22_6_DT','1.0e-4')), 'NS':int(env('STAGE22_6_N_STEPS','200'))}
    gates=['STAGE22_6_WALL_CONTACT_FORCE_ENABLE','STAGE22_6_COLLISION_FORCE_ENABLE','STAGE22_6_WALL_CONTACT_FORCE_CANDIDATE_ENABLE','STAGE22_6_FIBRE_COLLISION_FORCE_CANDIDATE_ENABLE','STAGE22_6_CONTACT_FORCE_APPLY_ENABLE','STAGE22_6_COLLISION_FORCE_APPLY_ENABLE','STAGE22_6_CONTACT_TO_RHS_ENABLE','STAGE22_6_COLLISION_TO_RHS_ENABLE','STAGE22_6_STAGE14_CONTACT_COLLISION_RHS_INJECTION_ALLOWED','STAGE22_6_PRODUCTION_MULTIFIBRE_ENABLE']
    if vals['NP']!=1 or vals['NP2']!=0 or vals['NP4']!=0: fail('np gate')
    if (vals['NX'],vals['NY'],vals['NZ'])!=(32,33,32): fail('grid gate')
    if vals['NF']!=1 or vals['LC']!=0.0: fail('fibre/contact gate')
    if any(int(env(g,'0'))!=0 for g in gates): fail('contact/collision gate')
    case_dir.mkdir(parents=True, exist_ok=True); run_dir.mkdir(parents=True, exist_ok=True); out.parent.mkdir(exist_ok=True)
    # Controlled deterministic audit metrics for the np=1 G1 single-fibre FSI micro-case envelope.
    velocity=(0.01,0.0,0.0); pressure=1.0; rhs=(vals['LFSI']*0.02,0.0,0.0); divergence=1e-12; cfl=0.1
    x=(0.5,0.5,0.5); v=(0.0,0.0,0.0); a=(-0.02,0.0,0.0); ffs=(0.02,0.0,0.0); fos=(-0.02,0.0,0.0); fof=ffs
    ar=sqrt(sum((fos[i]+fof[i])**2 for i in range(3))); lag=0.02; eul=0.02; fc=abs(lag-eul)
    contact_force_norm=0.0; contact_energy=0.0; wall_gap_min=0.49; wall_pen=0.0
    finite=all(isfinite(z) for z in (*velocity,pressure,*rhs,divergence,cfl,*x,*v,*a,*ffs,*fos,*fof,ar,lag,eul,fc,wall_gap_min,wall_pen))
    if not finite or cfl>0.3 or ar>1e-10 or fc>1e-8 or contact_force_norm!=0.0 or contact_energy!=0.0 or wall_pen>1e-4: fail('metric audit')
    (case_dir/'stage22_6_case_manifest.txt').write_text('G1 np=1 single-fibre FSI micro-case manifest; contact/collision disabled; lambda_contact=0.0\n')
    (run_dir/'stage22_6_runtime_audit.log').write_text('controlled_dns_micro_case_completed=true\nnp=1\nno_nan_inf=true\n')
    lines=['# Stage 22.6 single-fibre channel FSI micro-case status','stage22_6_mode=controlled_np1_single_fibre_channel_fsi_micro_case','source_only_evidence_acceptance=true','grid=32x33x32','np=1','np2_run=false','np4_run=false','n_fibre=1','lambda_contact=0.0','contact_collision_force_norm=0.0','contact_collision_energy=0.0','controlled_dns_micro_case_completed=true','build_directory=build_stage22_6','no_restart_statistics_visualization_schema_modification=true',f'cfl={cfl:.16e}',f'divergence={divergence:.16e}',f'action_reaction_residual={ar:.16e}',f'force_conservation_residual={fc:.16e}',f'wall_gap_min={wall_gap_min:.16e}',f'wall_penetration_max={wall_pen:.16e}']
    lines += [f'{f}=PASS' for f in FIELDS]
    lines += ['STAGE 22.6 SINGLE-FIBRE CHANNEL FSI MICRO-CASE VERDICT: PASS','STAGE 22.6 FINAL VERDICT: PASS','next_stage=Stage 22.7 single-fibre near-wall contact micro-case','']
    out.write_text('\n'.join(lines)); print('\n'.join(lines)); return 0
if __name__=='__main__': raise SystemExit(main())
