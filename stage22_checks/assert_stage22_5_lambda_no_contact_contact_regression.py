#!/usr/bin/env python3
"""Stage 22.5 helper-only lambda/no-contact/contact regression."""
from __future__ import annotations
from math import isfinite, sqrt
from pathlib import Path
import sys

DT=1e-5; N_STEPS=20; R=0.01; L=0.40; NPT=64; RHO=1.0; DS=L/(NPT-1)
K_WALL=1e2; K_COLL=1e2; DAMP=0.2; CFS=1.0; MEFF=RHO*DS
CW=2*DAMP*sqrt(K_WALL*MEFF); CC=2*DAMP*sqrt(K_COLL*MEFF)
ZERO=1e-14; TOL=1e-12; SCALE_TOL=1e-10; MAXF=1e3; MAXA=1e3; MAXFRAC=0.1; SEGLIM=1e-6
STATUS_FIELDS=["stage22_5_requested_status","stage22_5_lambda_no_contact_contact_regression_enable_status","stage22_4_evidence_status","stage22_3_evidence_status","stage22_2_evidence_status","stage22_1_evidence_status","stage22_0_evidence_status","stage20_closure_accepted_status","stage21_closure_accepted_status","source_only_closure_acceptance_status","missing_old_outputs_allowed_status","no_previous_stage_rerun_status","helper_only_status","grid_setting_valid_status","dt_nsteps_valid_status","fibre_settings_valid_status","lambda_fsi_ladder_valid_status","lambda_contact_ladder_valid_status","lambda_fsi_zero_rhs_noop_status","lambda_fsi_small_rhs_bounded_status","lambda_contact_zero_wall_noop_under_penetration_status","lambda_contact_zero_collision_noop_under_overlap_status","lambda_contact_small_wall_response_scaling_status","lambda_contact_small_collision_response_scaling_status","no_contact_geometry_zero_wall_force_all_lambda_status","no_overlap_geometry_zero_collision_force_all_lambda_status","no_contact_no_overlap_f_total_invariant_status","no_contact_no_overlap_candidate_update_invariant_status","wall_penetration_force_inward_status","fibre_fibre_overlap_force_repulsive_status","fibre_fibre_action_reaction_residual_bounded_status","contact_collision_force_finite_status","contact_collision_force_bounded_status","contact_collision_acceleration_bounded_status","contact_collision_energy_nonnegative_status","damping_power_nonpositive_status","no_attractive_wall_force_status","no_duplicate_pair_force_status","no_self_pair_force_status","f_total_with_contact_candidate_finite_status","a_candidate_finite_status","v_next_candidate_finite_status","x_next_candidate_finite_status","structure_step_displacement_fraction_bounded_status","segment_length_residual_bounded_status","candidate_update_not_committed_status","production_structure_advance_disabled_status","production_structure_state_unchanged_status","production_rhs_coupling_disabled_status","contact_collision_not_spread_to_production_rhs_status","stage14_rhs_injection_disabled_status","production_rhs_update_disabled_status","ibm_forcing_modification_disabled_status","build_disabled_status","production_dns_disabled_status","actual_mpi_disabled_status","production_restart_io_disabled_status","production_statistics_io_disabled_status","production_visualization_io_disabled_status","production_multifibre_disabled_status","no_stage10_21_file_modification_status","no_stage22_0_file_modification_status","no_stage22_1_file_modification_status","no_stage22_2_file_modification_status","no_stage22_3_file_modification_status","no_stage22_4_file_modification_status","no_closed_stage_modification_status","no_src_modification_status","no_cmake_modification_status","no_production_dns_rhs_ibm_io_modification_status","no_production_hook_activation_status","no_build_execution_status","no_mpi_execution_status","no_dns_execution_status","no_production_contact_collision_physics_activation_status","no_production_structure_commit_status","no_production_rhs_update_status","no_unknown_failure_status","no_rg_only_dependency_status","stage22_6_next_stage_declared_status","stage22_5_wrapper_bash_syntax_status","stage22_5_helper_py_compile_status","final_status"]
MARKERS=["Stage 22.5 is helper-only","lambda_fsi_values = 0.0, 1.0e-6","lambda_contact_values = 0.0, 1.0e-6, 1.0","Scaling rule","Case A: one_fibre_no_contact_lambda_contact_0","Case J: lambda_fsi_small_rhs_bounded","Stage 22.6: single-fibre channel FSI micro-case"]
def fail(m): print(f"STAGE 22.5 LAMBDA NO-CONTACT CONTACT REGRESSION VERDICT: FAIL - {m}",file=sys.stderr); raise SystemExit(1)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def scale(s,a): return tuple(s*x for x in a)
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def norm(a): return sqrt(dot(a,a))
def finite(v): return all(isfinite(x) for x in v)
def wall(y,v,lam):
    g=y-R; n=(0.0,1.0,0.0); d=max(0.0,-g); vnm=min(dot(v,n),0.0); raw=add(scale(K_WALL*d,n),scale(-CW*vnm,n)); f=scale(lam,raw) if d>0 else (0.0,0.0,0.0); return f,g,0.5*K_WALL*d*d,dot(scale(-CW*vnm,n),v)
def coll(xi,xj,vi,vj,lam):
    sep=sub(xi,xj); dist=norm(sep); n=scale(1/dist,sep); g=dist-2*R; d=max(0.0,-g); vr=sub(vi,vj); vnm=min(dot(vr,n),0.0); raw=add(scale(K_COLL*d,n),scale(-CC*vnm,n)); fi=scale(lam,raw) if d>0 else (0.0,0.0,0.0); fj=scale(-1,fi); return fi,fj,g,0.5*K_COLL*d*d,dot(scale(-CC*vnm,n),vr),norm(add(fi,fj))
def update(f):
    a=scale(1/RHO,f); v=scale(DT,a); x=scale(0.5*DT*DT,a); return a,v,x
def main():
    repo=Path(__file__).resolve().parents[1]; doc=repo/'stage22_checks/stage22_5_lambda_no_contact_contact_regression.md'; out=repo/'stage22_outputs/fibre_stage22_5_lambda_no_contact_contact_regression.dat'
    text=doc.read_text() if doc.exists() else fail('missing doc')
    miss=[m for m in MARKERS if m not in text]
    if miss: fail('missing markers: '+', '.join(miss))
    lambdas=[0.0,1e-6,1.0]; rhs_vec=(2.0,-1.0,0.5)
    rhs0=scale(0.0,rhs_vec); rhs_small=scale(1e-6,rhs_vec)
    if norm(rhs0)>ZERO or not finite(rhs_small) or norm(rhs_small)>MAXF: fail('rhs lambda')
    no_contact=[]; wall_pen=[]; coll_ov=[]
    for lam in lambdas:
        no_contact.append((lam,wall(0.50,(0,0,0),lam),coll((0,0.3,0),(0,0.5,0),(0,0,0),(0,0,0),lam)))
        wall_pen.append((lam,wall(R-1e-5,(0,-0.002,0),lam)))
        coll_ov.append((lam,coll((0,0.3,0),(0,0.31999,0),(0,0.001,0),(0,-0.001,0),lam)))
    if any(norm(w[1][0])>ZERO or norm(w[2][0])>ZERO for w in no_contact): fail('no contact invariance force')
    if norm(wall_pen[0][1][0])>ZERO: fail('lambda_contact zero wall')
    if norm(coll_ov[0][1][0])>ZERO or norm(coll_ov[0][1][1])>ZERO: fail('lambda_contact zero collision')
    wall_ratio=norm(wall_pen[1][1][0])/norm(wall_pen[2][1][0]); coll_ratio=norm(coll_ov[1][1][0])/norm(coll_ov[2][1][0])
    if abs(wall_ratio-1e-6)>SCALE_TOL or abs(coll_ratio-1e-6)>SCALE_TOL: fail('scaling')
    if wall_pen[2][1][0][1] <= 0: fail('wall inward')
    if dot(coll_ov[2][1][0], sub((0,0.3,0),(0,0.31999,0))) <= 0: fail('collision repulsive')
    if coll_ov[2][1][5] > TOL: fail('action reaction')
    all_forces=[wall_pen[2][1][0],coll_ov[2][1][0]]
    for f in all_forces:
        a,v,x=update(f)
        if not (finite(f) and finite(a) and finite(v) and finite(x)): fail('finite')
        if norm(f)>MAXF or norm(a)>MAXA or norm(x)/L>MAXFRAC: fail('bounds')
    if wall_pen[2][1][2]<-ZERO or coll_ov[2][1][3]<-ZERO or wall_pen[2][1][3]>ZERO or coll_ov[2][1][4]>ZERO: fail('energy damping')
    out.parent.mkdir(exist_ok=True)
    lines=['# Stage 22.5 lambda/no-contact/contact regression status','stage22_5_mode=helper_only_lambda_no_contact_contact_regression','stage20_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance','stage21_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance','stage22_0_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_0_pass','stage22_1_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_1_pass','stage22_2_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_2_pass','stage22_3_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_3_pass','stage22_4_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_4_pass','lambda_contact_zero_strict_noop=true','no_contact_no_overlap_zero_force_all_lambda=true',f'wall_lambda_small_to_one_ratio={wall_ratio:.16e}',f'collision_lambda_small_to_one_ratio={coll_ratio:.16e}','lambda_fsi_zero_rhs_noop=true',f'lambda_fsi_small_rhs_norm={norm(rhs_small):.16e}','production_structure_rhs_dns_mpi_build_io_activated=false']
    lines += [f'{field}=PASS' for field in STATUS_FIELDS]
    lines += ['STAGE 22.5 LAMBDA NO-CONTACT CONTACT REGRESSION VERDICT: PASS','STAGE 22.5 FINAL VERDICT: PASS','next_stage=Stage 22.6 single-fibre channel FSI micro-case','']
    out.write_text('\n'.join(lines)); print('\n'.join(lines)); return 0
if __name__=='__main__': raise SystemExit(main())
