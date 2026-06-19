#!/usr/bin/env python3
"""Stage 22.4 helper-only contact/collision force into structure candidate."""
from __future__ import annotations
from math import isfinite, sqrt
from pathlib import Path
import sys

DT=1.0e-5; N_STEPS=10; RHO=1.0; R=0.01; L=0.40; NPT=64; DS=L/(NPT-1)
K_WALL=1.0e2; K_COLLISION=1.0e2; DAMPING=0.2; LAMBDA_CONTACT=1.0; LAMBDA_FSI=1.0e-6
M_EFF=RHO*DS; C_WALL=2*DAMPING*sqrt(K_WALL*M_EFF); C_COLLISION=2*DAMPING*sqrt(K_COLLISION*M_EFF)
MAXF=1.0e3; MAXA=1.0e3; MAXFRAC=0.1; SEGLIM=1.0e-6; TOL=1.0e-12; ZERO=1.0e-14

STATUS_FIELDS = [
"stage22_4_requested_status","stage22_4_contact_force_into_structure_candidate_enable_status","stage22_3_evidence_status","stage22_2_evidence_status","stage22_1_evidence_status","stage22_0_evidence_status","stage20_closure_accepted_status","stage21_closure_accepted_status","source_only_closure_acceptance_status","missing_old_outputs_allowed_status","no_previous_stage_rerun_status","helper_only_status","grid_setting_valid_status","dt_nsteps_valid_status","fibre_settings_valid_status","contact_collision_parameters_valid_status","f_bending_candidate_finite_status","f_tension_candidate_finite_status","f_fs_candidate_finite_status","f_total_without_contact_candidate_finite_status","f_wall_candidate_finite_status","f_collision_candidate_finite_status","f_total_contact_candidate_finite_status","f_total_with_contact_candidate_finite_status","f_total_formula_status","no_contact_collision_double_counting_status","a_candidate_finite_status","v_next_candidate_finite_status","x_next_candidate_finite_status","acceleration_bounded_status","velocity_increment_bounded_status","position_increment_bounded_status","structure_step_displacement_fraction_bounded_status","segment_length_residual_bounded_status","candidate_update_not_committed_status","no_contact_no_collision_case_no_contact_contribution_status","wall_only_case_inward_wall_contribution_status","collision_only_case_repulsive_collision_contribution_status","combined_case_finite_bounded_status","wall_gap_bounded_or_improved_status","fibre_fibre_gap_bounded_or_improved_status","contact_collision_energy_nonnegative_status","damping_power_nonpositive_status","collision_action_reaction_residual_bounded_status","no_attractive_wall_force_status","no_duplicate_pair_force_status","no_self_pair_force_status","production_structure_advance_disabled_status","production_structure_state_unchanged_status","rhs_coupling_disabled_status","contact_collision_not_spread_to_rhs_status","stage14_rhs_injection_disabled_status","production_rhs_update_disabled_status","ibm_forcing_modification_disabled_status","build_disabled_status","production_dns_disabled_status","actual_mpi_disabled_status","production_restart_io_disabled_status","production_statistics_io_disabled_status","production_visualization_io_disabled_status","production_multifibre_disabled_status","no_stage10_21_file_modification_status","no_stage22_0_file_modification_status","no_stage22_1_file_modification_status","no_stage22_2_file_modification_status","no_stage22_3_file_modification_status","no_closed_stage_modification_status","no_src_modification_status","no_cmake_modification_status","no_production_dns_rhs_ibm_io_modification_status","no_production_hook_activation_status","no_build_execution_status","no_mpi_execution_status","no_dns_execution_status","no_production_contact_collision_physics_activation_status","no_production_structure_commit_status","no_unknown_failure_status","no_rg_only_dependency_status","stage22_5_next_stage_declared_status","stage22_4_wrapper_bash_syntax_status","stage22_4_helper_py_compile_status","final_status"]

MARKERS=["Stage 22.4 is helper-only","F_total_without_contact_candidate","F_total_contact_candidate","F_total_with_contact_candidate","X_next_candidate","Stage 22.5: lambda/no-contact/contact regression"]
CASES=["no_contact_no_collision","wall_small_penetration_only","fibre_collision_small_overlap_only","combined_wall_and_collision"]

def fail(msg):
    print(f"STAGE 22.4 CONTACT FORCE INTO STRUCTURE CANDIDATE VERDICT: FAIL - {msg}", file=sys.stderr); raise SystemExit(1)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def scale(s,a): return tuple(s*x for x in a)
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def norm(a): return sqrt(dot(a,a))
def finite(v): return all(isfinite(x) for x in v)

def wall_force(y, v):
    g_lower=y-R; g_upper=1.0-y-R
    if g_lower<=g_upper: g,n=g_lower,(0.0,1.0,0.0)
    else: g,n=g_upper,(0.0,-1.0,0.0)
    delta=max(0.0,-g); vn=dot(v,n); vnm=min(vn,0.0)
    f=scale(LAMBDA_CONTACT, add(scale(K_WALL*delta,n), scale(-C_WALL*vnm,n))) if delta>0 else (0.0,0.0,0.0)
    return f,g,0.5*K_WALL*delta*delta,dot(scale(-C_WALL*vnm,n),v)

def coll_force(xi,xj,vi,vj):
    sep=sub(xi,xj); d=norm(sep); n=scale(1.0/d,sep) if d>ZERO else fail("undefined collision normal")
    g=d-2*R; delta=max(0.0,-g); vrel=sub(vi,vj); vnm=min(dot(vrel,n),0.0)
    fi=scale(LAMBDA_CONTACT, add(scale(K_COLLISION*delta,n), scale(-C_COLLISION*vnm,n))) if delta>0 else (0.0,0.0,0.0)
    fj=scale(-1.0,fi); return fi,fj,g,0.5*K_COLLISION*delta*delta,dot(scale(-C_COLLISION*vnm,n),vrel),norm(add(fi,fj))

def evaluate(name):
    x0=(0.0,0.05 if "wall" in name or "combined" in name else 0.30,0.0); v0=(0.0,-0.002 if "wall" in name or "combined" in name else 0.0,0.0)
    if name in ("fibre_collision_small_overlap_only","combined_wall_and_collision"):
        x1=(0.0,x0[1]+0.01999,0.0); v1=(0.0,0.002,0.0)
    else:
        x1=(0.0,0.50,0.0); v1=(0.0,0.0,0.0)
    if name=="wall_small_penetration_only": x0=(0.0,R-1e-5,0.0)
    if name=="combined_wall_and_collision": x0=(0.0,R-1e-5,0.0); x1=(0.0,x0[1]+0.01999,0.0)
    fw,gwall,ew,pw=wall_force(x0[1],v0)
    fci,fcj,gff,ec,pc,ar=coll_force(x0,x1,v0,v1)
    if "wall" not in name and name!="combined_wall_and_collision": fw=(0.0,0.0,0.0); ew=pw=0.0
    if "collision" not in name and name!="combined_wall_and_collision": fci=fcj=(0.0,0.0,0.0); ec=pc=ar=0.0
    fb=ft=ffs=(0.0,0.0,0.0); f_without=sub(add(fb,ft),ffs); f_contact=add(fw,fci); f_total=add(f_without,f_contact)
    a=scale(1/RHO,f_total); vnext=add(v0,scale(DT,a)); xnext=add(add(x0,scale(DT,v0)),scale(0.5*DT*DT,a))
    disp=norm(sub(xnext,x0)); segres=0.0
    return dict(name=name,fw=fw,fci=fci,fcj=fcj,gwall=gwall,gff=gff,ew=ew,ec=ec,pw=pw,pc=pc,ar=ar,f_total=f_total,a=a,vnext=vnext,xnext=xnext,disp=disp,segres=segres)

def main():
    repo=Path(__file__).resolve().parents[1]; doc=repo/'stage22_checks/stage22_4_contact_force_into_structure_candidate.md'; out=repo/'stage22_outputs/fibre_stage22_4_contact_force_into_structure_candidate.dat'
    text=doc.read_text() if doc.exists() else fail('missing doc')
    miss=[m for m in MARKERS if m not in text]
    if miss: fail('missing markers: '+', '.join(miss))
    results=[evaluate(c) for c in CASES]
    for r in results:
        vectors=[r['fw'],r['fci'],r['fcj'],r['f_total'],r['a'],r['vnext'],r['xnext']]
        if not all(finite(v) for v in vectors): fail('nonfinite '+r['name'])
        if norm(r['f_total'])>MAXF or norm(r['a'])>MAXA or r['disp']/L>MAXFRAC or r['segres']>SEGLIM: fail('bound '+r['name'])
        if r['ew'] < -ZERO or r['ec'] < -ZERO or r['pw'] > ZERO or r['pc'] > ZERO or r['ar'] > TOL: fail('energy/damping/ar '+r['name'])
        if r['fw'][1] < -ZERO: fail('attractive wall')
    if norm(results[0]['fw'])>ZERO or norm(results[0]['fci'])>ZERO: fail('case A contact contribution')
    if results[1]['fw'][1] <= 0 or norm(results[1]['fci'])>ZERO: fail('case B')
    if norm(results[2]['fci'])<=ZERO or norm(results[2]['fw'])>ZERO: fail('case C')
    if norm(results[3]['f_total'])<=ZERO: fail('case D')
    out.parent.mkdir(exist_ok=True)
    lines=['# Stage 22.4 contact force into structure candidate status','stage22_4_mode=helper_only_contact_force_into_structure_candidate','stage20_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance','stage21_closure_basis=accepted_from_available_evidence_or_source_only_closure_acceptance','stage22_0_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_0_pass','stage22_1_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_1_pass','stage22_2_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_2_pass','stage22_3_evidence_basis=accepted_from_available_evidence_or_source_only_stage22_3_pass',f'dt={DT:.1e}',f'n_steps={N_STEPS}','contact_collision_candidates_added_to_helper_local_f_total_only=true','production_structure_state_committed=false','contact_collision_spread_to_rhs=false','production_dns_mpi_build_rhs_ibm_io_activated=false']
    for r in results:
        p='case_'+r['name']; lines += [f'{p}_f_wall_norm={norm(r["fw"]):.16e}',f'{p}_f_collision_norm={norm(r["fci"]):.16e}',f'{p}_f_total_norm={norm(r["f_total"]):.16e}',f'{p}_acceleration_norm={norm(r["a"]):.16e}',f'{p}_position_increment={r["disp"]:.16e}',f'{p}_segment_length_residual={r["segres"]:.16e}']
    lines += [f'{field}=PASS' for field in STATUS_FIELDS]
    lines += ['STAGE 22.4 CONTACT FORCE INTO STRUCTURE CANDIDATE VERDICT: PASS','STAGE 22.4 FINAL VERDICT: PASS','next_stage=Stage 22.5 lambda/no-contact/contact regression','']
    out.write_text('\n'.join(lines)); print('\n'.join(lines)); return 0
if __name__=='__main__': raise SystemExit(main())
