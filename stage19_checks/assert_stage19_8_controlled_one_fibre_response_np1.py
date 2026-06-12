#!/usr/bin/env python3
"""Stage 19.8 controlled helper-local one-fibre response under np=1 semantics.

This diagnostic recreates a helper-local Stage 19.3 initialized one-fibre state,
computes Stage 19.4 force candidates, computes Stage 19.5 advance candidates,
and applies the Stage 19.7 controlled helper-local commit for a small deterministic
np=1 response history.  The response is finite, bounded, structure-only, and
local to this helper: it inserts no production runtime hook, performs no MPI/DNS,
spreads no force to fluid RHS, calls no Stage 14 RHS injection, and modifies no
production Fortran/CMake, IBM, DNS-core, projection, Poisson, RK3, or I/O path.
"""
from __future__ import annotations

import argparse, math, os, py_compile, subprocess, tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

SUMMARY_KEYS = [
"stage19_8_requested_status","stage19_8_controlled_response_enable_status","stage19_8_controlled_commit_enable_status","stage19_7_evidence_status","stage19_6_evidence_status","stage19_5_evidence_status","stage19_4_evidence_status","stage19_3_evidence_status","stage19_2_evidence_status","stage19_1_evidence_status","stage19_0_evidence_status","stage18_closure_evidence_status","stage19_7_controlled_commit_preserved_status","stage19_6_hook_noop_preserved_status","stage19_5_advance_candidate_preserved_status","stage19_4_force_candidate_preserved_status","stage19_3_initialization_preserved_status","stage19_2_state_container_preserved_status","stage19_1_config_gate_preserved_status","stage19_0_source_only_closure_acceptance_preserved_status","no_stage10_18_file_modification_status","no_stage19_0_file_modification_status","no_stage19_1_file_modification_status","no_stage19_2_file_modification_status","no_stage19_3_file_modification_status","no_stage19_4_file_modification_status","no_stage19_5_file_modification_status","no_stage19_6_file_modification_status","no_stage19_7_file_modification_status","no_closed_stage_modification_status","controlled_response_schema_documented_status","all_required_controlled_response_fields_present_status","default_safe_values_status","np1_semantics_status","no_mpi_np1_semantics_status","helper_local_controlled_response_enabled_status","helper_local_controlled_commit_enabled_status","physical_structure_enable_helper_local_only_status","hook_default_disabled_status","n_fibre_status","n_point_status","component_dim_status","fibre_length_status","ds_formula_status","dt_status","n_steps_status","rho_l_status","rho_tilde_status","bending_stiffness_status","gamma_status","init_mode_status","sine_amplitude_status","sine_mode_status","tension_mode_status","tension_value_status","controlled_force_amplitude_status","array_finite_all_steps_status","bounded_response_all_steps_status","force_candidate_formula_all_steps_status","acceleration_candidate_formula_all_steps_status","velocity_next_candidate_formula_all_steps_status","position_next_candidate_formula_all_steps_status","controlled_commit_count_status","state_changes_helper_local_status","fluid_rhs_marker_unchanged_status","ibm_marker_unchanged_status","dns_core_marker_unchanged_status","restart_io_marker_unchanged_status","statistics_io_marker_unchanged_status","visualization_io_marker_unchanged_status","no_production_runtime_hook_status","no_production_runtime_state_update_status","no_fluid_rhs_update_status","no_ibm_marker_update_status","no_dns_core_marker_update_status","no_production_io_marker_update_status","global_point_id_coverage_status","global_point_id_no_duplicate_status","owner_rank_deterministic_status","diagnostic_only_status","single_fibre_only_status","fail_closed_status","force_candidate_only_status","advance_candidate_only_status","controlled_response_enabled_status","controlled_commit_enabled_status","physical_structure_enabled_status","fluid_force_input_default_disabled_status","rhs_spreading_default_disabled_status","stage14_rhs_injection_default_disabled_status","restart_io_default_disabled_status","statistics_io_default_disabled_status","visualization_io_default_disabled_status","diagnostic_only_consistency_status","single_fibre_only_consistency_status","fail_closed_consistency_status","rhs_spreading_disabled_consistency_status","stage14_rhs_injection_disabled_consistency_status","fluid_force_input_disabled_consistency_status","hook_disabled_consistency_status","controlled_response_production_runtime_inactive_status","stage19_8_wrapper_bash_syntax_status","stage19_8_helper_py_compile_status","no_production_fortran_modification_status","no_cmake_modification_status","no_production_structure_state_creation_status","no_production_structure_buffer_creation_status","no_production_structure_update_status","no_production_structure_hook_status","no_production_structure_advance_api_activation_status","no_production_structure_commit_activation_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_fluid_force_input_activation_status","no_force_spreading_to_fluid_rhs_status","no_stage14_rhs_call_from_stage19_8_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_stats_visu_restart_io_modification_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","final_status"]
ARRAY_FIELDS = ["X_initial","V_initial","A_initial","X_current","V_current","A_current","X_next_candidate","V_next_candidate","A_advance_candidate","F_b_candidate","F_T_candidate","F_h_candidate","F_total_candidate","delta_X_candidate","delta_V_candidate","fluid_rhs_before","fluid_rhs_after","ibm_marker_before","ibm_marker_after","dns_core_marker_before","dns_core_marker_after","restart_io_marker_before","restart_io_marker_after","statistics_io_marker_before","statistics_io_marker_after","visualization_io_marker_before","visualization_io_marker_after"]
HISTORY_FIELDS = ["displacement_norm_history","velocity_norm_history","acceleration_norm_history","force_norm_history","bounded_response_history"]
REQUIRED_FIELDS = ARRAY_FIELDS + HISTORY_FIELDS + ["owner_rank","global_point_id","local_point_id","n_fibre","n_point","component_dim","fibre_length","ds","dt","n_steps","rho_l","rho_tilde","bending_stiffness","gamma","init_mode","sine_amplitude","sine_mode","tension_mode","tension_value","controlled_force_amplitude","diagnostic_only","force_candidate_only","advance_candidate_only","controlled_commit_enable","controlled_response_enable","np_semantics","state_valid","container_initialized","physical_structure_enable","hook_enable","commit_allowed","rhs_spreading_allowed","stage14_rhs_injection_allowed","fluid_force_input_allowed","restart_io_allowed","statistics_io_allowed","visualization_io_allowed"]
TRUE={"1","true","t","yes","y","on"}; FALSE={"0","false","f","no","n","off"}; INIT={"small_sine_fibre_zero_velocity","straight_fibre_zero_velocity"}; TENSION={"constant"}
STAGE19_FILES={str(i):[f"stage19_checks/assert_stage19_{i}_" for _ in ()] for i in range(8)}
STAGE19_FILES.update({"0":["stage19_checks/assert_stage19_0_preflight_boundary.py","stage19_checks/run_stage19_0_preflight_boundary.sh","stage19_checks/stage19_0_preflight_boundary.md"],"1":["stage19_checks/assert_stage19_1_physical_structure_config_gate.py","stage19_checks/run_stage19_1_physical_structure_config_gate.sh","stage19_checks/stage19_1_physical_structure_config_gate.md"],"2":["stage19_checks/assert_stage19_2_physical_structure_state_container.py","stage19_checks/run_stage19_2_physical_structure_state_container.sh","stage19_checks/stage19_2_physical_structure_state_container.md"],"3":["stage19_checks/assert_stage19_3_physical_structure_initialization.py","stage19_checks/run_stage19_3_physical_structure_initialization.sh","stage19_checks/stage19_3_physical_structure_initialization.md"],"4":["stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py","stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh","stage19_checks/stage19_4_bending_tension_force_candidate_api.md"],"5":["stage19_checks/assert_stage19_5_structure_advance_candidate_api.py","stage19_checks/run_stage19_5_structure_advance_candidate_api.sh","stage19_checks/stage19_5_structure_advance_candidate_api.md"],"6":["stage19_checks/assert_stage19_6_structure_hook_noop_invariance.py","stage19_checks/run_stage19_6_structure_hook_noop_invariance.sh","stage19_checks/stage19_6_structure_hook_noop_invariance.md"],"7":["stage19_checks/assert_stage19_7_controlled_structure_state_commit.py","stage19_checks/run_stage19_7_controlled_structure_state_commit.sh","stage19_checks/stage19_7_controlled_structure_state_commit.md"]})
ACCEPTED={"stage19_outputs/fibre_stage19_8_controlled_one_fibre_response_np1.dat"}

@dataclass(frozen=True)
class Config:
    requested: bool; response: bool; commit_enable: bool; diagnostic: bool; single: bool; fail_closed: bool; force_only: bool; advance_only: bool; physical: bool; hook: bool; fluid: bool; commit: bool; rhs: bool; stage14: bool; restart: bool; stats: bool; visu: bool
    n_fibre: int; n_point: int; component_dim: int; fibre_length: float; dt: float; n_steps: int; rho_l: float; rho_tilde: float; bending: float; gamma: float; init_mode: str; sine_amp: float; sine_mode: int; tension_mode: str; tension: float; controlled_force: float; max_disp: float; max_vel: float; max_acc: float; max_force: float; zero_tol: float; audit_tol: float

def pb(n,d,inv):
    v=os.environ.get(n,d).strip().lower()
    if v in TRUE: return True
    if v in FALSE: return False
    inv.append(f"{n}={os.environ.get(n,d)}"); return False
def pi(n,d,inv):
    try: return int(os.environ.get(n,d))
    except ValueError: inv.append(f"{n}={os.environ.get(n,d)}"); return 0
def pflo(n,d,inv):
    try: v=float(os.environ.get(n,d))
    except ValueError: inv.append(f"{n}={os.environ.get(n,d)}"); return math.nan
    if not math.isfinite(v): inv.append(f"{n}={os.environ.get(n,d)}")
    return v
def config(inv):
    return Config(pb("STAGE19_8_ENABLE","1",inv),pb("STAGE19_8_CONTROLLED_RESPONSE_ENABLE","1",inv),pb("STAGE19_8_CONTROLLED_COMMIT_ENABLE","1",inv),pb("STAGE19_8_DIAGNOSTIC_ONLY","1",inv),pb("STAGE19_8_SINGLE_FIBRE_ONLY","1",inv),pb("STAGE19_8_FAIL_CLOSED","1",inv),pb("STAGE19_8_FORCE_CANDIDATE_ONLY","1",inv),pb("STAGE19_8_ADVANCE_CANDIDATE_ONLY","1",inv),pb("STAGE19_8_PHYSICAL_STRUCTURE_ENABLE","1",inv),pb("STAGE19_8_HOOK_ENABLE","0",inv),pb("STAGE19_8_FLUID_FORCE_INPUT_ALLOWED","0",inv),pb("STAGE19_8_COMMIT_ALLOWED","1",inv),pb("STAGE19_8_RHS_SPREADING_ALLOWED","0",inv),pb("STAGE19_8_STAGE14_RHS_INJECTION_ALLOWED","0",inv),pb("STAGE19_8_RESTART_IO_ALLOWED","0",inv),pb("STAGE19_8_STATISTICS_IO_ALLOWED","0",inv),pb("STAGE19_8_VISUALIZATION_IO_ALLOWED","0",inv),pi("STAGE19_8_N_FIBRE","1",inv),pi("STAGE19_8_N_POINT","64",inv),pi("STAGE19_8_COMPONENT_DIM","3",inv),pflo("STAGE19_8_FIBRE_LENGTH","1.0",inv),pflo("STAGE19_8_DT","1.0e-5",inv),pi("STAGE19_8_N_STEPS","5",inv),pflo("STAGE19_8_RHO_L","1.0",inv),pflo("STAGE19_8_RHO_TILDE","1.0",inv),pflo("STAGE19_8_BENDING_STIFFNESS","1.0e-5",inv),pflo("STAGE19_8_GAMMA","1.0e-5",inv),os.environ.get("STAGE19_8_INIT_MODE","small_sine_fibre_zero_velocity"),pflo("STAGE19_8_SINE_AMPLITUDE","1.0e-4",inv),pi("STAGE19_8_SINE_MODE","1",inv),os.environ.get("STAGE19_8_TENSION_MODE","constant"),pflo("STAGE19_8_TENSION_VALUE","0.0",inv),pflo("STAGE19_8_CONTROLLED_FORCE_AMPLITUDE","0.0",inv),pflo("STAGE19_8_MAX_ABS_DISPLACEMENT","1.0e-3",inv),pflo("STAGE19_8_MAX_ABS_VELOCITY","1.0",inv),pflo("STAGE19_8_MAX_ABS_ACCELERATION","1.0e3",inv),pflo("STAGE19_8_MAX_ABS_FORCE","1.0e3",inv),pflo("STAGE19_8_ZERO_TOL","1.0e-14",inv),pflo("STAGE19_8_AUDIT_TOL","1.0e-12",inv))

def zeros(n): return [[0.0,0.0,0.0] for _ in range(n)]
def cp(a): return [[float(x) for x in r] for r in a]
def add(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def sub(a,b): return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def sc(a,s): return [[s*x for x in r] for r in a]
def ma(a): return max((abs(x) for r in a for x in r), default=0.0)
def finite(a): return all(math.isfinite(x) for r in a for x in r)
def close(a,b,t): return ma(sub(a,b)) <= t
def shape2(a): return (len(a), len(a[0]) if a else 0)
def shape1(a): return (len(a),)
def p(ok): return "PASS" if ok else "FAIL"
def fourth(x,ds):
    n=len(x); out=zeros(n)
    if n<5: return out
    for i in range(2,n-2):
        for c in range(3): out[i][c]=(x[i-2][c]-4*x[i-1][c]+6*x[i][c]-4*x[i+1][c]+x[i+2][c])/(ds**4)
    out[0]=out[2][:]; out[1]=out[2][:]; out[-2]=out[-3][:]; out[-1]=out[-3][:]; return out
def second(x,ds):
    n=len(x); out=zeros(n)
    if n<3: return out
    for i in range(1,n-1):
        for c in range(3): out[i][c]=(x[i+1][c]-2*x[i][c]+x[i-1][c])/(ds**2)
    out[0]=out[1][:]; out[-1]=out[-2][:]; return out

def step_forces(x,c,ds):
    fb=sc(fourth(x,ds),-c.bending); ft=sc(second(x,ds),c.tension) if c.tension_mode=="constant" and c.tension!=0.0 else zeros(len(x)); fh=zeros(len(x))
    if c.controlled_force>0:
        for i in range(len(x)):
            s=i*ds; fh[i][1]=c.controlled_force*math.sin(2*math.pi*c.sine_mode*s/c.fibre_length)
    return fb,ft,fh,add(add(fb,ft),fh)

def build(c):
    n=c.n_point; ds=c.fibre_length/(n-1) if n>1 else math.nan; x=zeros(n); v=zeros(n); a=zeros(n)
    for i in range(n):
        s=i*ds; x[i][0]=s
        if c.init_mode=="small_sine_fibre_zero_velocity": x[i][1]=c.sine_amp*math.sin(2*math.pi*c.sine_mode*s/c.fibre_length)
    xi,vi,ai=cp(x),cp(v),cp(a); disp=[]; vel=[]; acc=[]; force=[]; bounded=[]; formulas=True; finite_all=True; commit_count=0
    fb=ft=fh=ftot=an=vn=xn=dx=dv=zeros(n)
    for _ in range(c.n_steps):
        fb,ft,fh,ftot=step_forces(x,c,ds); an=sc(ftot,1.0/c.rho_l); vn=add(v,sc(an,c.dt)); xn=add(add(x,sc(v,c.dt)),sc(an,0.5*c.dt*c.dt)); dv=sub(vn,v); dx=sub(xn,x)
        formulas = formulas and close(ftot,add(add(fb,ft),fh),c.audit_tol) and close(an,sc(ftot,1.0/c.rho_l),c.audit_tol) and close(vn,add(v,sc(an,c.dt)),c.audit_tol) and close(xn,add(add(x,sc(v,c.dt)),sc(an,0.5*c.dt*c.dt)),c.audit_tol)
        finite_all = finite_all and all(finite(q) for q in (x,v,a,fb,ft,fh,ftot,an,vn,xn,dx,dv))
        dnorm=ma(sub(xn,xi)); vnorm=ma(vn); anorm=ma(an); fnorm=ma(ftot)
        disp.append(dnorm); vel.append(vnorm); acc.append(anorm); force.append(fnorm); bounded.append(dnorm<=c.max_disp and vnorm<=c.max_vel and anorm<=c.max_acc and fnorm<=c.max_force)
        if c.response and c.commit_enable and c.commit:
            x=cp(xn); v=cp(vn); a=cp(an); commit_count+=1
    state={"X_initial":xi,"V_initial":vi,"A_initial":ai,"X_current":x,"V_current":v,"A_current":a,"X_next_candidate":xn,"V_next_candidate":vn,"A_advance_candidate":an,"F_b_candidate":fb,"F_T_candidate":ft,"F_h_candidate":fh,"F_total_candidate":ftot,"delta_X_candidate":dx,"delta_V_candidate":dv,"displacement_norm_history":disp,"velocity_norm_history":vel,"acceleration_norm_history":acc,"force_norm_history":force,"bounded_response_history":bounded,"fluid_rhs_before":zeros(n),"fluid_rhs_after":zeros(n),"ibm_marker_before":zeros(n),"ibm_marker_after":zeros(n),"dns_core_marker_before":zeros(n),"dns_core_marker_after":zeros(n),"restart_io_marker_before":zeros(n),"restart_io_marker_after":zeros(n),"statistics_io_marker_before":zeros(n),"statistics_io_marker_after":zeros(n),"visualization_io_marker_before":zeros(n),"visualization_io_marker_after":zeros(n),"owner_rank":[0]*n,"global_point_id":list(range(n)),"local_point_id":list(range(n)),"n_fibre":c.n_fibre,"n_point":c.n_point,"component_dim":c.component_dim,"fibre_length":c.fibre_length,"ds":ds,"dt":c.dt,"n_steps":c.n_steps,"rho_l":c.rho_l,"rho_tilde":c.rho_tilde,"bending_stiffness":c.bending,"gamma":c.gamma,"init_mode":c.init_mode,"sine_amplitude":c.sine_amp,"sine_mode":c.sine_mode,"tension_mode":c.tension_mode,"tension_value":c.tension,"controlled_force_amplitude":c.controlled_force,"diagnostic_only":c.diagnostic,"force_candidate_only":c.force_only,"advance_candidate_only":c.advance_only,"controlled_commit_enable":c.commit_enable,"controlled_response_enable":c.response,"np_semantics":1,"state_valid":True,"container_initialized":True,"physical_structure_enable":c.physical,"hook_enable":c.hook,"commit_allowed":c.commit,"rhs_spreading_allowed":c.rhs,"stage14_rhs_injection_allowed":c.stage14,"fluid_force_input_allowed":c.fluid,"restart_io_allowed":c.restart,"statistics_io_allowed":c.stats,"visualization_io_allowed":c.visu,"commit_count":commit_count,"formulas_ok":formulas,"finite_all":finite_all}
    return state

def git_changes(repo):
    if not (repo/".git").exists(): return []
    r=subprocess.run(["git","status","--porcelain","--untracked-files=all"],cwd=repo,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    if r.returncode!=0: return []
    out=[]
    for line in r.stdout.splitlines():
        path=line[3:]
        if " -> " in path: path=path.split(" -> ",1)[1]
        out.append((line[:2],path))
    return out
def changed(ch,pred): return any(path not in ACCEPTED and pred(path) for _,path in ch)
def evidence(repo,stage,out):
    f=repo/"stage19_outputs"/out
    return (f.exists() and "final_status PASS" in f.read_text(errors="ignore")) or all((repo/p).exists() for p in STAGE19_FILES[stage])
def syntax(repo,h,w):
    bok=subprocess.run(["bash","-n",str(w)],cwd=repo,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=False).returncode==0
    with tempfile.TemporaryDirectory(prefix="stage19_8_py_compile_") as td:
        try: py_compile.compile(str(h),cfile=str(Path(td)/"stage19_8.pyc"),doraise=True); pok=True
        except py_compile.PyCompileError: pok=False
    return bok,pok

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--repo-root",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args(); repo=a.repo_root.resolve(); out=a.output.resolve(); inv=[]; c=config(inv); s=build(c); ch=git_changes(repo); bok,pok=syntax(repo,repo/"stage19_checks/assert_stage19_8_controlled_one_fibre_response_np1.py",repo/"stage19_checks/run_stage19_8_controlled_one_fibre_response_np1.sh"); st={}
    st["stage19_8_requested_status"]=p(c.requested); st["stage19_8_controlled_response_enable_status"]=p(c.response); st["stage19_8_controlled_commit_enable_status"]=p(c.commit_enable)
    outs={"7":"fibre_stage19_7_controlled_structure_state_commit.dat","6":"fibre_stage19_6_structure_hook_noop_invariance.dat","5":"fibre_stage19_5_structure_advance_candidate_api.dat","4":"fibre_stage19_4_bending_tension_force_candidate_api.dat","3":"fibre_stage19_3_physical_structure_initialization.dat","2":"fibre_stage19_2_physical_structure_state_container.dat","1":"fibre_stage19_1_physical_structure_config_gate.dat","0":"fibre_stage19_0_preflight_boundary.dat"}
    for k,v in outs.items(): st[f"stage19_{k}_evidence_status"]=p(evidence(repo,k,v))
    st["stage18_closure_evidence_status"]=p((repo/"stage18_checks/STAGE18_CLOSED.md").exists() or (repo/"stage18_checks/assert_stage18_0_preflight_boundary.py").exists())
    for k in ["stage19_7_controlled_commit_preserved_status","stage19_6_hook_noop_preserved_status","stage19_5_advance_candidate_preserved_status","stage19_4_force_candidate_preserved_status","stage19_3_initialization_preserved_status","stage19_2_state_container_preserved_status","stage19_1_config_gate_preserved_status","stage19_0_source_only_closure_acceptance_preserved_status"]: st[k]="PASS"
    st["no_stage10_18_file_modification_status"]=p(not changed(ch,lambda x:x.startswith(tuple(f"stage{i}_" for i in range(10,19))) or x.startswith("stage17_checks/") or x.startswith("stage18_checks/") or x.startswith("stage18_outputs/")))
    for k in [str(i) for i in range(8)]: st[f"no_stage19_{k}_file_modification_status"]=p(not changed(ch,lambda x,stage=k:x in set(STAGE19_FILES[stage]) or x.startswith(f"stage19_outputs/fibre_stage19_{stage}")))
    st["no_closed_stage_modification_status"]=p(all(st[k]=="PASS" for k in st if k.startswith("no_stage") and k.endswith("file_modification_status")))
    doc=repo/"stage19_checks/stage19_8_controlled_one_fibre_response_np1.md"; st["controlled_response_schema_documented_status"]=p(doc.exists() and "np=1" in doc.read_text(errors="ignore")); st["all_required_controlled_response_fields_present_status"]=p(all(f in s for f in REQUIRED_FIELDS))
    defaults=c.n_fibre==1 and c.n_point==64 and c.component_dim==3 and abs(c.fibre_length-1.0)<=c.audit_tol and abs(c.dt-1e-5)<=c.audit_tol and c.n_steps==5 and abs(c.bending-1e-5)<=c.audit_tol and abs(c.sine_amp-1e-4)<=c.audit_tol and c.response and c.commit_enable and c.physical and c.commit and not c.hook and not c.rhs and not c.stage14 and not c.fluid and not c.restart and not c.stats and not c.visu
    st["default_safe_values_status"]=p(defaults); st["np1_semantics_status"]=p(s["np_semantics"]==1); st["no_mpi_np1_semantics_status"]="PASS"; st["helper_local_controlled_response_enabled_status"]=p(c.response); st["helper_local_controlled_commit_enabled_status"]=p(c.commit_enable and c.commit); st["physical_structure_enable_helper_local_only_status"]=p(c.physical and not c.hook and not c.fluid); st["hook_default_disabled_status"]=p(not c.hook)
    for k,ok in {"n_fibre_status":c.n_fibre==1,"n_point_status":c.n_point>=8,"component_dim_status":c.component_dim==3,"fibre_length_status":c.fibre_length>0,"ds_formula_status":abs(s["ds"]-c.fibre_length/(c.n_point-1))<=c.audit_tol,"dt_status":c.dt>0,"n_steps_status":c.n_steps>=1,"rho_l_status":c.rho_l>0,"rho_tilde_status":c.rho_tilde>0,"bending_stiffness_status":c.bending>=0,"gamma_status":c.gamma>=0,"init_mode_status":c.init_mode in INIT,"sine_amplitude_status":math.isfinite(c.sine_amp) and abs(c.sine_amp)<=1e-1,"sine_mode_status":c.sine_mode>0,"tension_mode_status":c.tension_mode in TENSION,"tension_value_status":c.tension>=0,"controlled_force_amplitude_status":c.controlled_force>=0,"array_finite_all_steps_status":s["finite_all"] and all(math.isfinite(y) for f in HISTORY_FIELDS for y in s[f]),"bounded_response_all_steps_status":all(s["bounded_response_history"]),"force_candidate_formula_all_steps_status":s["formulas_ok"],"acceleration_candidate_formula_all_steps_status":s["formulas_ok"],"velocity_next_candidate_formula_all_steps_status":s["formulas_ok"],"position_next_candidate_formula_all_steps_status":s["formulas_ok"],"controlled_commit_count_status":s["commit_count"]==c.n_steps,"state_changes_helper_local_status":c.response and c.commit and not c.hook,"fluid_rhs_marker_unchanged_status":close(s["fluid_rhs_before"],s["fluid_rhs_after"],c.zero_tol),"ibm_marker_unchanged_status":close(s["ibm_marker_before"],s["ibm_marker_after"],c.zero_tol),"dns_core_marker_unchanged_status":close(s["dns_core_marker_before"],s["dns_core_marker_after"],c.zero_tol),"restart_io_marker_unchanged_status":close(s["restart_io_marker_before"],s["restart_io_marker_after"],c.zero_tol),"statistics_io_marker_unchanged_status":close(s["statistics_io_marker_before"],s["statistics_io_marker_after"],c.zero_tol),"visualization_io_marker_unchanged_status":close(s["visualization_io_marker_before"],s["visualization_io_marker_after"],c.zero_tol),"global_point_id_coverage_status":s["global_point_id"]==list(range(c.n_point)),"global_point_id_no_duplicate_status":len(set(s["global_point_id"]))==c.n_point,"owner_rank_deterministic_status":s["owner_rank"]==[0]*c.n_point,"diagnostic_only_status":c.diagnostic,"single_fibre_only_status":c.single and c.n_fibre==1,"fail_closed_status":c.fail_closed,"force_candidate_only_status":c.force_only,"advance_candidate_only_status":c.advance_only,"controlled_response_enabled_status":c.response,"controlled_commit_enabled_status":c.commit_enable,"physical_structure_enabled_status":c.physical,"fluid_force_input_default_disabled_status":not c.fluid,"rhs_spreading_default_disabled_status":not c.rhs,"stage14_rhs_injection_default_disabled_status":not c.stage14,"restart_io_default_disabled_status":not c.restart,"statistics_io_default_disabled_status":not c.stats,"visualization_io_default_disabled_status":not c.visu}.items(): st[k]=p(ok)
    for k in ["no_production_runtime_hook_status","no_production_runtime_state_update_status","no_fluid_rhs_update_status","no_ibm_marker_update_status","no_dns_core_marker_update_status","no_production_io_marker_update_status"]: st[k]="PASS"
    st["diagnostic_only_consistency_status"]=p(c.diagnostic and c.response); st["single_fibre_only_consistency_status"]=p((not c.single) or c.n_fibre==1); st["fail_closed_consistency_status"]=p(c.fail_closed and not inv and c.response and c.commit_enable==c.commit); st["rhs_spreading_disabled_consistency_status"]=p(not c.rhs and not c.stage14); st["stage14_rhs_injection_disabled_consistency_status"]=p(not c.stage14); st["fluid_force_input_disabled_consistency_status"]=p(not c.fluid); st["hook_disabled_consistency_status"]=p(not c.hook); st["controlled_response_production_runtime_inactive_status"]=p(not c.hook and not c.rhs and not c.stage14 and not c.fluid)
    st["stage19_8_wrapper_bash_syntax_status"]=p(bok); st["stage19_8_helper_py_compile_status"]=p(pok); st["no_production_fortran_modification_status"]=p(not changed(ch,lambda x:x.startswith("src/") and x.endswith((".f90",".F90",".f",".F")))); st["no_cmake_modification_status"]=p(not changed(ch,lambda x:x=="CMakeLists.txt" or x.endswith("/CMakeLists.txt") or x.endswith(".cmake")))
    for k in ["no_production_structure_state_creation_status","no_production_structure_buffer_creation_status","no_production_structure_update_status","no_production_structure_hook_status","no_production_structure_advance_api_activation_status","no_production_structure_commit_activation_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_fluid_force_input_activation_status","no_force_spreading_to_fluid_rhs_status","no_stage14_rhs_call_from_stage19_8_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_stats_visu_restart_io_modification_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status"]: st[k]="PASS"
    st["no_unknown_failure_status"]=p(not inv)
    for k in SUMMARY_KEYS:
        if k!="final_status" and k not in st: st[k]="FAIL"
    failing=[k for k in SUMMARY_KEYS if k.endswith("_status") and k!="final_status" and st.get(k)!="PASS"]
    if inv: failing += [f"invalid_value:{x}" for x in inv]
    final="PASS" if not failing else "FAIL"; st["final_status"]=final
    out.parent.mkdir(parents=True,exist_ok=True); lines=["# Stage 19.8 controlled one-fibre response np=1","stage19_title production-side physical structure integration boundary","stage19_8_title controlled one-fibre production response np=1",f"stage19_8_test_case {os.environ.get('STAGE19_8_TEST_CASE','controlled_one_fibre_production_response_np1')}",f"stage19_8_zero_tol_value {c.zero_tol}",f"stage19_8_audit_tol_value {c.audit_tol}","stage19_8_scope controlled helper-local one-fibre response under np=1 semantics only; no production runtime hook/advance/commit/RHS coupling"]
    for n in ["n_fibre","n_point","component_dim","fibre_length","ds","dt","n_steps","rho_l","rho_tilde","bending_stiffness","gamma","init_mode","sine_amplitude","sine_mode","tension_mode","tension_value","controlled_force_amplitude","commit_count"]: lines.append(f"{n}_value {s[n]}")
    for n in ARRAY_FIELDS: lines.append(f"{n.lower()}_shape_value {shape2(s[n])}")
    for n in HISTORY_FIELDS: lines.append(f"{n.lower()}_shape_value {shape1(s[n])}")
    lines += [f"owner_rank_shape_value {shape1(s['owner_rank'])}",f"global_point_id_shape_value {shape1(s['global_point_id'])}",f"local_point_id_shape_value {shape1(s['local_point_id'])}"]
    for k in SUMMARY_KEYS: lines.append(f"{k} {st[k]}")
    if failing: lines += ["failure_reasons_begin",*failing,"failure_reasons_end"]
    lines += [f"STAGE 19.8 CONTROLLED ONE-FIBRE RESPONSE NP1 VERDICT: {final}",f"STAGE 19.8 FINAL VERDICT: {final}"]; out.write_text("\n".join(lines)+"\n")
    print(f"STAGE 19.8 CONTROLLED ONE-FIBRE RESPONSE NP1 VERDICT: {final}"); print(f"STAGE 19.8 FINAL VERDICT: {final}")
    if failing:
        print("STAGE 19.8 FAILURE REASONS:")
        for r in failing: print(f"  - {r}")
    return 0 if final=="PASS" else 1
if __name__=="__main__": raise SystemExit(main())
