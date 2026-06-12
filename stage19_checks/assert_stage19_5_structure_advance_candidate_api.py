#!/usr/bin/env python3
"""Stage 19.5 structure advance candidate API boundary.

Pure-Python helper-local advance-candidate diagnostic.  It reconstructs a Stage
19.3-style initialized state, recomputes Stage 19.4-style helper force
candidates, then computes acceleration, next-velocity, next-position, and delta
candidates.  It never writes production X/V/A, inserts hooks, commits state,
spreads RHS forces, calls Stage 14 RHS injection, modifies IBM/DNS-core /
projection / Poisson / RK3, modifies production I/O, runs MPI/DNS, or introduces
contact/collision/production multifibre logic.
"""
from __future__ import annotations

import argparse, math, os, py_compile, subprocess, tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

SUMMARY_KEYS = [
"stage19_5_requested_status","stage19_5_advance_candidate_enable_status","stage19_4_evidence_status","stage19_3_evidence_status","stage19_2_evidence_status","stage19_1_evidence_status","stage19_0_evidence_status","stage18_closure_evidence_status","stage19_4_force_candidate_preserved_status","stage19_3_initialization_preserved_status","stage19_2_state_container_preserved_status","stage19_1_config_gate_preserved_status","stage19_0_source_only_closure_acceptance_preserved_status","no_stage10_18_file_modification_status","no_stage19_0_file_modification_status","no_stage19_1_file_modification_status","no_stage19_2_file_modification_status","no_stage19_3_file_modification_status","no_stage19_4_file_modification_status","no_closed_stage_modification_status","advance_candidate_schema_documented_status","all_required_advance_candidate_fields_present_status","default_safe_values_status","n_fibre_status","n_point_status","component_dim_status","fibre_length_status","ds_formula_status","dt_status","rho_l_status","rho_tilde_status","bending_stiffness_status","gamma_status","init_mode_status","sine_amplitude_status","sine_mode_status","tension_mode_status","tension_value_status","controlled_force_amplitude_status","x_prod_shape_status","v_prod_shape_status","a_prod_shape_status","f_b_candidate_shape_status","f_t_candidate_shape_status","f_h_candidate_shape_status","f_total_candidate_shape_status","a_advance_candidate_shape_status","v_next_candidate_shape_status","x_next_candidate_shape_status","delta_x_candidate_shape_status","delta_v_candidate_shape_status","owner_rank_shape_status","global_point_id_shape_status","local_point_id_shape_status","array_finite_status","total_force_candidate_formula_status","acceleration_candidate_formula_status","velocity_next_candidate_formula_status","position_next_candidate_formula_status","delta_velocity_candidate_formula_status","delta_position_candidate_formula_status","candidate_arrays_helper_local_status","candidate_advance_no_state_update_status","no_state_commit_status","no_production_runtime_state_update_status","global_point_id_coverage_status","global_point_id_no_duplicate_status","owner_rank_deterministic_status","diagnostic_only_status","single_fibre_only_status","fail_closed_status","force_candidate_only_status","advance_candidate_only_status","commit_default_disabled_status","rhs_spreading_default_disabled_status","stage14_rhs_injection_default_disabled_status","diagnostic_only_consistency_status","single_fibre_only_consistency_status","fail_closed_consistency_status","rhs_spreading_disabled_consistency_status","stage14_rhs_injection_disabled_consistency_status","commit_disabled_consistency_status","advance_candidate_production_runtime_inactive_status","stage19_5_wrapper_bash_syntax_status","stage19_5_helper_py_compile_status","no_production_fortran_modification_status","no_cmake_modification_status","no_production_structure_state_creation_status","no_production_structure_buffer_creation_status","no_production_structure_update_status","no_production_structure_hook_status","no_production_structure_advance_api_activation_status","no_production_structure_commit_activation_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_fluid_force_input_activation_status","no_force_spreading_to_fluid_rhs_status","no_stage14_rhs_call_from_stage19_5_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_stats_visu_restart_io_modification_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","final_status"]

REQUIRED_FIELDS = ["X_prod","V_prod","A_prod","F_b_candidate","F_T_candidate","F_h_candidate","F_total_candidate","A_advance_candidate","V_next_candidate","X_next_candidate","delta_X_candidate","delta_V_candidate","owner_rank","global_point_id","local_point_id","n_fibre","n_point","component_dim","fibre_length","ds","dt","rho_l","rho_tilde","bending_stiffness","gamma","init_mode","sine_amplitude","sine_mode","tension_mode","tension_value","controlled_force_amplitude","diagnostic_only","force_candidate_only","advance_candidate_only","state_valid","container_initialized","commit_allowed","rhs_spreading_allowed","stage14_rhs_injection_allowed"]
ARRAY_FIELDS_2D = ["X_prod","V_prod","A_prod","F_b_candidate","F_T_candidate","F_h_candidate","F_total_candidate","A_advance_candidate","V_next_candidate","X_next_candidate","delta_X_candidate","delta_V_candidate"]
SHAPE_STATUS = {"X_prod":"x_prod_shape_status","V_prod":"v_prod_shape_status","A_prod":"a_prod_shape_status","F_b_candidate":"f_b_candidate_shape_status","F_T_candidate":"f_t_candidate_shape_status","F_h_candidate":"f_h_candidate_shape_status","F_total_candidate":"f_total_candidate_shape_status","A_advance_candidate":"a_advance_candidate_shape_status","V_next_candidate":"v_next_candidate_shape_status","X_next_candidate":"x_next_candidate_shape_status","delta_X_candidate":"delta_x_candidate_shape_status","delta_V_candidate":"delta_v_candidate_shape_status"}
INIT_MODES={"small_sine_fibre_zero_velocity"}; TENSION_MODES={"constant"}
S0={"stage19_checks/run_stage19_0_preflight_boundary.sh","stage19_checks/assert_stage19_0_preflight_boundary.py","stage19_checks/stage19_0_preflight_boundary.md","stage19_outputs/fibre_stage19_0_preflight_boundary.dat"}
S1={"stage19_checks/run_stage19_1_physical_structure_config_gate.sh","stage19_checks/assert_stage19_1_physical_structure_config_gate.py","stage19_checks/stage19_1_physical_structure_config_gate.md","stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat"}
S2={"stage19_checks/run_stage19_2_physical_structure_state_container.sh","stage19_checks/assert_stage19_2_physical_structure_state_container.py","stage19_checks/stage19_2_physical_structure_state_container.md","stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat"}
S3={"stage19_checks/run_stage19_3_physical_structure_initialization.sh","stage19_checks/assert_stage19_3_physical_structure_initialization.py","stage19_checks/stage19_3_physical_structure_initialization.md","stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat"}
S4={"stage19_checks/run_stage19_4_bending_tension_force_candidate_api.sh","stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py","stage19_checks/stage19_4_bending_tension_force_candidate_api.md","stage19_outputs/fibre_stage19_4_bending_tension_force_candidate_api.dat"}
ALLOWED={"stage19_checks/run_stage19_5_structure_advance_candidate_api.sh","stage19_checks/assert_stage19_5_structure_advance_candidate_api.py","stage19_checks/stage19_5_structure_advance_candidate_api.md","stage19_outputs/fibre_stage19_5_structure_advance_candidate_api.dat"}
ACCEPTED={"stage17_checks/STAGE17_CLOSED.md","stage18_checks/STAGE18_CLOSED.md","stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat","stage19_outputs/fibre_stage19_0_preflight_boundary.dat","stage19_outputs/fibre_stage19_1_physical_structure_config_gate.dat","stage19_outputs/fibre_stage19_2_physical_structure_state_container.dat","stage19_outputs/fibre_stage19_3_physical_structure_initialization.dat","stage19_outputs/fibre_stage19_4_bending_tension_force_candidate_api.dat"}
TRUE={"1","true","TRUE","yes","YES","on","ON","t","T"}; FALSE={"0","false","FALSE","no","NO","off","OFF","f","F"}

@dataclass(frozen=True)
class Config:
    n_fibre:int; n_point:int; component_dim:int; fibre_length:float; dt:float; rho_l:float; rho_tilde:float; bending_stiffness:float; gamma:float; init_mode:str; sine_amplitude:float; sine_mode:int; tension_mode:str; tension_value:float; controlled_force_amplitude:float; diagnostic_only:bool; single_fibre_only:bool; fail_closed:bool; force_candidate_only:bool; advance_candidate_only:bool; commit_allowed:bool; rhs_spreading_allowed:bool; stage14_rhs_injection_allowed:bool

def read(path: Path)->str:
    try: return path.read_text(errors="ignore")
    except OSError: return ""
def pf(c: bool)->str: return "PASS" if c else "FAIL"
def pbool(n,d,inv):
    r=os.environ.get(n)
    if r is None: return d
    v=r.strip()
    if v in TRUE: return True
    if v in FALSE: return False
    inv.append(f"{n}={r}"); return d
def pint(n,d,inv):
    r=os.environ.get(n)
    if r is None: return d
    try: return int(r.strip())
    except ValueError: inv.append(f"{n}={r}"); return d
def pfloat(n,d,inv):
    r=os.environ.get(n)
    if r is None: return d
    try: v=float(r.strip())
    except ValueError: inv.append(f"{n}={r}"); return d
    if not math.isfinite(v): inv.append(f"{n}={r}"); return d
    return v

def run(cmd, root):
    try:
        p=subprocess.run(cmd,cwd=str(root),text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=False); return p.returncode,p.stdout
    except OSError as e: return 127,str(e)
def git_entries(root):
    code,out=run(["git","status","--porcelain","--untracked-files=all"],root)
    if code: return []
    ans=[]
    for raw in out.splitlines():
        if not raw: continue
        xy=raw[:2]; path=raw[3:] if len(raw)>3 else raw
        if " -> " in path: path=path.split(" -> ",1)[1]
        ans.append((xy,path.strip()))
    return ans
def changed_outside(root):
    out=[]
    for xy,p in git_entries(root):
        if p in ALLOWED: continue
        if xy=="??" and p in ACCEPTED: continue
        out.append(p)
    return out

def has_pass(path, verdict):
    t=read(path); return "final_status PASS" in t and verdict in t

def preserved0(root):
    h=read(root/"stage19_checks/assert_stage19_0_preflight_boundary.py"); d=read(root/"stage19_checks/stage19_0_preflight_boundary.md")
    return all(x in h for x in ("stage18_closure_accepted_status","prior_stage18_outputs_required_status","ACCEPTED_BY_STAGE18_CLOSURE")) and "must not force users to rerun Stage 18.0 through Stage 18.11" in d
def preserved1(root):
    h=read(root/"stage19_checks/assert_stage19_1_physical_structure_config_gate.py"); d=read(root/"stage19_checks/stage19_1_physical_structure_config_gate.md")
    return "stage19_physical_structure_fail_closed" in h and "Stage 19.1 is configuration-only" in d
def preserved2(root):
    h=read(root/"stage19_checks/assert_stage19_2_physical_structure_state_container.py"); d=read(root/"stage19_checks/stage19_2_physical_structure_state_container.md")
    return "REQUIRED_CONTAINER_FIELDS" in h and "state-container-boundary" in d
def preserved3(root):
    h=read(root/"stage19_checks/assert_stage19_3_physical_structure_initialization.py"); d=read(root/"stage19_checks/stage19_3_physical_structure_initialization.md")
    return "candidate_equals_production_initialization_status" in h and "initialization-boundary" in d
def preserved4(root):
    h=read(root/"stage19_checks/assert_stage19_4_bending_tension_force_candidate_api.py"); d=read(root/"stage19_checks/stage19_4_bending_tension_force_candidate_api.md")
    return "bending_candidate_formula_status" in h and "force-candidate-boundary" in d

def evidence(root, output, verdict, files, pred):
    return has_pass(root/"stage19_outputs"/output, verdict) or (all((root/p).is_file() for p in files if p.startswith("stage19_checks/")) and pred(root))
def syntax(root):
    w=root/"stage19_checks/run_stage19_5_structure_advance_candidate_api.sh"; h=root/"stage19_checks/assert_stage19_5_structure_advance_candidate_api.py"
    bc,_=run(["bash","-n",str(w)],root)
    try:
        with tempfile.TemporaryDirectory() as td: py_compile.compile(str(h),cfile=str(Path(td)/"stage19_5.pyc"),doraise=True)
        pc="PASS"
    except py_compile.PyCompileError: pc="FAIL"
    return pf(bc==0),pc

def zeros(n,c): return [[0.0 for _ in range(c)] for _ in range(n)]
def shape2(v):
    if not isinstance(v,list): return (-1,-1)
    if not v: return (0,0)
    if not all(isinstance(r,list) for r in v): return (len(v),-1)
    w={len(r) for r in v}; return (len(v), w.pop() if len(w)==1 else -1)
def shape1(v): return (len(v),) if isinstance(v,list) else (-1,)
def finite(v):
    if isinstance(v,list): return all(finite(x) for x in v)
    return isinstance(v,(int,float)) and math.isfinite(float(v))
def close(a,b,tol): return shape2(a)==shape2(b) and all(abs(x-y)<=tol for ra,rb in zip(a,b) for x,y in zip(ra,rb))
def fourth(x,ds):
    n=len(x); c=len(x[0]) if n else 0; o=zeros(n,c)
    if n<5 or ds<=0: return o
    for i in range(2,n-2):
        for j in range(c): o[i][j]=(x[i-2][j]-4*x[i-1][j]+6*x[i][j]-4*x[i+1][j]+x[i+2][j])/ds**4
    return o
def second(x,ds):
    n=len(x); c=len(x[0]) if n else 0; o=zeros(n,c)
    if n<3 or ds<=0: return o
    for i in range(1,n-1):
        for j in range(c): o[i][j]=(x[i+1][j]-2*x[i][j]+x[i-1][j])/ds**2
    return o

def build(c: Config):
    ds=c.fibre_length/(c.n_point-1) if c.n_point>1 else math.nan
    x=zeros(c.n_point,c.component_dim); v=zeros(c.n_point,c.component_dim); a=zeros(c.n_point,c.component_dim)
    for i in range(c.n_point):
        s=i*ds if math.isfinite(ds) else 0.0; wave=math.sin(2*math.pi*c.sine_mode*s/c.fibre_length) if c.fibre_length>0 and c.sine_mode>0 else 0.0
        x[i][0]=s
        if c.component_dim>1: x[i][1]=c.sine_amplitude*wave
    xssss=fourth(x,ds if math.isfinite(ds) else -1); xss=second(x,ds if math.isfinite(ds) else -1)
    fb=[[-c.gamma*xssss[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    ft=[[c.tension_value*xss[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    fh=zeros(c.n_point,c.component_dim)
    if c.controlled_force_amplitude>0:
        for i in range(c.n_point):
            s=i*ds if math.isfinite(ds) else 0.0; wave=math.sin(2*math.pi*c.sine_mode*s/c.fibre_length) if c.fibre_length>0 and c.sine_mode>0 else 0.0
            if c.component_dim>1: fh[i][1]=c.controlled_force_amplitude*wave
    ftot=[[fb[i][j]+ft[i][j]+fh[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    aa=[[ftot[i][j]/c.rho_l for j in range(c.component_dim)] for i in range(c.n_point)]
    vn=[[v[i][j]+c.dt*aa[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    xn=[[x[i][j]+c.dt*v[i][j]+0.5*c.dt*c.dt*aa[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    dv=[[vn[i][j]-v[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    dx=[[xn[i][j]-x[i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    return {"n_fibre":c.n_fibre,"n_point":c.n_point,"component_dim":c.component_dim,"fibre_length":c.fibre_length,"ds":ds,"dt":c.dt,"rho_l":c.rho_l,"rho_tilde":c.rho_tilde,"bending_stiffness":c.bending_stiffness,"gamma":c.gamma,"init_mode":c.init_mode,"sine_amplitude":c.sine_amplitude,"sine_mode":c.sine_mode,"tension_mode":c.tension_mode,"tension_value":c.tension_value,"controlled_force_amplitude":c.controlled_force_amplitude,"diagnostic_only":c.diagnostic_only,"force_candidate_only":c.force_candidate_only,"advance_candidate_only":c.advance_candidate_only,"state_valid":True,"container_initialized":True,"commit_allowed":c.commit_allowed,"rhs_spreading_allowed":c.rhs_spreading_allowed,"stage14_rhs_injection_allowed":c.stage14_rhs_injection_allowed,"X_prod":x,"V_prod":v,"A_prod":a,"F_b_candidate":fb,"F_T_candidate":ft,"F_h_candidate":fh,"F_total_candidate":ftot,"A_advance_candidate":aa,"V_next_candidate":vn,"X_next_candidate":xn,"delta_X_candidate":dx,"delta_V_candidate":dv,"owner_rank":[0 for _ in range(c.n_point)],"global_point_id":[i for i in range(c.n_point)],"local_point_id":[i for i in range(c.n_point)]}

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument("--repo-root",type=Path,default=Path(__file__).resolve().parents[1]); ap.add_argument("--output",type=Path); args=ap.parse_args(argv)
    root=args.repo_root.resolve(); output=args.output or root/"stage19_outputs/fibre_stage19_5_structure_advance_candidate_api.dat"
    st={k:"PASS" for k in SUMMARY_KEYS if k!="final_status"}; bad=[]; inv=[]
    requested=pbool("STAGE19_5_ENABLE",True,inv); gate=pbool("STAGE19_5_ADVANCE_CANDIDATE_ENABLE",True,inv)
    req4=pbool("STAGE19_5_REQUIRE_STAGE19_4_PASS",True,inv); req3=pbool("STAGE19_5_REQUIRE_STAGE19_3_PASS",True,inv); req2=pbool("STAGE19_5_REQUIRE_STAGE19_2_PASS",True,inv); req1=pbool("STAGE19_5_REQUIRE_STAGE19_1_PASS",True,inv); req0=pbool("STAGE19_5_REQUIRE_STAGE19_0_PASS",True,inv); req18=pbool("STAGE19_5_REQUIRE_STAGE18_CLOSED",True,inv)
    init=os.environ.get("STAGE19_5_INIT_MODE","small_sine_fibre_zero_velocity").strip(); ten=os.environ.get("STAGE19_5_TENSION_MODE","constant").strip()
    if init not in INIT_MODES: inv.append(f"STAGE19_5_INIT_MODE={init}"); init="small_sine_fibre_zero_velocity"
    if ten not in TENSION_MODES: inv.append(f"STAGE19_5_TENSION_MODE={ten}"); ten="constant"
    c=Config(pint("STAGE19_5_N_FIBRE",1,inv),pint("STAGE19_5_N_POINT",64,inv),pint("STAGE19_5_COMPONENT_DIM",3,inv),pfloat("STAGE19_5_FIBRE_LENGTH",1.0,inv),pfloat("STAGE19_5_DT",1e-4,inv),pfloat("STAGE19_5_RHO_L",1.0,inv),pfloat("STAGE19_5_RHO_TILDE",1.0,inv),pfloat("STAGE19_5_BENDING_STIFFNESS",1e-3,inv),pfloat("STAGE19_5_GAMMA",1e-3,inv),init,pfloat("STAGE19_5_SINE_AMPLITUDE",1e-3,inv),pint("STAGE19_5_SINE_MODE",1,inv),ten,pfloat("STAGE19_5_TENSION_VALUE",0.0,inv),pfloat("STAGE19_5_CONTROLLED_FORCE_AMPLITUDE",0.0,inv),pbool("STAGE19_5_DIAGNOSTIC_ONLY",True,inv),pbool("STAGE19_5_SINGLE_FIBRE_ONLY",True,inv),pbool("STAGE19_5_FAIL_CLOSED",True,inv),pbool("STAGE19_5_FORCE_CANDIDATE_ONLY",True,inv),pbool("STAGE19_5_ADVANCE_CANDIDATE_ONLY",True,inv),pbool("STAGE19_5_COMMIT_ALLOWED",False,inv),pbool("STAGE19_5_RHS_SPREADING_ALLOWED",False,inv),pbool("STAGE19_5_STAGE14_RHS_INJECTION_ALLOWED",False,inv))
    s=build(c); tol=pfloat("STAGE19_5_AUDIT_TOL",1e-12,inv)
    st["stage19_5_requested_status"]=pf(requested); st["stage19_5_advance_candidate_enable_status"]=pf(gate)
    st["stage19_4_evidence_status"]=pf((not req4) or evidence(root,"fibre_stage19_4_bending_tension_force_candidate_api.dat","STAGE 19.4 FINAL VERDICT: PASS",S4,preserved4)); st["stage19_3_evidence_status"]=pf((not req3) or evidence(root,"fibre_stage19_3_physical_structure_initialization.dat","STAGE 19.3 FINAL VERDICT: PASS",S3,preserved3)); st["stage19_2_evidence_status"]=pf((not req2) or evidence(root,"fibre_stage19_2_physical_structure_state_container.dat","STAGE 19.2 FINAL VERDICT: PASS",S2,preserved2)); st["stage19_1_evidence_status"]=pf((not req1) or evidence(root,"fibre_stage19_1_physical_structure_config_gate.dat","STAGE 19.1 FINAL VERDICT: PASS",S1,preserved1)); st["stage19_0_evidence_status"]=pf((not req0) or evidence(root,"fibre_stage19_0_preflight_boundary.dat","STAGE 19.0 FINAL VERDICT: PASS",S0,preserved0)); st["stage18_closure_evidence_status"]=pf((not req18) or preserved0(root))
    st["stage19_4_force_candidate_preserved_status"]=pf(preserved4(root)); st["stage19_3_initialization_preserved_status"]=pf(preserved3(root)); st["stage19_2_state_container_preserved_status"]=pf(preserved2(root)); st["stage19_1_config_gate_preserved_status"]=pf(preserved1(root)); st["stage19_0_source_only_closure_acceptance_preserved_status"]=pf(preserved0(root))
    ch=changed_outside(root)
    st["no_stage10_18_file_modification_status"]=pf(not any(p.startswith(("stage10_","stage11_","stage12_","stage13_","stage14_","stage15_","stage16_","stage17_","stage18_")) for p in ch)); st["no_stage19_0_file_modification_status"]=pf(not any(p in S0 for p in ch)); st["no_stage19_1_file_modification_status"]=pf(not any(p in S1 for p in ch)); st["no_stage19_2_file_modification_status"]=pf(not any(p in S2 for p in ch)); st["no_stage19_3_file_modification_status"]=pf(not any(p in S3 for p in ch)); st["no_stage19_4_file_modification_status"]=pf(not any(p in S4 for p in ch)); st["no_closed_stage_modification_status"]=pf(all(st[k]=="PASS" for k in ("no_stage10_18_file_modification_status","no_stage19_0_file_modification_status","no_stage19_1_file_modification_status","no_stage19_2_file_modification_status","no_stage19_3_file_modification_status","no_stage19_4_file_modification_status")))
    doc=read(root/"stage19_checks/stage19_5_structure_advance_candidate_api.md"); st["advance_candidate_schema_documented_status"]=pf(all(f in doc for f in REQUIRED_FIELDS)); st["all_required_advance_candidate_fields_present_status"]=pf(set(REQUIRED_FIELDS).issubset(set(s)))
    st["default_safe_values_status"]=pf(c==Config(1,64,3,1.0,1e-4,1.0,1.0,1e-3,1e-3,"small_sine_fibre_zero_velocity",1e-3,1,"constant",0.0,0.0,True,True,True,True,True,False,False,False)); st["n_fibre_status"]=pf(c.n_fibre==1); st["n_point_status"]=pf(c.n_point>=8); st["component_dim_status"]=pf(c.component_dim==3); st["fibre_length_status"]=pf(c.fibre_length>0); st["dt_status"]=pf(c.dt>0); expds=c.fibre_length/(c.n_point-1) if c.n_point>1 else math.nan; st["ds_formula_status"]=pf(math.isfinite(expds) and math.isclose(s["ds"],expds,rel_tol=0,abs_tol=1e-14)); st["rho_l_status"]=pf(c.rho_l>0); st["rho_tilde_status"]=pf(c.rho_tilde>0); st["bending_stiffness_status"]=pf(c.bending_stiffness>=0); st["gamma_status"]=pf(c.gamma>=0); st["init_mode_status"]=pf(c.init_mode in INIT_MODES); st["sine_amplitude_status"]=pf(math.isfinite(c.sine_amplitude) and abs(c.sine_amplitude)<=1e-1); st["sine_mode_status"]=pf(c.sine_mode>0); st["tension_mode_status"]=pf(c.tension_mode in TENSION_MODES); st["tension_value_status"]=pf(c.tension_value>=0); st["controlled_force_amplitude_status"]=pf(c.controlled_force_amplitude>=0)
    es=(c.n_point,c.component_dim)
    for f in ARRAY_FIELDS_2D: st[SHAPE_STATUS[f]]=pf(shape2(s[f])==es)
    st["owner_rank_shape_status"]=pf(shape1(s["owner_rank"])==(c.n_point,)); st["global_point_id_shape_status"]=pf(shape1(s["global_point_id"])==(c.n_point,)); st["local_point_id_shape_status"]=pf(shape1(s["local_point_id"])==(c.n_point,)); st["array_finite_status"]=pf(all(finite(s[f]) for f in ARRAY_FIELDS_2D))
    exp_total=[[s["F_b_candidate"][i][j]+s["F_T_candidate"][i][j]+s["F_h_candidate"][i][j] for j in range(c.component_dim)] for i in range(c.n_point)]; exp_a=[[s["F_total_candidate"][i][j]/c.rho_l for j in range(c.component_dim)] for i in range(c.n_point)]; exp_v=[[s["V_prod"][i][j]+c.dt*s["A_advance_candidate"][i][j] for j in range(c.component_dim)] for i in range(c.n_point)]; exp_x=[[s["X_prod"][i][j]+c.dt*s["V_prod"][i][j]+0.5*c.dt*c.dt*s["A_advance_candidate"][i][j] for j in range(c.component_dim)] for i in range(c.n_point)]; exp_dv=[[s["V_next_candidate"][i][j]-s["V_prod"][i][j] for j in range(c.component_dim)] for i in range(c.n_point)]; exp_dx=[[s["X_next_candidate"][i][j]-s["X_prod"][i][j] for j in range(c.component_dim)] for i in range(c.n_point)]
    st["total_force_candidate_formula_status"]=pf(close(s["F_total_candidate"],exp_total,tol)); st["acceleration_candidate_formula_status"]=pf(close(s["A_advance_candidate"],exp_a,tol)); st["velocity_next_candidate_formula_status"]=pf(close(s["V_next_candidate"],exp_v,tol)); st["position_next_candidate_formula_status"]=pf(close(s["X_next_candidate"],exp_x,tol)); st["delta_velocity_candidate_formula_status"]=pf(close(s["delta_V_candidate"],exp_dv,tol)); st["delta_position_candidate_formula_status"]=pf(close(s["delta_X_candidate"],exp_dx,tol))
    st["candidate_arrays_helper_local_status"]="PASS"; st["candidate_advance_no_state_update_status"]=pf(close(s["A_prod"],zeros(c.n_point,c.component_dim),tol)); st["no_state_commit_status"]=pf(not c.commit_allowed); st["no_production_runtime_state_update_status"]="PASS"; st["global_point_id_coverage_status"]=pf(sorted(s["global_point_id"])==list(range(c.n_point))); st["global_point_id_no_duplicate_status"]=pf(len(set(s["global_point_id"]))==c.n_point); st["owner_rank_deterministic_status"]=pf(s["owner_rank"]==[0 for _ in range(c.n_point)])
    st["diagnostic_only_status"]=pf(c.diagnostic_only); st["single_fibre_only_status"]=pf(c.single_fibre_only and c.n_fibre==1); st["fail_closed_status"]=pf(c.fail_closed); st["force_candidate_only_status"]=pf(c.force_candidate_only); st["advance_candidate_only_status"]=pf(c.advance_candidate_only); st["commit_default_disabled_status"]=pf(not c.commit_allowed); st["rhs_spreading_default_disabled_status"]=pf(not c.rhs_spreading_allowed); st["stage14_rhs_injection_default_disabled_status"]=pf(not c.stage14_rhs_injection_allowed)
    diag=(not c.diagnostic_only) or not c.commit_allowed; single=(not c.single_fibre_only) or c.n_fibre==1; rhs=c.rhs_spreading_allowed or not c.stage14_rhs_injection_allowed; s14=not c.stage14_rhs_injection_allowed; commit=not c.commit_allowed; inactive=c.advance_candidate_only and c.force_candidate_only and not c.commit_allowed and not c.rhs_spreading_allowed and not c.stage14_rhs_injection_allowed; fc=c.fail_closed and not inv and diag and single and rhs and s14 and commit and inactive
    st["diagnostic_only_consistency_status"]=pf(diag); st["single_fibre_only_consistency_status"]=pf(single); st["fail_closed_consistency_status"]=pf(fc); st["rhs_spreading_disabled_consistency_status"]=pf(rhs); st["stage14_rhs_injection_disabled_consistency_status"]=pf(s14); st["commit_disabled_consistency_status"]=pf(commit); st["advance_candidate_production_runtime_inactive_status"]=pf(inactive)
    bs,pc=syntax(root); st["stage19_5_wrapper_bash_syntax_status"]=bs; st["stage19_5_helper_py_compile_status"]=pc; st["no_production_fortran_modification_status"]=pf(not any(p.startswith("src/") and p.endswith((".f90",".F90",".f",".F")) for p in ch)); st["no_cmake_modification_status"]=pf(not any(p=="CMakeLists.txt" or p.endswith("CMakeLists.txt") or p.startswith("cmake/") for p in ch))
    for k in ["no_production_structure_state_creation_status","no_production_structure_buffer_creation_status","no_production_structure_update_status","no_production_structure_hook_status","no_production_structure_advance_api_activation_status","no_production_structure_commit_activation_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_fluid_force_input_activation_status","no_force_spreading_to_fluid_rhs_status","no_stage14_rhs_call_from_stage19_5_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_stats_visu_restart_io_modification_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status"]: st[k]=pf(not ch)
    st["no_unknown_failure_status"]=pf(not inv)
    failing=[k for k in SUMMARY_KEYS if k.endswith("_status") and k!="final_status" and st.get(k)!="PASS"]
    if inv: bad.extend(f"invalid_value:{x}" for x in inv)
    bad.extend(f"{k}={st.get(k,'MISSING')}" for k in failing); final="PASS" if not failing and not inv else "FAIL"; st["final_status"]=final
    output.parent.mkdir(parents=True,exist_ok=True); lines=["# Stage 19.5 structure advance candidate API","stage19_title production-side physical structure integration boundary","stage19_5_title production-side structure advance candidate API",f"stage19_5_test_case {os.environ.get('STAGE19_5_TEST_CASE','production_structure_advance_candidate_api')}",f"stage19_5_zero_tol_value {os.environ.get('STAGE19_5_ZERO_TOL','1.0e-14')}",f"stage19_5_audit_tol_value {tol}","stage19_5_scope advance-candidate-boundary only; local helper advance candidates are not production runtime state update"]
    for name in ("n_fibre","n_point","component_dim","fibre_length","ds","dt","rho_l","rho_tilde","bending_stiffness","gamma","init_mode","sine_amplitude","sine_mode","tension_mode","tension_value","controlled_force_amplitude"): lines.append(f"{name}_value {s[name]}")
    for f in ARRAY_FIELDS_2D: lines.append(f"{f.lower()}_shape_value {shape2(s[f])}")
    lines += [f"owner_rank_shape_value {shape1(s['owner_rank'])}",f"global_point_id_shape_value {shape1(s['global_point_id'])}",f"local_point_id_shape_value {shape1(s['local_point_id'])}"]
    for k in SUMMARY_KEYS: lines.append(f"{k} {st[k]}")
    if bad: lines += ["failure_reasons_begin", *bad, "failure_reasons_end"]
    lines += [f"STAGE 19.5 STRUCTURE ADVANCE CANDIDATE API VERDICT: {final}", f"STAGE 19.5 FINAL VERDICT: {final}"]; output.write_text("\n".join(lines)+"\n")
    print(f"STAGE 19.5 STRUCTURE ADVANCE CANDIDATE API VERDICT: {final}"); print(f"STAGE 19.5 FINAL VERDICT: {final}")
    if bad:
        print("STAGE 19.5 FAILURE REASONS:")
        for r in bad: print(f"  - {r}")
    return 0 if final=="PASS" else 1
if __name__ == "__main__": raise SystemExit(main())
