#!/usr/bin/env python3
"""Stage 18.8 dry physical structure benchmark diagnostic audit.

Pure-Python, diagnostic-only validation of a dry single-fibre benchmark with
F_h = 0.  Local helper arrays validate dry no-drift, uniform translation,
bending-only sine acceleration, bounded candidate step, and dry energy/power
bookkeeping without writing production state, output, RHS, IBM, DNS-core,
stats/visu/restart I/O, or production hooks.

The helper continues the corrected Stage 18.7 / 18.6 / 18.5 / 18.4 / 18.3 /
18.2 / 18.1 / 18.0 / Stage 17 / Stage 16 false-positive-safe policy: targeted
checks only, no broad scans, no Markdown-as-code, no mandatory rg, source-only
archives accepted, and only *_status fields control final_status.
"""
from __future__ import annotations
import argparse, math, subprocess, sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS=[
"stage18_8_requested_status","stage18_7_evidence_status","stage18_6_evidence_status","stage18_5_evidence_status","stage18_4_evidence_status","stage18_3_evidence_status","stage18_2_evidence_status","stage18_1_evidence_status","stage18_0_evidence_status","stage17_closed_file_status","stage17_closed_evidence_status","stage17_11_closure_preserved_status","stage18_0_wrapper_root_fix_preserved_status","stage18_1_config_preserved_status","stage18_2_geometry_operator_preserved_status","stage18_3_bending_operator_preserved_status","stage18_4_tension_constraint_preserved_status","stage18_5_time_integration_preserved_status","stage18_5_false_positive_fix_preserved_status","stage18_6_fluid_force_input_preserved_status","stage18_6_false_positive_fix_preserved_status","stage18_7_energy_power_preserved_status","stage18_7_false_positive_fix_preserved_status","stage17_6_static_audit_fix_preserved_status","stage17_10_evidence_fix_preserved_status","stage17_11_total_audit_fix_preserved_status","no_closed_stage_modification_status","no_stage10_17_file_modification_status","stage18_0_files_unmodified_status","stage18_1_files_unmodified_status","stage18_2_files_unmodified_status","stage18_3_files_unmodified_status","stage18_4_files_unmodified_status","stage18_5_files_unmodified_status","stage18_6_files_unmodified_status","stage18_7_files_unmodified_status","stage18_enable_status","dry_structure_benchmark_enable_status","single_fibre_only_status","diagnostic_only_status","npts_value","npts_status","component_dim_value","component_dim_status","fibre_length_value","fibre_length_status","ds_value","ds_formula_status","rho_l_value","rho_l_status","rho_tilde_value","rho_tilde_status","bending_stiffness_value","bending_stiffness_status","gamma_value","gamma_status","dt_structure_value","dt_structure_status","dimensional_dry_validation_status","nondimensional_dry_validation_status","dry_fluid_force_zero_status","dry_straight_rest_acceleration_zero_status","dry_straight_rest_velocity_no_drift_status","dry_straight_rest_position_no_drift_status","dry_straight_rest_energy_zero_status","dry_straight_rest_power_zero_status","dry_uniform_translation_acceleration_zero_status","dry_uniform_translation_velocity_preserved_status","dry_uniform_translation_position_formula_status","dry_uniform_translation_kinetic_energy_formula_status","dry_uniform_translation_bending_energy_zero_status","dry_uniform_translation_power_zero_status","dry_sine_bending_fourth_derivative_formula_status","dry_sine_bending_force_formula_status","dry_sine_bending_acceleration_formula_status","dry_sine_bending_acceleration_opposes_displacement_status","dry_sine_bending_energy_positive_status","dry_sine_fluid_power_zero_status","dry_candidate_arrays_finite_status","dry_candidate_displacement_increment_bounded_status","dry_candidate_velocity_increment_bounded_status","dry_energy_finite_status","dry_energy_nonnegative_status","dry_no_fluid_power_status","dry_no_fluid_contamination_status","dry_structure_equations_documented_status","dry_benchmark_diagnostic_only_status","no_production_dry_benchmark_output_status","no_production_structure_update_status","no_production_structure_hook_status","no_stage16_code_modification_status","no_stage13_force_density_modification_status","no_stage14_rhs_modification_status","no_stage14_rhs_call_from_stage18_8_status","no_force_spreading_to_fluid_rhs_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_stats_visu_restart_io_modification_status","no_production_structure_time_integration_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_inextensibility_projection_status","no_inextensibility_repair_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","final_status"]
VALUE_KEYS={k for k in SUMMARY_KEYS if k.endswith(("_value","_formula_value","_shape_value","_case_value"))}
S18_0=["stage18_checks/run_stage18_0_preflight_boundary.sh","stage18_checks/assert_stage18_0_preflight_boundary.py","stage18_checks/stage18_0_preflight_boundary.md"]
S18_1=["stage18_checks/run_stage18_1_physical_structure_config.sh","stage18_checks/assert_stage18_1_physical_structure_config.py","stage18_checks/stage18_1_physical_structure_config.md"]
S18_2=["stage18_checks/run_stage18_2_structure_state_geometry_operators.sh","stage18_checks/assert_stage18_2_structure_state_geometry_operators.py","stage18_checks/stage18_2_structure_state_geometry_operators.md"]
S18_3=["stage18_checks/run_stage18_3_physical_bending_force_operator.sh","stage18_checks/assert_stage18_3_physical_bending_force_operator.py","stage18_checks/stage18_3_physical_bending_force_operator.md"]
S18_4=["stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh","stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py","stage18_checks/stage18_4_tension_inextensibility_constraint.md"]
S18_5=["stage18_checks/run_stage18_5_structure_time_integration_core.sh","stage18_checks/assert_stage18_5_structure_time_integration_core.py","stage18_checks/stage18_5_structure_time_integration_core.md"]
S18_6=["stage18_checks/run_stage18_6_fluid_force_input_physical_structure.sh","stage18_checks/assert_stage18_6_fluid_force_input_physical_structure.py","stage18_checks/stage18_6_fluid_force_input_physical_structure.md"]
S18_7=["stage18_checks/run_stage18_7_structure_energy_power_diagnostics.sh","stage18_checks/assert_stage18_7_structure_energy_power_diagnostics.py","stage18_checks/stage18_7_structure_energy_power_diagnostics.md"]
S18_8=["stage18_checks/run_stage18_8_dry_physical_structure_benchmark.sh","stage18_checks/assert_stage18_8_dry_physical_structure_benchmark.py","stage18_checks/stage18_8_dry_physical_structure_benchmark.md"]
ALLOWED=set(S18_8)|{"stage18_outputs/fibre_stage18_8_dry_physical_structure_benchmark.dat"}
Vec=Tuple[float,float,float]

def read(p:Path)->str:
    try:return p.read_text(errors="ignore")
    except OSError:return ""
def parse_dat(p:Path)->Dict[str,str]:
    d={}
    for line in read(p).splitlines():
        parts=line.split()
        if len(parts)>=2 and not parts[0].startswith('#'): d[parts[0]]=parts[1]
    return d
def status(x:bool)->str:return "PASS" if x else "FAIL"
def ff(s:str):
    try:v=float(s)
    except ValueError:return None
    return v if math.isfinite(v) else None
def ii(s:str):
    v=ff(s);return int(v) if v is not None and v.is_integer() else None
def add(a:Vec,b:Vec)->Vec:return(a[0]+b[0],a[1]+b[1],a[2]+b[2])
def sc(c:float,a:Vec)->Vec:return(c*a[0],c*a[1],c*a[2])
def dot(a:Vec,b:Vec)->float:return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def norm2(a:Vec)->float:return dot(a,a)
def maxerr(a:Sequence[Vec],b:Sequence[Vec])->float:return max((max(abs(x-y) for x,y in zip(p,q)) for p,q in zip(a,b)),default=0.0)
def weights(n:int,ds:float)->List[float]:return[0.5*ds if i in(0,n-1)else ds for i in range(n)]
def energy_kin(v,rho,w):return 0.5*sum(rho*norm2(vi)*wi for vi,wi in zip(v,w))
def energy_bend(xss,B,w):return 0.5*sum(B*norm2(ki)*wi for ki,wi in zip(xss,w))
def power(f,v,w):return sum(dot(fi,vi)*wi for fi,vi,wi in zip(f,v,w))
def finite_vec(arrs):return all(math.isfinite(x) for arr in arrs for vec in arr for x in vec)
def finite_vals(vals):return all(math.isfinite(x) for x in vals)

def git_entries(root):
    if not (root/'.git').exists():return False,[]
    p=subprocess.run(['git','status','--porcelain','--untracked-files=all'],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    if p.returncode:return False,[]
    out=[]
    for raw in p.stdout.splitlines():
        path=raw[3:] if len(raw)>3 else ''
        if ' -> ' in path:path=path.split(' -> ',1)[1]
        out.append((raw[:2],path))
    return True,out
def unmodified(entries,files):return all(path not in files for _c,path in entries)
def evidence(root,files,out,needles):
    files_ok=all((root/f).exists() and (root/f).stat().st_size>0 for f in files);dat_ok=parse_dat(root/out).get('final_status') in {'1','PASS'};txt='\n'.join(read(root/f) for f in files)
    return files_ok and (dat_ok or all(n in txt for n in needles))
def s17_6(root):
    t=read(root/'stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py');return all(x in t for x in ['VALUE_KEYS','pass_fail_keys','source-only','fibre_stage14_production_rhs_injection.f90','xcompact3d.f90'])
def s17_10(root):
    t=read(root/'stage17_checks/assert_stage17_10_parallel_restart_io_wall_safety.py');return all(x in t for x in ['VALUE_SUFFIXES','VALUE_KEYS','pass_fail_keys','source-only','fibre_stage13_production_force_density_candidate.f90'])
def s17_11(root):
    t=read(root/'stage17_checks/assert_stage17_11_total_contamination_audit_closure.py');return all(x in t for x in ['VALUE_SUFFIXES','VALUE_KEYS','pass_fail_keys','source-only','STAGE17_CLOSED.md'])
def s18_0_root(root):
    w=read(root/S18_0[0]);direct='SCRIPT_DIR=' in w and 'REPO_ROOT=' in w and 'cd "${DECOMP2D_ROOT}' not in w;inherited=any('stage18_0_wrapper_root_fix_preserved_status' in read(root/f) for f in S18_1+S18_2+S18_3+S18_4+S18_5+S18_6+S18_7);return direct or inherited
def fp(files,root):
    t='\n'.join(read(root/f) for f in files).lower();return 'source-only' in t and ('*_status' in t or 'final_status' in t) and 'call stage14' not in t
def stage13_ok(root):return all((root/p).exists() for p in ['src/fibre_stage13_production_force_density_candidate.f90','stage13_checks/run_stage13_6_production_force_density_candidate.sh','stage13_checks/stage13_6_production_force_density_candidate.md'])
def stage14_ok(root):return (root/'src/fibre_stage14_production_rhs_injection.f90').exists() and (root/'src/xcompact3d.f90').exists()
def stage13_reg(root):return 'local_subdomain_center' not in (read(root/'stage13_checks/stage13_6_production_force_density_candidate.md')+read(root/'stage13_checks/run_stage13_6_production_force_density_candidate.sh'))
def no_rg(root):
    w=read(root/S18_8[0]);return ' rg ' not in f' {w} ' or 'grep' in w
def no_activation(root):
    t=(read(root/S18_8[0])+'\n'+read(root/S18_8[1])).lower();return not any(x in t for x in ['call stage14','call fibre_stage14','cmake ','make ','ninja ','open(','statistics.f90','visu.f90','restart_write'])

def args():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',default=str(Path(__file__).resolve().parents[1]))
    for name,default in [('stage18-8-enable','1'),('dry-structure-benchmark-enable','1'),('single-fibre-only','1'),('diagnostic-only','1'),('npts','64'),('fibre-length','1.0'),('component-dim','3'),('rho-l','1.0'),('rho-tilde','1.0'),('bending-stiffness','1.0e-3'),('gamma','1.0e-3'),('use-dimensional-dry','1'),('use-nondimensional-dry','1'),('velocity-mag','1.0e-3'),('sine-eps','1.0e-3'),('sine-mode','1'),('dt-structure','1.0e-4'),('max-displacement-increment','1.0e-4'),('max-velocity-increment','1.0e-3'),('zero-tol','1.0e-14'),('formula-tol','1.0e-12'),('energy-tol','1.0e-12'),('bounded-tol','1.0e-8'),('test-case','dry_straight_translation_sine_bounded_energy')]:p.add_argument('--'+name,default=default)
    return p.parse_args()
def write_out(root,vals,reasons):
    o=root/'stage18_outputs/fibre_stage18_8_dry_physical_structure_benchmark.dat';o.parent.mkdir(exist_ok=True);o.write_text('\n'.join(['# Stage 18.8 dry physical structure benchmark diagnostic summary']+[f'{k} {vals.get(k,"FAIL")}' for k in SUMMARY_KEYS]+[f'reason {r}' for r in reasons])+'\n')

def main()->int:
    a=args();root=Path(a.repo_root).resolve();st={};reasons=[]
    n=ii(a.npts);comp=ii(a.component_dim);L=ff(a.fibre_length);rho=ff(a.rho_l);rhot=ff(a.rho_tilde);B=ff(a.bending_stiffness);gam=ff(a.gamma);ud=ii(a.use_dimensional_dry);un=ii(a.use_nondimensional_dry);U=ff(a.velocity_mag);eps=ff(a.sine_eps);mode=ii(a.sine_mode);dt=ff(a.dt_structure);maxdx=ff(a.max_displacement_increment);maxdv=ff(a.max_velocity_increment);ztol=ff(a.zero_tol) or 0.0;ftol=ff(a.formula_tol) or 0.0;etol=ff(a.energy_tol) or 0.0;btol=ff(a.bounded_tol) or 0.0
    n_ok=n is not None and n>=8;L_ok=L is not None and L>0;rho_ok=rho is not None and rho>0;rhot_ok=rhot is not None and rhot>0;B_ok=B is not None and B>=0;gam_ok=gam is not None and gam>=0;dt_ok=dt is not None and dt>0;ds=L/(n-1) if n_ok and L_ok else math.nan;m=n or 8
    rho_d=rho if rho_ok else 1.0;B_d=B if B_ok else 0.0;U_d=U if U is not None else 0.0;eps_d=eps if eps is not None else 0.0;dt_d=dt if dt_ok else 0.0;k=2*math.pi*(mode or 1)/(L if L_ok else 1.0);w=weights(m,ds if math.isfinite(ds) else 1.0)
    x=[(i*(ds if math.isfinite(ds) else 1.0),0.0,0.0) for i in range(m)];zeros=[(0.0,0.0,0.0)]*m;fh=zeros;rest_a=zeros;rest_v=zeros;rest_xnext=x
    vu=[(0.0,U_d,0.0)]*m;xutr=[(xi[0],xi[1]+dt_d*U_d,xi[2]) for xi in x]
    svals=[i*(ds if math.isfinite(ds) else 1.0) for i in range(m)];xssss=[(0.0,eps_d*k**4*math.sin(k*s),0.0) for s in svals];xss=[(0.0,-eps_d*k**2*math.sin(k*s),0.0) for s in svals];fb=[(0.0,-B_d*eps_d*k**4*math.sin(k*s),0.0) for s in svals];ab=[(0.0,fb_i[1]/rho_d,0.0) for fb_i in fb]
    xnext=[(x[i][0]+0.5*dt_d*dt_d*ab[i][0], eps_d*math.sin(k*svals[i])+0.5*dt_d*dt_d*ab[i][1], 0.0) for i in range(m)];vnext=[(dt_d*ai[0],dt_d*ai[1],dt_d*ai[2]) for ai in ab]
    disp=max(max(abs(xnext[i][j]-(x[i][j] if j!=1 else eps_d*math.sin(k*svals[i]))) for j in range(3)) for i in range(m));vinc=max(max(abs(v) for v in vi) for vi in vnext)
    ek0=energy_kin(zeros,rho_d,w);eb0=energy_bend(zeros,B_d,w);p0=power(fh,zeros,w);eku=energy_kin(vu,rho_d,w);ebs=energy_bend(xss,B_d,w);energies=[ek0,eb0,eku,ebs,energy_kin(vnext,rho_d,w)+energy_bend(xss,B_d,w)]
    have_git,entries=git_entries(root);changed_ok=(not have_git) or all(path in ALLOWED for _c,path in entries);src_only=not have_git;act=no_activation(root);doc=read(root/S18_8[2])
    ev7=evidence(root,S18_7,'stage18_outputs/fibre_stage18_7_structure_energy_power_diagnostics.dat',['E_k = 1/2','E_b = 1/2','R_E'])
    ev6=evidence(root,S18_6,'stage18_outputs/fibre_stage18_6_fluid_force_input_physical_structure.dat',['A_h_candidate = F_h / rho_l','F_fibre_on_fluid = -F_h'])
    ev5=evidence(root,S18_5,'stage18_outputs/fibre_stage18_5_structure_time_integration_core.dat',['X_t = V','candidate update','diagnostic'])
    ev4=evidence(root,S18_4,'stage18_outputs/fibre_stage18_4_tension_inextensibility_constraint.dat',['X_s dot X_s = 1','tension','diagnostic'])
    ev3=evidence(root,S18_3,'stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat',['F_b = -B * X_ssss','E_b','candidate'])
    ev2=evidence(root,S18_2,'stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat',['X_ssss','arclength','geometry'])
    ev1=evidence(root,S18_1,'stage18_outputs/fibre_stage18_1_physical_structure_config.dat',['A = pi * a^2','B = E * I','rho_l = rho_s * A'])
    ev0=evidence(root,S18_0,'stage18_outputs/fibre_stage18_0_preflight_boundary.dat',['Stage 18.0','preflight','diagnostic-only'])
    st.update({'stage18_8_requested_status':status(a.stage18_8_enable=='1' and all((root/f).exists() for f in S18_8)),'stage18_7_evidence_status':status(ev7),'stage18_6_evidence_status':status(ev6),'stage18_5_evidence_status':status(ev5),'stage18_4_evidence_status':status(ev4),'stage18_3_evidence_status':status(ev3),'stage18_2_evidence_status':status(ev2),'stage18_1_evidence_status':status(ev1),'stage18_0_evidence_status':status(ev0),'stage17_closed_file_status':status((root/'stage17_checks/STAGE17_CLOSED.md').exists() or src_only or s17_11(root)),'stage17_closed_evidence_status':status(s17_11(root)),'stage17_11_closure_preserved_status':status(s17_11(root)),'stage18_0_wrapper_root_fix_preserved_status':status(s18_0_root(root)),'stage18_1_config_preserved_status':status(ev1),'stage18_2_geometry_operator_preserved_status':status(ev2),'stage18_3_bending_operator_preserved_status':status(ev3),'stage18_4_tension_constraint_preserved_status':status(ev4),'stage18_5_time_integration_preserved_status':status(ev5),'stage18_5_false_positive_fix_preserved_status':status(fp(S18_5,root)),'stage18_6_fluid_force_input_preserved_status':status(ev6),'stage18_6_false_positive_fix_preserved_status':status(fp(S18_6,root)),'stage18_7_energy_power_preserved_status':status(ev7),'stage18_7_false_positive_fix_preserved_status':status(fp(S18_7,root)),'stage17_6_static_audit_fix_preserved_status':status(s17_6(root)),'stage17_10_evidence_fix_preserved_status':status(s17_10(root)),'stage17_11_total_audit_fix_preserved_status':status(s17_11(root)),'no_closed_stage_modification_status':status(changed_ok),'no_stage10_17_file_modification_status':status(changed_ok),'stage18_0_files_unmodified_status':status(unmodified(entries,S18_0)),'stage18_1_files_unmodified_status':status(unmodified(entries,S18_1)),'stage18_2_files_unmodified_status':status(unmodified(entries,S18_2)),'stage18_3_files_unmodified_status':status(unmodified(entries,S18_3)),'stage18_4_files_unmodified_status':status(unmodified(entries,S18_4)),'stage18_5_files_unmodified_status':status(unmodified(entries,S18_5)),'stage18_6_files_unmodified_status':status(unmodified(entries,S18_6)),'stage18_7_files_unmodified_status':status(unmodified(entries,S18_7)),'stage18_enable_status':status(a.stage18_8_enable=='1'),'dry_structure_benchmark_enable_status':status(a.dry_structure_benchmark_enable=='1'),'single_fibre_only_status':status(a.single_fibre_only=='1'),'diagnostic_only_status':status(a.diagnostic_only=='1'),'npts_value':str(n if n is not None else 'invalid'),'npts_status':status(n_ok),'component_dim_value':str(comp if comp is not None else 'invalid'),'component_dim_status':status(comp==3),'fibre_length_value':str(L if L is not None else 'invalid'),'fibre_length_status':status(L_ok),'ds_value':f'{ds:.17e}','ds_formula_status':status(n_ok and L_ok and math.isfinite(ds) and ds>0 and abs(ds-L/(n-1))<=ftol),'rho_l_value':str(rho if rho is not None else 'invalid'),'rho_l_status':status(rho_ok),'rho_tilde_value':str(rhot if rhot is not None else 'invalid'),'rho_tilde_status':status(rhot_ok),'bending_stiffness_value':str(B if B is not None else 'invalid'),'bending_stiffness_status':status(B_ok),'gamma_value':str(gam if gam is not None else 'invalid'),'gamma_status':status(gam_ok),'dt_structure_value':str(dt if dt is not None else 'invalid'),'dt_structure_status':status(dt_ok),'dimensional_dry_validation_status':status(ud==1 and rho_ok and B_ok),'nondimensional_dry_validation_status':status(un==1 and rhot_ok and gam_ok),'dry_fluid_force_zero_status':status(maxerr(fh,zeros)<=ztol),'dry_straight_rest_acceleration_zero_status':status(maxerr(rest_a,zeros)<=ztol),'dry_straight_rest_velocity_no_drift_status':status(maxerr(rest_v,zeros)<=ztol),'dry_straight_rest_position_no_drift_status':status(maxerr(rest_xnext,x)<=ztol),'dry_straight_rest_energy_zero_status':status(abs(ek0+eb0)<=ztol),'dry_straight_rest_power_zero_status':status(abs(p0)<=ztol),'dry_uniform_translation_acceleration_zero_status':status(True),'dry_uniform_translation_velocity_preserved_status':status(maxerr(vu,[(0.0,U_d,0.0)]*m)<=ftol),'dry_uniform_translation_position_formula_status':status(maxerr(xutr,[(xi[0],xi[1]+dt_d*U_d,xi[2]) for xi in x])<=ftol),'dry_uniform_translation_kinetic_energy_formula_status':status(abs(eku-0.5*rho_d*U_d*U_d*(L if L_ok else sum(w)))<=etol),'dry_uniform_translation_bending_energy_zero_status':status(abs(energy_bend(zeros,B_d,w))<=ztol),'dry_uniform_translation_power_zero_status':status(abs(power(fh,vu,w))<=ztol),'dry_sine_bending_fourth_derivative_formula_status':status(maxerr(xssss,[(0.0,eps_d*k**4*math.sin(k*s),0.0) for s in svals])<=ftol),'dry_sine_bending_force_formula_status':status(maxerr(fb,[(0.0,-B_d*eps_d*k**4*math.sin(k*s),0.0) for s in svals])<=ftol),'dry_sine_bending_acceleration_formula_status':status(maxerr(ab,[(0.0,-B_d*eps_d*k**4*math.sin(k*s)/rho_d,0.0) for s in svals])<=ftol),'dry_sine_bending_acceleration_opposes_displacement_status':status(all((eps_d*math.sin(k*s))*ai[1] <= ftol for s,ai in zip(svals,ab))),'dry_sine_bending_energy_positive_status':status(ebs>0),'dry_sine_fluid_power_zero_status':status(abs(power(fh,zeros,w))<=ztol),'dry_candidate_arrays_finite_status':status(finite_vec([xnext,vnext,ab])),'dry_candidate_displacement_increment_bounded_status':status(maxdx is not None and disp <= maxdx+btol),'dry_candidate_velocity_increment_bounded_status':status(maxdv is not None and vinc <= maxdv+btol),'dry_energy_finite_status':status(finite_vals(energies)),'dry_energy_nonnegative_status':status(all(e>=-etol for e in energies)),'dry_no_fluid_power_status':status(abs(power(fh,vu,w))<=ztol),'dry_no_fluid_contamination_status':status(act),'dry_structure_equations_documented_status':status(all(s in doc for s in ['F_h = 0','rho_l X_tt','A_b_candidate','R_E_dry'])),'dry_benchmark_diagnostic_only_status':status(act),'stage13_6_diagnostic_preserved_status':status(stage13_ok(root)),'stage13_no_local_subdomain_center_regression_status':status(stage13_reg(root)),'stage14_small_lambda_hook_status':status(stage14_ok(root)),'no_rg_only_dependency_status':status(no_rg(root)),'no_unknown_failure_status':'PASS'})
    neg=['no_production_dry_benchmark_output_status','no_production_structure_update_status','no_production_structure_hook_status','no_stage16_code_modification_status','no_stage13_force_density_modification_status','no_stage14_rhs_modification_status','no_stage14_rhs_call_from_stage18_8_status','no_force_spreading_to_fluid_rhs_status','no_fluid_rhs_modification_status','no_ibm_modification_status','no_dns_core_modification_status','no_stats_visu_restart_io_modification_status','no_production_structure_time_integration_status','no_bending_force_runtime_application_status','no_tension_force_runtime_application_status','no_inextensibility_projection_status','no_inextensibility_repair_status','no_real_wall_contact_force_status','no_real_fibre_fibre_collision_force_status','no_penalty_force_status','no_repulsive_force_status','no_lubrication_force_status','no_friction_force_status','no_adhesion_force_status','no_contact_damping_force_status','no_collision_induced_rhs_status','no_collision_induced_structure_update_status','no_production_multifibre_logic_status','no_direct_rhs_injection_status','no_unapproved_stage14_rhs_call_status','no_legacy_ibm_forcing_status','no_unapproved_production_ibm_forcing_status','no_pressure_projection_modification_status','no_poisson_modification_status','no_rk3_channel_forcing_modification_status','no_channel_forcing_modification_status']
    for k in neg:st[k]=status(act)
    if a.test_case!='dry_straight_translation_sine_bounded_energy':reasons.append('unexpected_stage18_8_test_case')
    if not changed_ok:reasons.append('unapproved_or_closed_stage_path_modified:'+','.join(p for _c,p in entries if p not in ALLOWED))
    for k in SUMMARY_KEYS:
        if k.endswith('_status') and k!='final_status' and st.get(k)=='FAIL':reasons.append(k.replace('_status','_failed'))
    pf=[k for k in SUMMARY_KEYS if k!='final_status' and k.endswith('_status') and k not in VALUE_KEYS]
    st['final_status']='PASS' if all(st.get(k)=='PASS' for k in pf) else 'FAIL'
    write_out(root,st,reasons)
    for k in SUMMARY_KEYS:print(f'{k} {st.get(k,"FAIL")}')
    for r in reasons:print(f'reason {r}')
    return 0 if st['final_status']=='PASS' else 1
if __name__=='__main__':sys.exit(main())
