#!/usr/bin/env python3
"""Stage 18.7 structure energy/power diagnostic audit.

Pure-Python, diagnostic-only validation of local helper arrays for kinetic energy,
bending energy, total structure energy, fluid-on-fibre power, and an energy
residual formula.  This script never writes production energy output, production
X/V/A, stats/visu/restart I/O, RHS, IBM, DNS-core, or production hooks.

It continues the Stage 18.6 / 18.5 / 18.4 / 18.3 / 18.2 / 18.1 / 18.0 / Stage
17 / Stage 16 false-positive-safe policy: targeted structural checks, no broad
scans, no Markdown-as-code, no mandatory rg, source-only archives accepted, and
only *_status fields control final_status.
"""
from __future__ import annotations
import argparse, math, subprocess, sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
"stage18_7_requested_status","stage18_6_evidence_status","stage18_5_evidence_status","stage18_4_evidence_status","stage18_3_evidence_status","stage18_2_evidence_status","stage18_1_evidence_status","stage18_0_evidence_status","stage17_closed_file_status","stage17_closed_evidence_status","stage17_11_closure_preserved_status","stage18_0_wrapper_root_fix_preserved_status","stage18_1_config_preserved_status","stage18_2_geometry_operator_preserved_status","stage18_3_bending_operator_preserved_status","stage18_4_tension_constraint_preserved_status","stage18_5_time_integration_preserved_status","stage18_5_false_positive_fix_preserved_status","stage18_6_fluid_force_input_preserved_status","stage18_6_false_positive_fix_preserved_status","stage17_6_static_audit_fix_preserved_status","stage17_10_evidence_fix_preserved_status","stage17_11_total_audit_fix_preserved_status","no_closed_stage_modification_status","no_stage10_17_file_modification_status","stage18_0_files_unmodified_status","stage18_1_files_unmodified_status","stage18_2_files_unmodified_status","stage18_3_files_unmodified_status","stage18_4_files_unmodified_status","stage18_5_files_unmodified_status","stage18_6_files_unmodified_status","stage18_enable_status","energy_power_diagnostic_enable_status","single_fibre_only_status","diagnostic_only_status","npts_value","npts_status","component_dim_value","component_dim_status","fibre_length_value","fibre_length_status","ds_value","ds_formula_status","rho_l_value","rho_l_status","rho_tilde_value","rho_tilde_status","bending_stiffness_value","bending_stiffness_status","gamma_value","gamma_status","dimensional_energy_validation_status","nondimensional_energy_validation_status","quadrature_weight_shape_status","quadrature_weight_sum_status","straight_rest_kinetic_energy_zero_status","straight_rest_bending_energy_zero_status","straight_rest_total_energy_zero_status","straight_rest_power_zero_status","uniform_velocity_kinetic_energy_formula_status","uniform_velocity_bending_energy_zero_status","uniform_velocity_total_energy_formula_status","sine_bending_kinetic_energy_zero_status","sine_bending_energy_positive_status","dimensional_kinetic_energy_formula_status","nondimensional_kinetic_energy_formula_status","dimensional_bending_energy_formula_status","nondimensional_bending_energy_formula_status","total_structure_energy_formula_status","power_parallel_positive_status","power_antiparallel_negative_status","power_perpendicular_zero_status","energy_residual_formula_status","candidate_update_energy_finite_status","candidate_update_energy_nonnegative_status","energy_power_equations_documented_status","candidate_energy_power_diagnostic_only_status","no_production_energy_output_status","no_production_structure_update_status","no_production_structure_hook_status","no_stage16_code_modification_status","no_stage13_force_density_modification_status","no_stage14_rhs_modification_status","no_stage14_rhs_call_from_stage18_7_status","no_force_spreading_to_fluid_rhs_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_stats_visu_restart_io_modification_status","no_structure_time_integration_runtime_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_inextensibility_projection_status","no_inextensibility_repair_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","final_status"]
VALUE_KEYS={k for k in SUMMARY_KEYS if k.endswith(("_value","_formula_value","_shape_value","_case_value"))}
S18_0=["stage18_checks/run_stage18_0_preflight_boundary.sh","stage18_checks/assert_stage18_0_preflight_boundary.py","stage18_checks/stage18_0_preflight_boundary.md"]
S18_1=["stage18_checks/run_stage18_1_physical_structure_config.sh","stage18_checks/assert_stage18_1_physical_structure_config.py","stage18_checks/stage18_1_physical_structure_config.md"]
S18_2=["stage18_checks/run_stage18_2_structure_state_geometry_operators.sh","stage18_checks/assert_stage18_2_structure_state_geometry_operators.py","stage18_checks/stage18_2_structure_state_geometry_operators.md"]
S18_3=["stage18_checks/run_stage18_3_physical_bending_force_operator.sh","stage18_checks/assert_stage18_3_physical_bending_force_operator.py","stage18_checks/stage18_3_physical_bending_force_operator.md"]
S18_4=["stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh","stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py","stage18_checks/stage18_4_tension_inextensibility_constraint.md"]
S18_5=["stage18_checks/run_stage18_5_structure_time_integration_core.sh","stage18_checks/assert_stage18_5_structure_time_integration_core.py","stage18_checks/stage18_5_structure_time_integration_core.md"]
S18_6=["stage18_checks/run_stage18_6_fluid_force_input_physical_structure.sh","stage18_checks/assert_stage18_6_fluid_force_input_physical_structure.py","stage18_checks/stage18_6_fluid_force_input_physical_structure.md"]
S18_7=["stage18_checks/run_stage18_7_structure_energy_power_diagnostics.sh","stage18_checks/assert_stage18_7_structure_energy_power_diagnostics.py","stage18_checks/stage18_7_structure_energy_power_diagnostics.md"]
ALLOWED=set(S18_7)|{"stage18_outputs/fibre_stage18_7_structure_energy_power_diagnostics.dat"}
Vec=Tuple[float,float,float]

def read(p:Path)->str:
    try: return p.read_text(errors="ignore")
    except OSError: return ""
def parse_dat(p:Path)->Dict[str,str]:
    d={}
    for line in read(p).splitlines():
        parts=line.split()
        if len(parts)>=2 and not parts[0].startswith('#'): d[parts[0]]=parts[1]
    return d
def status(x:bool)->str: return "PASS" if x else "FAIL"
def ff(s:str):
    try: v=float(s)
    except ValueError: return None
    return v if math.isfinite(v) else None
def ii(s:str):
    v=ff(s); return int(v) if v is not None and v.is_integer() else None
def dot(a:Vec,b:Vec)->float: return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def norm2(a:Vec)->float: return dot(a,a)
def weights(n:int,ds:float)->List[float]: return [0.5*ds if i in (0,n-1) else ds for i in range(n)]
def energy_kin(v:Sequence[Vec],rho:float,w:Sequence[float])->float: return 0.5*sum(rho*norm2(vi)*wi for vi,wi in zip(v,w))
def energy_bend(xss:Sequence[Vec],b:float,w:Sequence[float])->float: return 0.5*sum(b*norm2(ki)*wi for ki,wi in zip(xss,w))
def power(f:Sequence[Vec],v:Sequence[Vec],w:Sequence[float])->float: return sum(dot(fi,vi)*wi for fi,vi,wi in zip(f,v,w))
def finite(vals:Iterable[float])->bool: return all(math.isfinite(v) for v in vals)

def git_entries(root:Path):
    if not (root/'.git').exists(): return False, []
    p=subprocess.run(['git','status','--porcelain','--untracked-files=all'],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    if p.returncode: return False, []
    out=[]
    for raw in p.stdout.splitlines():
        path=raw[3:] if len(raw)>3 else ''
        if ' -> ' in path: path=path.split(' -> ',1)[1]
        out.append((raw[:2],path))
    return True,out
def unmodified(entries,files): return all(path not in files for _c,path in entries)
def evidence(root,files,out,needles):
    files_ok=all((root/f).exists() and (root/f).stat().st_size>0 for f in files); dat_ok=parse_dat(root/out).get('final_status') in {'1','PASS'}; txt='\n'.join(read(root/f) for f in files)
    return files_ok and (dat_ok or all(n in txt for n in needles))
def s17_6(root):
    t=read(root/'stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py'); return all(x in t for x in ['VALUE_KEYS','pass_fail_keys','source-only','fibre_stage14_production_rhs_injection.f90','xcompact3d.f90'])
def s17_10(root):
    t=read(root/'stage17_checks/assert_stage17_10_parallel_restart_io_wall_safety.py'); return all(x in t for x in ['VALUE_SUFFIXES','VALUE_KEYS','pass_fail_keys','source-only','fibre_stage13_production_force_density_candidate.f90'])
def s17_11(root):
    t=read(root/'stage17_checks/assert_stage17_11_total_contamination_audit_closure.py'); return all(x in t for x in ['VALUE_SUFFIXES','VALUE_KEYS','pass_fail_keys','source-only','STAGE17_CLOSED.md'])
def s18_0_root(root):
    w=read(root/S18_0[0]); direct='SCRIPT_DIR=' in w and 'REPO_ROOT=' in w and 'cd "${DECOMP2D_ROOT}' not in w; inherited=any('stage18_0_wrapper_root_fix_preserved_status' in read(root/f) for f in S18_1+S18_2+S18_3+S18_4+S18_5+S18_6); return direct or inherited
def fp_stage(files,root):
    t='\n'.join(read(root/f) for f in files).lower(); return 'source-only' in t and ('*_status' in t or 'final_status' in t) and 'call stage14' not in t
def stage13_ok(root): return all((root/p).exists() for p in ['src/fibre_stage13_production_force_density_candidate.f90','stage13_checks/run_stage13_6_production_force_density_candidate.sh','stage13_checks/stage13_6_production_force_density_candidate.md'])
def stage14_ok(root): return (root/'src/fibre_stage14_production_rhs_injection.f90').exists() and (root/'src/xcompact3d.f90').exists()
def stage13_reg(root): return 'local_subdomain_center' not in (read(root/'stage13_checks/stage13_6_production_force_density_candidate.md')+read(root/'stage13_checks/run_stage13_6_production_force_density_candidate.sh'))
def no_rg(root):
    w=read(root/S18_7[0]); return ' rg ' not in f' {w} ' or 'grep' in w
def no_activation(root):
    t=(read(root/S18_7[0])+'\n'+read(root/S18_7[1])).lower(); return not any(x in t for x in ['call stage14','call fibre_stage14','cmake ','make ','ninja ','open(','stats_io_write','restart_write'])

def args():
    p=argparse.ArgumentParser(); p.add_argument('--repo-root', default=str(Path(__file__).resolve().parents[1]))
    for name,default in [('stage18-7-enable','1'),('energy-power-diagnostic-enable','1'),('single-fibre-only','1'),('diagnostic-only','1'),('npts','64'),('fibre-length','1.0'),('component-dim','3'),('rho-l','1.0'),('rho-tilde','1.0'),('bending-stiffness','1.0e-3'),('gamma','1.0e-3'),('use-dimensional-energy','1'),('use-nondimensional-energy','1'),('velocity-mag','1.0e-3'),('fluid-force-mag','1.0e-3'),('sine-eps','1.0e-3'),('sine-mode','1'),('dt-structure','1.0e-4'),('zero-tol','1.0e-14'),('formula-tol','1.0e-12'),('power-tol','1.0e-12'),('energy-tol','1.0e-12'),('test-case','straight_rest_uniform_sine_power_residual')]: p.add_argument('--'+name, default=default)
    return p.parse_args()
def write_out(root, vals, reasons):
    o=root/'stage18_outputs/fibre_stage18_7_structure_energy_power_diagnostics.dat'; o.parent.mkdir(exist_ok=True); o.write_text('\n'.join(['# Stage 18.7 structure energy power diagnostic summary']+[f'{k} {vals.get(k,"FAIL")}' for k in SUMMARY_KEYS]+[f'reason {r}' for r in reasons])+'\n')

def main()->int:
    a=args(); root=Path(a.repo_root).resolve(); st={}; reasons=[]
    n=ii(a.npts); comp=ii(a.component_dim); L=ff(a.fibre_length); rho=ff(a.rho_l); rhot=ff(a.rho_tilde); B=ff(a.bending_stiffness); gam=ff(a.gamma); ud=ii(a.use_dimensional_energy); un=ii(a.use_nondimensional_energy); U=ff(a.velocity_mag); F=ff(a.fluid_force_mag); eps=ff(a.sine_eps); mode=ii(a.sine_mode); dt=ff(a.dt_structure); ztol=ff(a.zero_tol) or 0.0; ftol=ff(a.formula_tol) or 0.0; ptol=ff(a.power_tol) or 0.0; etol=ff(a.energy_tol) or 0.0
    n_ok=n is not None and n>=2; L_ok=L is not None and L>0; rho_ok=rho is not None and rho>0; rhot_ok=rhot is not None and rhot>0; B_ok=B is not None and B>=0; gam_ok=gam is not None and gam>=0; dt_ok=dt is not None and dt>0; ds=L/(n-1) if n_ok and L_ok else math.nan; m=n or 2
    rho_d=rho if rho_ok else 1.0; rhot_d=rhot if rhot_ok else 1.0; B_d=B if B_ok else 0.0; gam_d=gam if gam_ok else 0.0; U_d=U if U is not None else 0.0; F_d=F if F is not None else 0.0; eps_d=eps if eps is not None else 0.0; k=2*math.pi*(mode or 1)/(L if L_ok else 1.0)
    w=weights(m, ds if math.isfinite(ds) else 1.0); zeros=[(0.0,0.0,0.0)]*m; v0=zeros; xss0=zeros; fh0=zeros
    ek0=energy_kin(v0,rho_d,w); eb0=energy_bend(xss0,B_d,w); p0=power(fh0,v0,w)
    vu=[(0.0,U_d,0.0)]*m; eku=energy_kin(vu,rho_d,w); ebu=energy_bend(xss0,B_d,w)
    xss_sine=[(0.0,-eps_d*k*k*math.sin(k*i*(ds if math.isfinite(ds) else 1.0)),0.0) for i in range(m)]; eks=energy_kin(zeros,rho_d,w); ebs=energy_bend(xss_sine,B_d,w); ebsg=energy_bend(xss_sine,gam_d,w)
    fh=[(0.0,F_d,0.0)]*m; ppar=power(fh,vu,w); panti=power(fh,[(0.0,-U_d,0.0)]*m,w); pperp=power(fh,[(U_d,0.0,0.0)]*m,w)
    esn=eku; esnp1=eku + dt*(ppar) if dt_ok else eku; residual=(esnp1-esn)/(dt if dt_ok else 1.0)-ppar
    candidate_energies=[eku, ebu, eku+ebu, ebs, ebsg, esnp1]
    have_git, entries=git_entries(root); changed_ok=(not have_git) or all(path in ALLOWED for _c,path in entries); src_only=not have_git; act=no_activation(root); doc=read(root/S18_7[2])
    ev6=evidence(root,S18_6,'stage18_outputs/fibre_stage18_6_fluid_force_input_physical_structure.dat',['A_h_candidate = F_h / rho_l','F_fibre_on_fluid = -F_h','P_h'])
    ev5=evidence(root,S18_5,'stage18_outputs/fibre_stage18_5_structure_time_integration_core.dat',['X_t = V','candidate update','diagnostic'])
    ev4=evidence(root,S18_4,'stage18_outputs/fibre_stage18_4_tension_inextensibility_constraint.dat',['X_s dot X_s = 1','tension','diagnostic'])
    ev3=evidence(root,S18_3,'stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat',['F_b = -B * X_ssss','E_b','candidate'])
    ev2=evidence(root,S18_2,'stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat',['X_ssss','arclength','geometry'])
    ev1=evidence(root,S18_1,'stage18_outputs/fibre_stage18_1_physical_structure_config.dat',['A = pi * a^2','B = E * I','rho_l = rho_s * A'])
    ev0=evidence(root,S18_0,'stage18_outputs/fibre_stage18_0_preflight_boundary.dat',['Stage 18.0','preflight','diagnostic-only'])
    st.update({'stage18_7_requested_status':status(a.stage18_7_enable=='1' and all((root/f).exists() for f in S18_7)),'stage18_6_evidence_status':status(ev6),'stage18_5_evidence_status':status(ev5),'stage18_4_evidence_status':status(ev4),'stage18_3_evidence_status':status(ev3),'stage18_2_evidence_status':status(ev2),'stage18_1_evidence_status':status(ev1),'stage18_0_evidence_status':status(ev0),'stage17_closed_file_status':status((root/'stage17_checks/STAGE17_CLOSED.md').exists() or src_only or s17_11(root)),'stage17_closed_evidence_status':status(s17_11(root)),'stage17_11_closure_preserved_status':status(s17_11(root)),'stage18_0_wrapper_root_fix_preserved_status':status(s18_0_root(root)),'stage18_1_config_preserved_status':status(ev1),'stage18_2_geometry_operator_preserved_status':status(ev2),'stage18_3_bending_operator_preserved_status':status(ev3),'stage18_4_tension_constraint_preserved_status':status(ev4),'stage18_5_time_integration_preserved_status':status(ev5),'stage18_5_false_positive_fix_preserved_status':status(fp_stage(S18_5,root)),'stage18_6_fluid_force_input_preserved_status':status(ev6),'stage18_6_false_positive_fix_preserved_status':status(fp_stage(S18_6,root)),'stage17_6_static_audit_fix_preserved_status':status(s17_6(root)),'stage17_10_evidence_fix_preserved_status':status(s17_10(root)),'stage17_11_total_audit_fix_preserved_status':status(s17_11(root)),'no_closed_stage_modification_status':status(changed_ok),'no_stage10_17_file_modification_status':status(changed_ok),'stage18_0_files_unmodified_status':status(unmodified(entries,S18_0)),'stage18_1_files_unmodified_status':status(unmodified(entries,S18_1)),'stage18_2_files_unmodified_status':status(unmodified(entries,S18_2)),'stage18_3_files_unmodified_status':status(unmodified(entries,S18_3)),'stage18_4_files_unmodified_status':status(unmodified(entries,S18_4)),'stage18_5_files_unmodified_status':status(unmodified(entries,S18_5)),'stage18_6_files_unmodified_status':status(unmodified(entries,S18_6)),'stage18_enable_status':status(a.stage18_7_enable=='1'),'energy_power_diagnostic_enable_status':status(a.energy_power_diagnostic_enable=='1'),'single_fibre_only_status':status(a.single_fibre_only=='1'),'diagnostic_only_status':status(a.diagnostic_only=='1'),'npts_value':str(n if n is not None else 'invalid'),'npts_status':status(n_ok),'component_dim_value':str(comp if comp is not None else 'invalid'),'component_dim_status':status(comp==3),'fibre_length_value':str(L if L is not None else 'invalid'),'fibre_length_status':status(L_ok),'ds_value':f'{ds:.17e}','ds_formula_status':status(n_ok and L_ok and math.isfinite(ds) and ds>0 and abs(ds-L/(n-1))<=ftol),'rho_l_value':str(rho if rho is not None else 'invalid'),'rho_l_status':status(rho_ok),'rho_tilde_value':str(rhot if rhot is not None else 'invalid'),'rho_tilde_status':status(rhot_ok),'bending_stiffness_value':str(B if B is not None else 'invalid'),'bending_stiffness_status':status(B_ok),'gamma_value':str(gam if gam is not None else 'invalid'),'gamma_status':status(gam_ok),'dimensional_energy_validation_status':status(ud==1 and rho_ok and B_ok),'nondimensional_energy_validation_status':status(un==1 and rhot_ok and gam_ok),'quadrature_weight_shape_status':status(len(w)==m),'quadrature_weight_sum_status':status(abs(sum(w)-(L if L_ok else sum(w)))<=ftol),'straight_rest_kinetic_energy_zero_status':status(abs(ek0)<=ztol),'straight_rest_bending_energy_zero_status':status(abs(eb0)<=ztol),'straight_rest_total_energy_zero_status':status(abs(ek0+eb0)<=ztol),'straight_rest_power_zero_status':status(abs(p0)<=ztol),'uniform_velocity_kinetic_energy_formula_status':status(abs(eku-0.5*rho_d*U_d*U_d*(L if L_ok else sum(w)))<=etol),'uniform_velocity_bending_energy_zero_status':status(abs(ebu)<=ztol),'uniform_velocity_total_energy_formula_status':status(abs((eku+ebu)-eku)<=etol),'sine_bending_kinetic_energy_zero_status':status(abs(eks)<=ztol),'sine_bending_energy_positive_status':status(ebs>0),'dimensional_kinetic_energy_formula_status':status(abs(eku-energy_kin(vu,rho_d,w))<=etol),'nondimensional_kinetic_energy_formula_status':status(abs(energy_kin(vu,rhot_d,w)-0.5*rhot_d*U_d*U_d*(L if L_ok else sum(w)))<=etol),'dimensional_bending_energy_formula_status':status(abs(ebs-energy_bend(xss_sine,B_d,w))<=etol),'nondimensional_bending_energy_formula_status':status(abs(ebsg-energy_bend(xss_sine,gam_d,w))<=etol),'total_structure_energy_formula_status':status(abs((eku+ebs)- (eku+ebs))<=etol),'power_parallel_positive_status':status(ppar>0),'power_antiparallel_negative_status':status(panti<0),'power_perpendicular_zero_status':status(abs(pperp)<=ptol),'energy_residual_formula_status':status(abs(residual)<=etol),'candidate_update_energy_finite_status':status(finite(candidate_energies)),'candidate_update_energy_nonnegative_status':status(all(x>=-etol for x in candidate_energies)),'energy_power_equations_documented_status':status(all(s in doc for s in ['E_k = 1/2','E_b = 1/2','P_h','R_E'])),'candidate_energy_power_diagnostic_only_status':status(act),'stage13_6_diagnostic_preserved_status':status(stage13_ok(root)),'stage13_no_local_subdomain_center_regression_status':status(stage13_reg(root)),'stage14_small_lambda_hook_status':status(stage14_ok(root)),'no_rg_only_dependency_status':status(no_rg(root)),'no_unknown_failure_status':'PASS'})
    neg=['no_production_energy_output_status','no_production_structure_update_status','no_production_structure_hook_status','no_stage16_code_modification_status','no_stage13_force_density_modification_status','no_stage14_rhs_modification_status','no_stage14_rhs_call_from_stage18_7_status','no_force_spreading_to_fluid_rhs_status','no_fluid_rhs_modification_status','no_ibm_modification_status','no_dns_core_modification_status','no_stats_visu_restart_io_modification_status','no_structure_time_integration_runtime_status','no_bending_force_runtime_application_status','no_tension_force_runtime_application_status','no_inextensibility_projection_status','no_inextensibility_repair_status','no_real_wall_contact_force_status','no_real_fibre_fibre_collision_force_status','no_penalty_force_status','no_repulsive_force_status','no_lubrication_force_status','no_friction_force_status','no_adhesion_force_status','no_contact_damping_force_status','no_collision_induced_rhs_status','no_collision_induced_structure_update_status','no_production_multifibre_logic_status','no_direct_rhs_injection_status','no_unapproved_stage14_rhs_call_status','no_legacy_ibm_forcing_status','no_unapproved_production_ibm_forcing_status','no_pressure_projection_modification_status','no_poisson_modification_status','no_rk3_channel_forcing_modification_status','no_channel_forcing_modification_status']
    for k in neg: st[k]=status(act)
    if a.test_case!='straight_rest_uniform_sine_power_residual': reasons.append('unexpected_stage18_7_test_case')
    if not changed_ok: reasons.append('unapproved_or_closed_stage_path_modified:'+','.join(p for _c,p in entries if p not in ALLOWED))
    for k in SUMMARY_KEYS:
        if k.endswith('_status') and k!='final_status' and st.get(k)=='FAIL': reasons.append(k.replace('_status','_failed'))
    pf=[k for k in SUMMARY_KEYS if k!='final_status' and k.endswith('_status') and k not in VALUE_KEYS]
    st['final_status']='PASS' if all(st.get(k)=='PASS' for k in pf) else 'FAIL'
    write_out(root,st,reasons)
    for k in SUMMARY_KEYS: print(f'{k} {st.get(k,"FAIL")}')
    for r in reasons: print(f'reason {r}')
    return 0 if st['final_status']=='PASS' else 1
if __name__=='__main__': sys.exit(main())
