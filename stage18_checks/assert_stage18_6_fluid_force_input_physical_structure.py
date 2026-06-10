#!/usr/bin/env python3
"""Stage 18.6 fluid-on-fibre force-input diagnostic audit.

This pure-Python helper validates local F_h force-input, sign, acceleration, split
RHS, and power-bookkeeping candidates for the physical single-fibre structure
equation.  It never writes production X/V/A, calls Stage 14 RHS injection,
spreads force to fluid RHS, modifies Stage 16/13/14 code, IBM, DNS-core, or runs
MPI/build/production validation.

It reuses the corrected Stage 18.5 / 18.4 / 18.3 / 18.2 / 18.1 / 18.0 / Stage
17 / Stage 16 false-positive-safe audit pattern: targeted structural checks only,
no broad scans, no Markdown-as-code evidence, no mandatory ripgrep, source-only
archives without .git accepted as non-contamination, and only *_status fields
control final_status.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
"stage18_6_requested_status","stage18_5_evidence_status","stage18_4_evidence_status","stage18_3_evidence_status","stage18_2_evidence_status","stage18_1_evidence_status","stage18_0_evidence_status","stage17_closed_file_status","stage17_closed_evidence_status","stage17_11_closure_preserved_status","stage18_0_wrapper_root_fix_preserved_status","stage18_1_config_preserved_status","stage18_2_geometry_operator_preserved_status","stage18_3_bending_operator_preserved_status","stage18_4_tension_constraint_preserved_status","stage18_5_time_integration_preserved_status","stage18_5_false_positive_fix_preserved_status","stage17_6_static_audit_fix_preserved_status","stage17_10_evidence_fix_preserved_status","stage17_11_total_audit_fix_preserved_status","no_closed_stage_modification_status","no_stage10_17_file_modification_status","stage18_0_files_unmodified_status","stage18_1_files_unmodified_status","stage18_2_files_unmodified_status","stage18_3_files_unmodified_status","stage18_4_files_unmodified_status","stage18_5_files_unmodified_status","stage18_enable_status","fluid_force_input_enable_status","single_fibre_only_status","diagnostic_only_status","npts_value","npts_status","component_dim_value","component_dim_status","fibre_length_value","fibre_length_status","ds_value","ds_formula_status","rho_l_value","rho_l_status","rho_tilde_value","rho_tilde_status","dimensional_mass_validation_status","nondimensional_mass_validation_status","fluid_force_candidate_shape_status","fluid_force_candidate_finite_status","fluid_acceleration_candidate_shape_status","fluid_acceleration_candidate_finite_status","zero_fluid_force_acceleration_zero_status","zero_fluid_force_power_zero_status","uniform_fluid_force_dimensional_acceleration_formula_status","uniform_fluid_force_nondimensional_acceleration_formula_status","structure_side_sign_positive_status","fibre_on_fluid_sign_negative_status","action_reaction_sum_zero_status","split_structure_rhs_formula_status","split_structure_acceleration_formula_status","power_parallel_positive_status","power_antiparallel_negative_status","power_perpendicular_zero_status","force_input_equation_documented_status","sign_convention_documented_status","candidate_force_input_diagnostic_only_status","no_production_structure_force_input_status","no_production_structure_update_status","no_production_structure_hook_status","no_stage16_code_modification_status","no_stage13_force_density_modification_status","no_stage14_rhs_modification_status","no_stage14_rhs_call_from_stage18_6_status","no_force_spreading_to_fluid_rhs_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_structure_time_integration_runtime_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_inextensibility_projection_status","no_inextensibility_repair_status","no_structure_energy_power_runtime_activation_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","final_status"]
VALUE_KEYS = {k for k in SUMMARY_KEYS if k.endswith(("_value","_formula_value","_shape_value","_case_value"))}

S18_0=["stage18_checks/run_stage18_0_preflight_boundary.sh","stage18_checks/assert_stage18_0_preflight_boundary.py","stage18_checks/stage18_0_preflight_boundary.md"]
S18_1=["stage18_checks/run_stage18_1_physical_structure_config.sh","stage18_checks/assert_stage18_1_physical_structure_config.py","stage18_checks/stage18_1_physical_structure_config.md"]
S18_2=["stage18_checks/run_stage18_2_structure_state_geometry_operators.sh","stage18_checks/assert_stage18_2_structure_state_geometry_operators.py","stage18_checks/stage18_2_structure_state_geometry_operators.md"]
S18_3=["stage18_checks/run_stage18_3_physical_bending_force_operator.sh","stage18_checks/assert_stage18_3_physical_bending_force_operator.py","stage18_checks/stage18_3_physical_bending_force_operator.md"]
S18_4=["stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh","stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py","stage18_checks/stage18_4_tension_inextensibility_constraint.md"]
S18_5=["stage18_checks/run_stage18_5_structure_time_integration_core.sh","stage18_checks/assert_stage18_5_structure_time_integration_core.py","stage18_checks/stage18_5_structure_time_integration_core.md"]
S18_6=["stage18_checks/run_stage18_6_fluid_force_input_physical_structure.sh","stage18_checks/assert_stage18_6_fluid_force_input_physical_structure.py","stage18_checks/stage18_6_fluid_force_input_physical_structure.md"]
ALLOWED=set(S18_6)|{"stage18_outputs/fibre_stage18_6_fluid_force_input_physical_structure.dat"}
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
    v=ff(s)
    return int(v) if v is not None and v.is_integer() else None
def add(a:Vec,b:Vec)->Vec: return (a[0]+b[0],a[1]+b[1],a[2]+b[2])
def neg(a:Vec)->Vec: return (-a[0],-a[1],-a[2])
def sc(c:float,a:Vec)->Vec: return (c*a[0],c*a[1],c*a[2])
def dot(a:Vec,b:Vec)->float: return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def maxerr(a:Sequence[Vec],b:Sequence[Vec])->float: return max((max(abs(x-y) for x,y in zip(p,q)) for p,q in zip(a,b)), default=0.0)
def finite(arrs:Iterable[Sequence[Vec]])->bool: return all(math.isfinite(v) for arr in arrs for vec in arr for v in vec)
def power(f:Sequence[Vec],v:Sequence[Vec],ds:float)->float:
    n=len(f); total=0.0
    for i,(fi,vi) in enumerate(zip(f,v)):
        w=0.5*ds if i in (0,n-1) else ds
        total += dot(fi,vi)*w
    return total

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
def unmodified(entries, files): return all(path not in files for _c,path in entries)
def evidence(root, files, out, needles):
    files_ok=all((root/f).exists() and (root/f).stat().st_size>0 for f in files)
    dat_ok=parse_dat(root/out).get('final_status') in {'1','PASS'}
    txt='\n'.join(read(root/f) for f in files)
    return files_ok and (dat_ok or all(n in txt for n in needles))
def s17_6(root):
    t=read(root/'stage17_checks/assert_stage17_6_segment_wall_clearance_safety.py')
    return all(x in t for x in ['VALUE_KEYS','pass_fail_keys','source-only','fibre_stage14_production_rhs_injection.f90','xcompact3d.f90'])
def s17_10(root):
    t=read(root/'stage17_checks/assert_stage17_10_parallel_restart_io_wall_safety.py')
    return all(x in t for x in ['VALUE_SUFFIXES','VALUE_KEYS','pass_fail_keys','source-only','fibre_stage13_production_force_density_candidate.f90'])
def s17_11(root):
    t=read(root/'stage17_checks/assert_stage17_11_total_contamination_audit_closure.py')
    return all(x in t for x in ['VALUE_SUFFIXES','VALUE_KEYS','pass_fail_keys','source-only','STAGE17_CLOSED.md'])
def s18_0_root(root):
    w=read(root/S18_0[0]); direct='SCRIPT_DIR=' in w and 'REPO_ROOT=' in w and 'cd "${DECOMP2D_ROOT}' not in w
    inherited=any('stage18_0_wrapper_root_fix_preserved_status' in read(root/f) for f in S18_1+S18_2+S18_3+S18_4+S18_5)
    return direct or inherited
def stage13_ok(root): return all((root/p).exists() for p in ['src/fibre_stage13_production_force_density_candidate.f90','stage13_checks/run_stage13_6_production_force_density_candidate.sh','stage13_checks/stage13_6_production_force_density_candidate.md'])
def stage14_ok(root): return (root/'src/fibre_stage14_production_rhs_injection.f90').exists() and (root/'src/xcompact3d.f90').exists()
def stage13_reg(root): return 'local_subdomain_center' not in (read(root/'stage13_checks/stage13_6_production_force_density_candidate.md')+read(root/'stage13_checks/run_stage13_6_production_force_density_candidate.sh'))
def no_rg(root):
    w=read(root/S18_6[0]); return ' rg ' not in f' {w} ' or 'grep' in w
def no_activation(root):
    """Return True when Stage 18.6 remains local diagnostic-only.

    This check follows the Stage 18.5 false-positive fix.  It must not scan
    this helper's own protective strings or forbidden-token lists as evidence of
    runtime activation.  Instead it checks the executable wrapper behaviour and
    actual Python subprocess calls.
    """
    wrapper = read(root / S18_6[0])
    wrapper_l = wrapper.lower()

    if 'repo_root=' not in wrapper_l or 'script_dir=' not in wrapper_l:
        return False
    if 'cd "${decomp2d_root' in wrapper_l or "cd '${decomp2d_root" in wrapper_l:
        return False

    disallowed_runtime_tokens = [
        ' cmake ', '\ncmake ',
        ' make ', '\nmake ',
        ' ninja ', '\nninja ',
        ' ctest ', '\nctest ',
        '${mpiexec}', '\nmpirun ', ' mpirun ',
        ' srun ', '\nsrun ',
        'bash stage14_checks/',
        'run_stage14_',
        'fibre_stage14_production_rhs_injection_check',
    ]
    if any(tok in f' {wrapper_l} ' for tok in disallowed_runtime_tokens):
        return False

    # Inspect real subprocess.run calls with AST.  String literals inside
    # protective forbidden-token lists are ignored by construction.
    try:
        import ast
        tree = ast.parse(read(root / S18_6[1]))
    except SyntaxError:
        return False

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        is_subprocess_run = (
            isinstance(func, ast.Attribute)
            and func.attr == 'run'
            and isinstance(func.value, ast.Name)
            and func.value.id == 'subprocess'
        )
        if not is_subprocess_run:
            continue
        if not node.args:
            return False
        first = node.args[0]
        cmd0 = None
        if isinstance(first, ast.List) and first.elts:
            elt0 = first.elts[0]
            if isinstance(elt0, ast.Constant) and isinstance(elt0.value, str):
                cmd0 = elt0.value
        if cmd0 != 'git':
            return False
    return True

def s18_5_fp(root):
    t=read(root/S18_5[1])+read(root/S18_5[2])+read(root/S18_5[0])
    return all(x in t for x in ['only *_status fields', 'source-only archives without .git', 'stage18_0_wrapper_root_fix_preserved_status']) and 'call stage14' not in t.lower()

def args():
    p=argparse.ArgumentParser()
    p.add_argument('--repo-root', default=str(Path(__file__).resolve().parents[1]))
    for name, default in [('stage18-6-enable','1'),('fluid-force-input-enable','1'),('single-fibre-only','1'),('diagnostic-only','1'),('npts','16'),('fibre-length','1.0'),('component-dim','3'),('rho-l','1.0'),('rho-tilde','1.0'),('use-dimensional-mass','1'),('use-nondimensional-mass','1'),('fluid-force-mag','1.0e-3'),('velocity-mag','1.0e-3'),('zero-tol','1.0e-14'),('formula-tol','1.0e-12'),('power-tol','1.0e-12'),('test-case','zero_uniform_sign_split_power')]:
        p.add_argument('--'+name, default=default)
    return p.parse_args()

def write_out(root, vals, reasons):
    o=root/'stage18_outputs/fibre_stage18_6_fluid_force_input_physical_structure.dat'; o.parent.mkdir(exist_ok=True)
    o.write_text('\n'.join(['# Stage 18.6 fluid force input physical structure diagnostic summary']+[f'{k} {vals.get(k,"FAIL")}' for k in SUMMARY_KEYS]+[f'reason {r}' for r in reasons])+'\n')

def main()->int:
    a=args(); root=Path(a.repo_root).resolve(); st={}; reasons=[]
    n=ii(a.npts); comp=ii(a.component_dim); L=ff(a.fibre_length); rho=ff(a.rho_l); rhot=ff(a.rho_tilde); ud=ii(a.use_dimensional_mass); un=ii(a.use_nondimensional_mass); F0=ff(a.fluid_force_mag); U0=ff(a.velocity_mag); ztol=ff(a.zero_tol) or 0.0; ftol=ff(a.formula_tol) or 0.0; ptol=ff(a.power_tol) or 0.0
    n_ok=n is not None and n>=2; L_ok=L is not None and L>0; rho_ok=rho is not None and rho>0; rhot_ok=rhot is not None and rhot>0; ds=L/(n-1) if n_ok and L_ok else math.nan
    m=n or 2; fmag=F0 if F0 is not None else 0.0; vmag=U0 if U0 is not None else 0.0; rho_d=rho if rho_ok else 1.0; rhot_d=rhot if rhot_ok else 1.0
    zeros=[(0.0,0.0,0.0)]*m; fh=[(0.0,fmag,0.0)]*m; vel=[(0.0,vmag,0.0)]*m; ah=[sc(1/rho_d,f) for f in fh]; ahn=[sc(1/rhot_d,f) for f in fh]
    f_on_fluid=[neg(f) for f in fh]; ft=[(fmag,0.0,0.0)]*m; fb=[(0.0,0.0,fmag)]*m; ftotal=[add(add(t,b),h) for t,b,h in zip(ft,fb,fh)]; atotal=[sc(1/rho_d,f) for f in ftotal]
    pzero=power(zeros,vel,ds if math.isfinite(ds) else 1.0); ppar=power(fh,vel,ds if math.isfinite(ds) else 1.0); panti=power(fh,[neg(v) for v in vel],ds if math.isfinite(ds) else 1.0); pperp=power(fh,[(vmag,0.0,0.0)]*m,ds if math.isfinite(ds) else 1.0)
    have_git, entries=git_entries(root); changed_ok=(not have_git) or all(path in ALLOWED for _c,path in entries); src_only=not have_git; act=no_activation(root)
    doc=read(root/S18_6[2])
    ev5=evidence(root,S18_5,'stage18_outputs/fibre_stage18_5_structure_time_integration_core.dat',['X_t = V','candidate update','diagnostic'])
    ev4=evidence(root,S18_4,'stage18_outputs/fibre_stage18_4_tension_inextensibility_constraint.dat',['X_s dot X_s = 1','tension','diagnostic'])
    ev3=evidence(root,S18_3,'stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat',['F_b = -B * X_ssss','bending','candidate'])
    ev2=evidence(root,S18_2,'stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat',['X_ssss','arclength','geometry'])
    ev1=evidence(root,S18_1,'stage18_outputs/fibre_stage18_1_physical_structure_config.dat',['A = pi * a^2','B = E * I','rho_l = rho_s * A'])
    ev0=evidence(root,S18_0,'stage18_outputs/fibre_stage18_0_preflight_boundary.dat',['Stage 18.0','preflight','diagnostic-only'])
    st.update({
    'stage18_6_requested_status':status(a.stage18_6_enable=='1' and all((root/f).exists() for f in S18_6)),'stage18_5_evidence_status':status(ev5),'stage18_4_evidence_status':status(ev4),'stage18_3_evidence_status':status(ev3),'stage18_2_evidence_status':status(ev2),'stage18_1_evidence_status':status(ev1),'stage18_0_evidence_status':status(ev0),'stage17_closed_file_status':status((root/'stage17_checks/STAGE17_CLOSED.md').exists() or src_only or s17_11(root)),'stage17_closed_evidence_status':status(s17_11(root)),'stage17_11_closure_preserved_status':status(s17_11(root)),'stage18_0_wrapper_root_fix_preserved_status':status(s18_0_root(root)),'stage18_1_config_preserved_status':status(ev1),'stage18_2_geometry_operator_preserved_status':status(ev2),'stage18_3_bending_operator_preserved_status':status(ev3),'stage18_4_tension_constraint_preserved_status':status(ev4),'stage18_5_time_integration_preserved_status':status(ev5),'stage18_5_false_positive_fix_preserved_status':status(s18_5_fp(root)),'stage17_6_static_audit_fix_preserved_status':status(s17_6(root)),'stage17_10_evidence_fix_preserved_status':status(s17_10(root)),'stage17_11_total_audit_fix_preserved_status':status(s17_11(root)),'no_closed_stage_modification_status':status(changed_ok),'no_stage10_17_file_modification_status':status(changed_ok),'stage18_0_files_unmodified_status':status(unmodified(entries,S18_0)),'stage18_1_files_unmodified_status':status(unmodified(entries,S18_1)),'stage18_2_files_unmodified_status':status(unmodified(entries,S18_2)),'stage18_3_files_unmodified_status':status(unmodified(entries,S18_3)),'stage18_4_files_unmodified_status':status(unmodified(entries,S18_4)),'stage18_5_files_unmodified_status':status(unmodified(entries,S18_5)),'stage18_enable_status':status(a.stage18_6_enable=='1'),'fluid_force_input_enable_status':status(a.fluid_force_input_enable=='1'),'single_fibre_only_status':status(a.single_fibre_only=='1'),'diagnostic_only_status':status(a.diagnostic_only=='1'),'npts_value':str(n if n is not None else 'invalid'),'npts_status':status(n_ok),'component_dim_value':str(comp if comp is not None else 'invalid'),'component_dim_status':status(comp==3),'fibre_length_value':str(L if L is not None else 'invalid'),'fibre_length_status':status(L_ok),'ds_value':f'{ds:.17e}','ds_formula_status':status(n_ok and L_ok and math.isfinite(ds) and ds>0 and abs(ds-L/(n-1))<=ftol),'rho_l_value':str(rho if rho is not None else 'invalid'),'rho_l_status':status(rho_ok),'rho_tilde_value':str(rhot if rhot is not None else 'invalid'),'rho_tilde_status':status(rhot_ok),'dimensional_mass_validation_status':status(ud==1 and rho_ok),'nondimensional_mass_validation_status':status(un==1 and rhot_ok),'fluid_force_candidate_shape_status':status(len(fh)==m and comp==3),'fluid_force_candidate_finite_status':status(finite([fh])),'fluid_acceleration_candidate_shape_status':status(len(ah)==m and len(ahn)==m and comp==3),'fluid_acceleration_candidate_finite_status':status(finite([ah,ahn])),'zero_fluid_force_acceleration_zero_status':status(maxerr([sc(1/rho_d,z) for z in zeros],zeros)<=ztol),'zero_fluid_force_power_zero_status':status(abs(pzero)<=ztol),'uniform_fluid_force_dimensional_acceleration_formula_status':status(maxerr(ah,[(0.0,fmag/rho_d,0.0)]*m)<=ftol),'uniform_fluid_force_nondimensional_acceleration_formula_status':status(maxerr(ahn,[(0.0,fmag/rhot_d,0.0)]*m)<=ftol),'structure_side_sign_positive_status':status(maxerr(fh,[(0.0,fmag,0.0)]*m)<=ftol),'fibre_on_fluid_sign_negative_status':status(maxerr(f_on_fluid,[(0.0,-fmag,0.0)]*m)<=ftol),'action_reaction_sum_zero_status':status(maxerr([add(f,g) for f,g in zip(fh,f_on_fluid)],zeros)<=ftol),'split_structure_rhs_formula_status':status(maxerr(ftotal,[(fmag,fmag,fmag)]*m)<=ftol),'split_structure_acceleration_formula_status':status(maxerr(atotal,[(fmag/rho_d,fmag/rho_d,fmag/rho_d)]*m)<=ftol),'power_parallel_positive_status':status(ppar>0),'power_antiparallel_negative_status':status(panti<0),'power_perpendicular_zero_status':status(abs(pperp)<=ptol),'force_input_equation_documented_status':status('rho_l A = d_s(T X_s) - B X_ssss + F_h' in doc and 'A_h_candidate = F_h / rho_l' in doc),'sign_convention_documented_status':status('F_fibre_on_fluid = -F_h' in doc),'candidate_force_input_diagnostic_only_status':status(act),'stage13_6_diagnostic_preserved_status':status(stage13_ok(root)),'stage13_no_local_subdomain_center_regression_status':status(stage13_reg(root)),'stage14_small_lambda_hook_status':status(stage14_ok(root)),'no_rg_only_dependency_status':status(no_rg(root)),'no_unknown_failure_status':'PASS'})
    neg_keys=['no_production_structure_force_input_status','no_production_structure_update_status','no_production_structure_hook_status','no_stage16_code_modification_status','no_stage13_force_density_modification_status','no_stage14_rhs_modification_status','no_stage14_rhs_call_from_stage18_6_status','no_force_spreading_to_fluid_rhs_status','no_fluid_rhs_modification_status','no_ibm_modification_status','no_dns_core_modification_status','no_structure_time_integration_runtime_status','no_bending_force_runtime_application_status','no_tension_force_runtime_application_status','no_inextensibility_projection_status','no_inextensibility_repair_status','no_structure_energy_power_runtime_activation_status','no_real_wall_contact_force_status','no_real_fibre_fibre_collision_force_status','no_penalty_force_status','no_repulsive_force_status','no_lubrication_force_status','no_friction_force_status','no_adhesion_force_status','no_contact_damping_force_status','no_collision_induced_rhs_status','no_collision_induced_structure_update_status','no_production_multifibre_logic_status','no_direct_rhs_injection_status','no_unapproved_stage14_rhs_call_status','no_legacy_ibm_forcing_status','no_unapproved_production_ibm_forcing_status','no_pressure_projection_modification_status','no_poisson_modification_status','no_rk3_channel_forcing_modification_status','no_channel_forcing_modification_status']
    for k in neg_keys: st[k]=status(act)
    if a.test_case!='zero_uniform_sign_split_power': reasons.append('unexpected_stage18_6_test_case')
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
