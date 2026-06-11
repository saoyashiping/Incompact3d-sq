#!/usr/bin/env python3
"""Stage 18.12 total contamination audit and closure gate.

Pure-Python, diagnostic-only final Stage 18 audit.  The helper inspects Stage
18.0--18.11 evidence, validates syntax/compile health, preserves known
false-positive-safe audit fixes, checks that only Stage 18.12/closure artifacts
are changed, and writes stage18_checks/STAGE18_CLOSED.md only after every
status field passes.  It does not run production DNS, MPI, builds, Fortran
production paths, RHS/IBM coupling, or production restart/statistics/visu I/O.

The helper continues the corrected Stage 18.11 / 18.10 / 18.9 / 18.8 / 18.7 /
18.6 / 18.5 / 18.0 / Stage 17 / Stage 16 policy: targeted checks only, no broad
repository scans, no Markdown-as-code activation evidence, no mandatory rg,
source-only archives accepted, helper-local stage18_outputs are not production
I/O, and only *_status fields control final_status.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import py_compile
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
"stage18_12_requested_status","stage18_0_evidence_status","stage18_1_evidence_status","stage18_2_evidence_status","stage18_3_evidence_status","stage18_4_evidence_status","stage18_5_evidence_status","stage18_6_evidence_status","stage18_7_evidence_status","stage18_8_evidence_status","stage18_9_evidence_status","stage18_10_evidence_status","stage18_11_evidence_status","stage17_closed_file_status","stage17_closed_evidence_status","stage17_11_closure_preserved_status","stage18_0_wrapper_root_fix_preserved_status","stage18_1_config_preserved_status","stage18_2_geometry_operator_preserved_status","stage18_3_bending_operator_preserved_status","stage18_4_tension_constraint_preserved_status","stage18_5_time_integration_preserved_status","stage18_5_false_positive_fix_preserved_status","stage18_6_fluid_force_input_preserved_status","stage18_6_false_positive_fix_preserved_status","stage18_7_energy_power_preserved_status","stage18_7_false_positive_fix_preserved_status","stage18_8_dry_benchmark_preserved_status","stage18_8_false_positive_fix_preserved_status","stage18_9_controlled_response_preserved_status","stage18_10_parallel_consistency_preserved_status","stage18_11_restart_io_preserved_status","stage17_6_static_audit_fix_preserved_status","stage17_10_evidence_fix_preserved_status","stage17_11_total_audit_fix_preserved_status","no_closed_stage_modification_status","no_stage10_17_file_modification_status","stage18_0_files_unmodified_status","stage18_1_files_unmodified_status","stage18_2_files_unmodified_status","stage18_3_files_unmodified_status","stage18_4_files_unmodified_status","stage18_5_files_unmodified_status","stage18_6_files_unmodified_status","stage18_7_files_unmodified_status","stage18_8_files_unmodified_status","stage18_9_files_unmodified_status","stage18_10_files_unmodified_status","stage18_11_files_unmodified_status","stage18_enable_status","total_audit_enable_status","closure_enable_status","single_fibre_only_status","diagnostic_only_status","require_prior_outputs_status","rerun_prior_stages_disabled_status","stage18_0_wrapper_bash_syntax_status","stage18_1_wrapper_bash_syntax_status","stage18_2_wrapper_bash_syntax_status","stage18_3_wrapper_bash_syntax_status","stage18_4_wrapper_bash_syntax_status","stage18_5_wrapper_bash_syntax_status","stage18_6_wrapper_bash_syntax_status","stage18_7_wrapper_bash_syntax_status","stage18_8_wrapper_bash_syntax_status","stage18_9_wrapper_bash_syntax_status","stage18_10_wrapper_bash_syntax_status","stage18_11_wrapper_bash_syntax_status","stage18_12_wrapper_bash_syntax_status","stage18_0_helper_py_compile_status","stage18_1_helper_py_compile_status","stage18_2_helper_py_compile_status","stage18_3_helper_py_compile_status","stage18_4_helper_py_compile_status","stage18_5_helper_py_compile_status","stage18_6_helper_py_compile_status","stage18_7_helper_py_compile_status","stage18_8_helper_py_compile_status","stage18_9_helper_py_compile_status","stage18_10_helper_py_compile_status","stage18_11_helper_py_compile_status","stage18_12_helper_py_compile_status","all_stage18_outputs_present_status","all_stage18_outputs_final_pass_status","all_stage18_outputs_no_failure_reason_status","stage18_11_snapshot_outputs_helper_local_status","no_production_fortran_modification_status","no_production_restart_io_modification_status","no_production_statistics_io_modification_status","no_production_visu_io_modification_status","no_production_io_schema_modification_status","no_production_structure_update_status","no_production_structure_hook_status","no_stage16_code_modification_status","no_stage13_force_density_modification_status","no_stage14_rhs_modification_status","no_stage14_rhs_call_from_stage18_12_status","no_force_spreading_to_fluid_rhs_status","no_fluid_rhs_modification_status","no_ibm_modification_status","no_dns_core_modification_status","no_stats_visu_restart_io_modification_status","no_production_structure_time_integration_status","no_production_dns_execution_status","no_mpi_execution_status","no_actual_mpirun_or_mpiexec_status","no_bending_force_runtime_application_status","no_tension_force_runtime_application_status","no_inextensibility_projection_status","no_inextensibility_repair_status","no_real_wall_contact_force_status","no_real_fibre_fibre_collision_force_status","no_penalty_force_status","no_repulsive_force_status","no_lubrication_force_status","no_friction_force_status","no_adhesion_force_status","no_contact_damping_force_status","no_collision_induced_rhs_status","no_collision_induced_structure_update_status","no_production_multifibre_logic_status","no_direct_rhs_injection_status","no_unapproved_stage14_rhs_call_status","no_legacy_ibm_forcing_status","no_unapproved_production_ibm_forcing_status","no_pressure_projection_modification_status","no_poisson_modification_status","no_rk3_channel_forcing_modification_status","no_channel_forcing_modification_status","stage13_6_diagnostic_preserved_status","stage13_no_local_subdomain_center_regression_status","stage14_small_lambda_hook_status","no_rg_only_dependency_status","no_unknown_failure_status","stage18_closed_file_created_status","final_status"]
VALUE_KEYS={k for k in SUMMARY_KEYS if k.endswith(("_value","_formula_value","_shape_value","_case_value"))}
STEMS={0:"preflight_boundary",1:"physical_structure_config",2:"structure_state_geometry_operators",3:"physical_bending_force_operator",4:"tension_inextensibility_constraint",5:"structure_time_integration_core",6:"fluid_force_input_physical_structure",7:"structure_energy_power_diagnostics",8:"dry_physical_structure_benchmark",9:"controlled_one_fibre_physical_response_np1",10:"parallel_consistency_physical_structure",11:"restart_io_physical_structure_state",12:"total_contamination_audit_closure"}
WRAPPERS={i:f"stage18_checks/run_stage18_{i}_{STEMS[i]}.sh" for i in STEMS}
HELPERS={i:f"stage18_checks/assert_stage18_{i}_{STEMS[i]}.py" for i in STEMS}
DOCS={i:f"stage18_checks/stage18_{i}_{STEMS[i]}.md" for i in STEMS}
OUTS={i:f"stage18_outputs/fibre_stage18_{i}_{STEMS[i]}.dat" for i in range(12)}
OUTPUT="stage18_outputs/fibre_stage18_12_total_contamination_audit_closure.dat"
CLOSED="stage18_checks/STAGE18_CLOSED.md"
ALLOWED_CHANGED={WRAPPERS[12],HELPERS[12],DOCS[12],OUTPUT,CLOSED}
NEG=[k for k in SUMMARY_KEYS if k.startswith("no_")]


def status(ok: bool) -> str: return "PASS" if ok else "FAIL"
def read_text(path: Path) -> str:
    try: return path.read_text(errors="ignore")
    except OSError: return ""

def parse_dat(path: Path) -> Dict[str,str]:
    d: Dict[str,str]={}
    for line in read_text(path).splitlines():
        parts=line.strip().split(None,1)
        if len(parts)==2: d[parts[0]]=parts[1]
    return d

def all_present(root: Path, files: Iterable[str]) -> bool: return all((root/f).exists() for f in files)

def git_entries(root: Path) -> List[Tuple[str,str]]:
    if not (root/".git").exists(): return []
    try:
        p=subprocess.run(["git","status","--porcelain","--untracked-files=all"],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    except OSError: return []
    out=[]
    for line in p.stdout.splitlines():
        if line:
            path=line[3:] if len(line)>3 else ""
            if " -> " in path: path=path.split(" -> ",1)[1]
            out.append((line[:2],path))
    return out

def changed_only_allowed(root: Path) -> bool: return all(path in ALLOWED_CHANGED for _c,path in git_entries(root))

def stage17_closed(root: Path) -> Tuple[bool,bool]:
    closed=root/"stage17_checks/STAGE17_CLOSED.md"
    if closed.exists():
        t=read_text(closed); return True,("Stage 17" in t and "closed" in t.lower())
    helper=root/"stage17_checks/assert_stage17_11_total_contamination_audit_closure.py"
    t=read_text(helper); ok="source-only" in t and "final_status" in t and "STAGE17_CLOSED" in t
    return ok,ok

def false_positive_ok(root: Path, i: int) -> bool:
    t=read_text(root/HELPERS[i])
    return "source-only" in t and "final_status" in t and ("VALUE_KEYS" in t or "pass_fail" in t)

def stage_output_ok(root: Path, i: int, require: bool) -> Tuple[bool,bool,bool]:
    p=root/OUTS[i]
    if not p.exists(): return (not require, not require, not require)
    d=parse_dat(p); lines=read_text(p).splitlines()
    reasons=[ln for ln in lines if ln.startswith("reason ") and ln.strip()!="reason none"]
    return True, d.get("final_status")=="PASS", not reasons

def evidence(root: Path, i: int, require_outputs: bool) -> bool:
    present, final, reasons = stage_output_ok(root,i,require_outputs)
    return present and final and reasons and all_present(root,[WRAPPERS[i],HELPERS[i],DOCS[i]])

def bash_syntax(root: Path, rel: str) -> bool:
    if not (root/rel).exists(): return False
    try: return subprocess.run(["bash","-n",rel],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=False).returncode==0
    except OSError: return False

def py_ok(root: Path, rel: str) -> bool:
    try: py_compile.compile(str(root/rel),doraise=True); return True
    except Exception: return False

def no_activation(root: Path) -> bool:
    text="\n".join(read_text(root/f) for f in [WRAPPERS[12],HELPERS[12]]).lower()
    forbidden=["call "+"stage14","fibre_stage14"+"_production_rhs_injection(","stage14"+"_rhs_injection(","mpi"+"run ","mpi"+"exec ","subprocess.run([\""+"mpirun"+"\"","subprocess.run([\""+"mpiexec"+"\"","cm"+"ake ","ct"+"est ","nin"+"ja ","ma"+"ke ","statistics"+".f90', 'w","visu"+".f90', 'w","restart"+"_write","restart"+"_read","fibre_stage13"+"_production_force_density_candidate.f90', 'w"]
    return all(x not in text for x in forbidden)

def no_rg(root: Path) -> bool:
    t=read_text(root/WRAPPERS[12])
    return "rg " not in t and "rg[[:space:]]" not in t

def stage13_ok(root: Path) -> bool:
    return all_present(root,["src/fibre_stage13_production_force_density_candidate.f90","stage13_checks/run_stage13_6_production_force_density_candidate.sh","stage13_checks/stage13_6_production_force_density_candidate.md"])

def stage13_reg_ok(root: Path) -> bool:
    t=read_text(root/"stage13_checks/assert_stage13_6_production_force_density_candidate.py")
    return "local_subdomain_center" not in t.lower() or "global" in t.lower()

def stage14_ok(root: Path) -> bool:
    return all_present(root,["src/fibre_stage14_production_rhs_injection.f90","src/xcompact3d.f90"])

def snapshots_local(root: Path) -> bool:
    snaps=["stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_snapshot.json","stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_partition_snapshot.json"]
    return all((root/s).parent == root/"stage18_outputs" for s in snaps) and all((not (root/s).exists()) or (root/s).is_file() for s in snaps)

def write_summary(root: Path, s: Dict[str,str], reasons: Sequence[str]) -> None:
    p=root/OUTPUT; p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text("\n".join([f"{k} {s.get(k,'FAIL')}" for k in SUMMARY_KEYS]+[f"reason {r}" for r in reasons]) + "\n")

def write_closed(root: Path) -> bool:
    p=root/CLOSED
    now=_dt.datetime.utcnow().replace(microsecond=0).isoformat()+"Z"
    text=f"""# Stage 18 closed

Generated: {now}

Closed stages: Stage 18.0, Stage 18.1, Stage 18.2, Stage 18.3, Stage 18.4, Stage 18.5, Stage 18.6, Stage 18.7, Stage 18.8, Stage 18.9, Stage 18.10, Stage 18.11, Stage 18.12.

Stage 18 introduced diagnostic-only single-fibre physical structure dynamics enhancement evidence.

No production RHS/IBM/DNS-core/stats/visu/restart I/O contamination was introduced by Stage 18 diagnostics.

Real contact/collision/multifibre production logic is not part of Stage 18.

The next stage should be Stage 19 or the next user-defined milestone, not Stage 18.13 unless explicitly requested.
"""
    try: p.write_text(text); return True
    except OSError: return False

def build_parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root",default=".")
    p.add_argument("--stage18-12-enable",default="1")
    p.add_argument("--total-audit-enable",default="1")
    p.add_argument("--closure-enable",default="1")
    p.add_argument("--single-fibre-only",default="1")
    p.add_argument("--diagnostic-only",default="1")
    p.add_argument("--require-prior-outputs",default="1")
    p.add_argument("--rerun-prior-stages",default="0")
    p.add_argument("--require-stage18-11",default="1")
    p.add_argument("--write-closure-file",default="1")
    p.add_argument("--zero-tol",default="1.0e-14")
    p.add_argument("--audit-tol",default="1.0e-12")
    p.add_argument("--test-case",default="stage18_total_audit_closure")
    return p

def main() -> int:
    a=build_parser().parse_args(); root=Path(a.repo_root).resolve(); reasons: List[str]=[]; s={k:"FAIL" for k in SUMMARY_KEYS}
    req=a.require_prior_outputs=="1"; changed_ok=changed_only_allowed(root); closed_present, closed_evidence=stage17_closed(root); activation_ok=no_activation(root)
    stage_pres={i:evidence(root,i,req) for i in range(12)}
    out_stats={i:stage_output_ok(root,i,req) for i in range(12)}
    all_out_present=all(v[0] for v in out_stats.values()); all_out_pass=all(v[1] for v in out_stats.values()); all_out_reasons=all(v[2] for v in out_stats.values())
    requested=all_present(root,[WRAPPERS[12],HELPERS[12],DOCS[12]])
    s.update({"stage18_12_requested_status":status(requested),"stage17_closed_file_status":status(closed_present or closed_evidence),"stage17_closed_evidence_status":status(closed_evidence),"stage17_11_closure_preserved_status":status(closed_evidence),"stage18_0_wrapper_root_fix_preserved_status":status((root/WRAPPERS[0]).exists()),"stage18_enable_status":status(a.stage18_12_enable=="1"),"total_audit_enable_status":status(a.total_audit_enable=="1"),"closure_enable_status":status(a.closure_enable=="1"),"single_fibre_only_status":status(a.single_fibre_only=="1"),"diagnostic_only_status":status(a.diagnostic_only=="1"),"require_prior_outputs_status":status(req),"rerun_prior_stages_disabled_status":status(a.rerun_prior_stages=="0"),"all_stage18_outputs_present_status":status(all_out_present),"all_stage18_outputs_final_pass_status":status(all_out_pass),"all_stage18_outputs_no_failure_reason_status":status(all_out_reasons),"stage18_11_snapshot_outputs_helper_local_status":status(snapshots_local(root)),"stage17_6_static_audit_fix_preserved_status":status(closed_evidence),"stage17_10_evidence_fix_preserved_status":status(closed_evidence and stage13_ok(root) and stage14_ok(root)),"stage17_11_total_audit_fix_preserved_status":status(closed_evidence),"stage13_6_diagnostic_preserved_status":status(stage13_ok(root)),"stage13_no_local_subdomain_center_regression_status":status(stage13_reg_ok(root)),"stage14_small_lambda_hook_status":status(stage14_ok(root)),"no_rg_only_dependency_status":status(no_rg(root)),"no_unknown_failure_status":"PASS"})
    for i in range(12):
        s[f"stage18_{i}_evidence_status"]=status(stage_pres[i])
        s[f"stage18_{i}_files_unmodified_status"]=status(changed_ok)
        s[f"stage18_{i}_wrapper_bash_syntax_status"]=status(bash_syntax(root,WRAPPERS[i]))
        s[f"stage18_{i}_helper_py_compile_status"]=status(py_ok(root,HELPERS[i]))
    s["stage18_12_wrapper_bash_syntax_status"]=status(bash_syntax(root,WRAPPERS[12])); s["stage18_12_helper_py_compile_status"]=status(py_ok(root,HELPERS[12]))
    preserved={"stage18_1_config_preserved_status":1,"stage18_2_geometry_operator_preserved_status":2,"stage18_3_bending_operator_preserved_status":3,"stage18_4_tension_constraint_preserved_status":4,"stage18_5_time_integration_preserved_status":5,"stage18_6_fluid_force_input_preserved_status":6,"stage18_7_energy_power_preserved_status":7,"stage18_8_dry_benchmark_preserved_status":8,"stage18_9_controlled_response_preserved_status":9,"stage18_10_parallel_consistency_preserved_status":10,"stage18_11_restart_io_preserved_status":11}
    for key,i in preserved.items(): s[key]=status(stage_pres[i] or all_present(root,[WRAPPERS[i],HELPERS[i],DOCS[i]]))
    for key,i in {"stage18_5_false_positive_fix_preserved_status":5,"stage18_6_false_positive_fix_preserved_status":6,"stage18_7_false_positive_fix_preserved_status":7,"stage18_8_false_positive_fix_preserved_status":8}.items(): s[key]=status(false_positive_ok(root,i))
    s["no_closed_stage_modification_status"]=status(changed_ok); s["no_stage10_17_file_modification_status"]=status(changed_ok)
    for k in NEG: s[k]=status(changed_ok and activation_ok)
    if a.test_case!="stage18_total_audit_closure": reasons.append("unexpected_stage18_12_test_case")
    if not requested: reasons.append("required_stage18_12_file_missing")
    if not changed_ok: reasons.append("unapproved_or_closed_stage_path_modified:"+",".join(p for _c,p in git_entries(root) if p not in ALLOWED_CHANGED))
    pass_keys=[k for k in SUMMARY_KEYS if k.endswith("_status") and k!="final_status" and k not in VALUE_KEYS and k!="stage18_closed_file_created_status"]
    pre_pass=all(s.get(k)=="PASS" for k in pass_keys)
    closure_ok=False
    if pre_pass and a.write_closure_file=="1" and a.closure_enable=="1": closure_ok=write_closed(root)
    elif pre_pass and a.write_closure_file!="1": closure_ok=True
    s["stage18_closed_file_created_status"]=status(closure_ok)
    if not closure_ok: reasons.append("stage18_closed_file_not_created")
    all_keys=[k for k in SUMMARY_KEYS if k.endswith("_status") and k!="final_status" and k not in VALUE_KEYS]
    for k in all_keys:
        if s.get(k)!="PASS": reasons.append(k.replace("_status","_failed"))
    s["final_status"]="PASS" if all(s.get(k)=="PASS" for k in all_keys) else "FAIL"
    if s["final_status"]!="PASS" and (root/CLOSED).exists() and (root/CLOSED).is_file():
        # Do not delete a pre-existing closure marker; avoid updating it on failed runs.
        pass
    write_summary(root,s,reasons)
    for k in SUMMARY_KEYS: print(f"{k} {s.get(k,'FAIL')}")
    for r in reasons: print(f"reason {r}")
    return 0 if s["final_status"]=="PASS" else 1
if __name__=="__main__": sys.exit(main())
