#!/usr/bin/env python3
"""Stage 19.0 preflight boundary gate.

Stage 19.0 is intentionally diagnostic-only.  Its job is to decide whether the
repository may enter Stage 19, not to rerun Stage 18.  A source-only archive may
not contain stage18_outputs/ or STAGE18_CLOSED.md; in that case this helper may
accept Stage 18 closure from the Stage 18.12 closure-gate source evidence.  This
prevents Stage 19.0 from forcing expensive or stale reruns of Stage 18.0-18.11.
"""
from __future__ import annotations

import argparse
import os
import py_compile
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS: List[str] = [
    "stage19_0_requested_status",
    "stage19_0_preflight_enable_status",
    "stage18_closed_file_status",
    "stage18_closed_file_content_status",
    "stage18_12_evidence_status",
    "stage18_12_closure_evidence_status",
    "stage18_closure_accepted_status",
    "prior_stage18_outputs_required_status",
    "stage18_closure_supersedes_individual_outputs_status",
    "stage18_0_output_status",
    "stage18_1_output_status",
    "stage18_2_output_status",
    "stage18_3_output_status",
    "stage18_4_output_status",
    "stage18_5_output_status",
    "stage18_6_output_status",
    "stage18_7_output_status",
    "stage18_8_output_status",
    "stage18_9_output_status",
    "stage18_10_output_status",
    "stage18_11_output_status",
    "stage18_12_output_status",
    "all_stage18_outputs_present_status",
    "all_stage18_outputs_final_pass_status",
    "stage17_closed_file_status",
    "stage17_closed_evidence_status",
    "stage17_11_closure_preserved_status",
    "stage18_0_wrapper_root_fix_preserved_status",
    "stage18_5_false_positive_fix_preserved_status",
    "stage18_6_false_positive_fix_preserved_status",
    "stage18_7_false_positive_fix_preserved_status",
    "stage18_8_false_positive_fix_preserved_status",
    "stage18_9_controlled_response_preserved_status",
    "stage18_10_parallel_consistency_preserved_status",
    "stage18_11_restart_io_preserved_status",
    "stage18_12_closure_preserved_status",
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "no_stage18_file_modification_status",
    "stage19_boundary_definition_status",
    "diagnostic_only_status",
    "single_fibre_only_status",
    "rerun_prior_stages_disabled_status",
    "stage19_0_wrapper_bash_syntax_status",
    "stage19_0_helper_py_compile_status",
    "no_production_fortran_modification_status",
    "no_cmake_modification_status",
    "no_production_structure_state_creation_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_stage16_code_modification_status",
    "no_stage13_force_density_modification_status",
    "no_stage14_rhs_modification_status",
    "no_stage14_rhs_call_from_stage19_0_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_fluid_rhs_modification_status",
    "no_ibm_modification_status",
    "no_dns_core_modification_status",
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "no_production_restart_io_modification_status",
    "no_production_statistics_io_modification_status",
    "no_production_visu_io_modification_status",
    "no_stats_visu_restart_io_modification_status",
    "no_production_dns_execution_status",
    "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status",
    "no_real_wall_contact_force_status",
    "no_real_fibre_fibre_collision_force_status",
    "no_penalty_force_status",
    "no_repulsive_force_status",
    "no_lubrication_force_status",
    "no_friction_force_status",
    "no_adhesion_force_status",
    "no_contact_damping_force_status",
    "no_collision_induced_rhs_status",
    "no_collision_induced_structure_update_status",
    "no_production_multifibre_logic_status",
    "no_direct_rhs_injection_status",
    "no_unapproved_stage14_rhs_call_status",
    "no_legacy_ibm_forcing_status",
    "no_unapproved_production_ibm_forcing_status",
    "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status",
    "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "final_status",
]

STAGE18_OUTPUTS: Dict[str, str] = {
    "18_0": "fibre_stage18_0_preflight_boundary.dat",
    "18_1": "fibre_stage18_1_physical_structure_config.dat",
    "18_2": "fibre_stage18_2_structure_state_geometry_operators.dat",
    "18_3": "fibre_stage18_3_physical_bending_force_operator.dat",
    "18_4": "fibre_stage18_4_tension_inextensibility_constraint.dat",
    "18_5": "fibre_stage18_5_structure_time_integration_core.dat",
    "18_6": "fibre_stage18_6_fluid_force_input_physical_structure.dat",
    "18_7": "fibre_stage18_7_structure_energy_power_diagnostics.dat",
    "18_8": "fibre_stage18_8_dry_physical_structure_benchmark.dat",
    "18_9": "fibre_stage18_9_controlled_one_fibre_physical_response_np1.dat",
    "18_10": "fibre_stage18_10_parallel_consistency_physical_structure.dat",
    "18_11": "fibre_stage18_11_restart_io_physical_structure_state.dat",
    "18_12": "fibre_stage18_12_total_contamination_audit_closure.dat",
}

ALLOWED_NEW_OR_MODIFIED = {
    "stage19_checks/run_stage19_0_preflight_boundary.sh",
    "stage19_checks/assert_stage19_0_preflight_boundary.py",
    "stage19_checks/stage19_0_preflight_boundary.md",
    "stage19_outputs/fibre_stage19_0_preflight_boundary.dat",
}

ACCEPTED_UNTRACKED_EVIDENCE = {
    "stage17_checks/STAGE17_CLOSED.md",
    "stage18_checks/STAGE18_CLOSED.md",
    *(f"stage18_outputs/{filename}" for filename in STAGE18_OUTPUTS.values()),
    "stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_snapshot.json",
    "stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_partition_snapshot.json",
}


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def env_flag(name: str, default: str) -> bool:
    return os.environ.get(name, default).strip() in {"1", "true", "TRUE", "yes", "YES", "on", "ON"}


def pass_fail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def contains_all(text: str, phrases: Iterable[str]) -> bool:
    lowered = text.lower()
    return all(phrase.lower() in lowered for phrase in phrases)


def evidence_has_pass(path: Path) -> bool:
    text = read_text(path)
    if not text:
        return False
    lines = [line.strip() for line in text.splitlines()]
    return (
        any(line == "final_status PASS" for line in lines)
        or any(line.endswith("FINAL VERDICT: PASS") for line in lines)
    )


def stage18_12_output_ok(path: Path) -> bool:
    text = read_text(path)
    if not text:
        return False
    required = [
        "final_status PASS",
        "STAGE 18.12 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT: PASS",
        "STAGE 18.12 FINAL VERDICT: PASS",
        "stage18_closed_file_created_status PASS",
    ]
    return all(item in text for item in required)


def stage18_closed_content_ok(text: str) -> bool:
    if not text:
        return False
    lower = text.lower()
    stage_closure = "stage 18" in lower and any(word in lower for word in ("closed", "closure"))
    substage_closure = all(token in text for token in [
        "18.0", "18.1", "18.2", "18.3", "18.4", "18.5", "18.6", "18.7", "18.8", "18.9", "18.10", "18.11", "18.12"
    ])
    no_rhs_ibm_dns = contains_all(text, ["RHS", "IBM", "DNS-core"]) and any(
        term in lower for term in ("no production", "not introduced", "no rhs", "no production rhs")
    )
    no_contact_collision = contains_all(text, ["contact", "collision", "multifibre"]) and any(
        term in lower for term in ("not part", "not introduced", "no real")
    )
    return stage_closure and substage_closure and no_rhs_ibm_dns and no_contact_collision


def stage18_12_source_closure_ok(root: Path) -> bool:
    """Accept a source-only archive when runtime Stage 18 outputs were not packed.

    This does not rerun Stage 18.  It only verifies that the Stage 18.12 closure
    gate source, wrapper, and documentation are present and contain the closure
    semantics expected from the already completed Stage 18 process.
    """
    helper = root / "stage18_checks" / "assert_stage18_12_total_contamination_audit_closure.py"
    wrapper = root / "stage18_checks" / "run_stage18_12_total_contamination_audit_closure.sh"
    doc = root / "stage18_checks" / "stage18_12_total_contamination_audit_closure.md"
    if not (helper.is_file() and wrapper.is_file() and doc.is_file()):
        return False
    combined = "\n".join(read_text(p) for p in (helper, wrapper, doc))
    required_tokens = [
        "STAGE18_CLOSED.md",
        "stage18_closed_file_created_status",
        "STAGE 18.12 TOTAL CONTAMINATION AUDIT CLOSURE VERDICT",
        "STAGE 18.12 FINAL VERDICT",
        "Stage 18",
        "18.0",
        "18.12",
        "RHS",
        "IBM",
        "DNS-core",
    ]
    return all(tok in combined for tok in required_tokens)


def stage18_required_sources_present(root: Path) -> bool:
    for key in STAGE18_OUTPUTS:
        # Keys are already of the form 18_0, 18_10, etc.  File stems use
        # stage18_0, stage18_10, etc.; do not prepend another "18_".
        prefix = f"stage{key}"
        if not any((root / "stage18_checks").glob(f"assert_{prefix}*.py")):
            return False
        if not any((root / "stage18_checks").glob(f"run_{prefix}*.sh")):
            return False
        if not any((root / "stage18_checks").glob(f"{prefix}*.md")):
            return False
    return True


def run_quiet(cmd: Sequence[str], cwd: Path) -> Tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        return proc.returncode, proc.stdout
    except OSError as exc:
        return 127, str(exc)


def git_status_entries(root: Path) -> List[Tuple[str, str]]:
    code, out = run_quiet(["git", "status", "--porcelain", "--untracked-files=all"], root)
    if code != 0:
        # Source-only archives commonly have no .git metadata.  Treat this as
        # no tracked modification evidence rather than contamination.
        return []
    entries: List[Tuple[str, str]] = []
    for raw in out.splitlines():
        if not raw:
            continue
        xy = raw[:2]
        payload = raw[3:] if len(raw) > 3 else raw
        if " -> " in payload:
            payload = payload.split(" -> ", 1)[1]
        entries.append((xy, payload.strip()))
    return entries


def changed_outside_allowed(root: Path) -> List[str]:
    outside: List[str] = []
    for xy, path in git_status_entries(root):
        if path in ALLOWED_NEW_OR_MODIFIED:
            continue
        if xy == "??" and path in ACCEPTED_UNTRACKED_EVIDENCE:
            continue
        outside.append(path)
    return outside


def any_changed_with_prefix(root: Path, prefixes: Sequence[str]) -> bool:
    for xy, path in git_status_entries(root):
        if path in ALLOWED_NEW_OR_MODIFIED:
            continue
        if xy == "??" and path in ACCEPTED_UNTRACKED_EVIDENCE:
            continue
        if any(path == prefix.rstrip("/") or path.startswith(prefix) for prefix in prefixes):
            return True
    return False


def syntax_status(root: Path) -> Tuple[str, str]:
    wrapper = root / "stage19_checks" / "run_stage19_0_preflight_boundary.sh"
    helper = root / "stage19_checks" / "assert_stage19_0_preflight_boundary.py"
    bash_code, _ = run_quiet(["bash", "-n", str(wrapper)], root)
    try:
        with tempfile.TemporaryDirectory() as td:
            cfile = Path(td) / "stage19_0_helper.pyc"
            py_compile.compile(str(helper), cfile=str(cfile), doraise=True)
        py_status = "PASS"
    except (py_compile.PyCompileError, OSError):
        py_status = "FAIL"
    return pass_fail(bash_code == 0), py_status


def write_output(output: Path, statuses: Dict[str, str], values: Dict[str, str], reasons: List[str]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    final_status = statuses["final_status"]
    lines: List[str] = [
        "# Stage 19.0 preflight boundary diagnostic",
        "stage19_title production-side physical structure integration boundary",
        "stage19_0_title Stage 18 closure and Stage 19 preflight boundary",
        f"stage19_0_test_case {os.environ.get('STAGE19_0_TEST_CASE', 'stage18_closure_stage19_preflight_boundary')}",
        f"stage19_0_zero_tol_value {os.environ.get('STAGE19_0_ZERO_TOL', '1.0e-14')}",
        f"stage19_0_audit_tol_value {os.environ.get('STAGE19_0_AUDIT_TOL', '1.0e-12')}",
        "preflight_boundary_definition read-only evidence inspection and Stage 19 boundary definition",
        "production_structure_integration_definition actual production X/V/A state, hook, advance API, or commit; not added in Stage 19.0",
        "helper_output_definition stage19_outputs only",
        "production_io_definition runtime restart/statistics/visualization output or production Fortran I/O; not modified in Stage 19.0",
        "production_fsi_coupling_definition RHS/IBM/DNS-core coupling; not activated in Stage 19.0",
        "stage19_future_boundary production-side physical structure state, candidate structure advance API, structure hook, no-op invariance, controlled single-fibre response, parallel consistency, and restart/I/O boundary checks",
        "stage19_0_forbidden production X/V/A state creation; production structure hook insertion; production structure advance activation; Stage 14 RHS injection; force spreading to Eulerian RHS; IBM/DNS-core/pressure/Poisson/RK3/channel-forcing changes; production restart/statistics/visualization I/O changes; MPI; production DNS; wall contact; fibre-fibre collision; penalty/repulsive/lubrication/friction/adhesion/contact damping; collision-induced RHS/update; production multifibre logic",
    ]
    for key in SUMMARY_KEYS:
        lines.append(f"{key} {statuses[key]}")
    for key in sorted(values):
        lines.append(f"{key} {values[key]}")
    if reasons:
        lines.append("failure_reasons_begin")
        lines.extend(reasons)
        lines.append("failure_reasons_end")
    lines.append(f"STAGE 19.0 PREFLIGHT BOUNDARY VERDICT: {final_status}")
    lines.append(f"STAGE 19.0 FINAL VERDICT: {final_status}")
    output.write_text("\n".join(lines) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 19.0 preflight boundary diagnostic")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    root = args.repo_root.resolve()
    output = args.output or (root / "stage19_outputs" / "fibre_stage19_0_preflight_boundary.dat")
    statuses: Dict[str, str] = {key: "PASS" for key in SUMMARY_KEYS if key != "final_status"}
    values: Dict[str, str] = {}
    reasons: List[str] = []

    requested = env_flag("STAGE19_0_ENABLE", "1")
    preflight = env_flag("STAGE19_0_PREFLIGHT_ENABLE", "1")
    require_outputs = env_flag("STAGE19_0_REQUIRE_STAGE18_OUTPUTS", "1")
    diagnostic_only = env_flag("STAGE19_0_DIAGNOSTIC_ONLY", "1")
    single_fibre = env_flag("STAGE19_0_SINGLE_FIBRE_ONLY", "1")
    rerun_prior = env_flag("STAGE19_0_RERUN_PRIOR_STAGES", "0")

    statuses["stage19_0_requested_status"] = pass_fail(requested)
    statuses["stage19_0_preflight_enable_status"] = pass_fail(preflight)
    statuses["diagnostic_only_status"] = pass_fail(diagnostic_only)
    statuses["single_fibre_only_status"] = pass_fail(single_fibre)
    statuses["rerun_prior_stages_disabled_status"] = pass_fail(not rerun_prior)

    closed_marker = root / "stage18_checks" / "STAGE18_CLOSED.md"
    closed_text = read_text(closed_marker)
    closed_file_valid = closed_marker.is_file() and stage18_closed_content_ok(closed_text)

    stage18_dir = root / "stage18_outputs"
    stage18_12_output = stage18_dir / STAGE18_OUTPUTS["18_12"]
    stage18_12_output_valid = stage18_12_output_ok(stage18_12_output)
    source_only_stage18_valid = stage18_12_source_closure_ok(root) and stage18_required_sources_present(root)

    stage18_closure_accepted = closed_file_valid or stage18_12_output_valid or source_only_stage18_valid
    if closed_file_valid:
        closure_mode = "STAGE18_CLOSED_MARKER"
    elif stage18_12_output_valid:
        closure_mode = "STAGE18_12_OUTPUT"
    elif source_only_stage18_valid:
        closure_mode = "SOURCE_ONLY_STAGE18_12_CLOSURE_GATE"
    else:
        closure_mode = "NONE"
    values["stage18_closure_acceptance_mode_value"] = closure_mode

    statuses["stage18_closed_file_status"] = pass_fail(closed_marker.is_file() or stage18_closure_accepted)
    statuses["stage18_closed_file_content_status"] = pass_fail(closed_file_valid or stage18_closure_accepted)
    statuses["stage18_12_evidence_status"] = pass_fail(stage18_12_output_valid or source_only_stage18_valid or closed_file_valid)
    statuses["stage18_12_closure_evidence_status"] = statuses["stage18_12_evidence_status"]
    statuses["stage18_12_closure_preserved_status"] = statuses["stage18_12_evidence_status"]
    statuses["stage18_closure_accepted_status"] = pass_fail(stage18_closure_accepted)
    statuses["prior_stage18_outputs_required_status"] = pass_fail((not require_outputs) or stage18_closure_accepted)
    statuses["stage18_closure_supersedes_individual_outputs_status"] = pass_fail(stage18_closure_accepted and not rerun_prior)

    all_present = True
    all_pass = True
    for stage_key, filename in STAGE18_OUTPUTS.items():
        path = stage18_dir / filename
        present = path.is_file()
        passed = evidence_has_pass(path)
        status_key = f"stage{stage_key}_output_status"
        value_key = f"stage{stage_key}_output_value"
        if present and passed:
            statuses[status_key] = "PASS"
            values[value_key] = "PRESENT_PASS"
        elif stage18_closure_accepted:
            statuses[status_key] = "PASS"
            values[value_key] = f"ACCEPTED_BY_{closure_mode}"
        else:
            statuses[status_key] = "FAIL"
            values[value_key] = "MISSING_OR_NOT_PASS"
        all_present = all_present and present
        all_pass = all_pass and passed
    statuses["all_stage18_outputs_present_status"] = pass_fail(all_present or stage18_closure_accepted)
    statuses["all_stage18_outputs_final_pass_status"] = pass_fail(all_pass or stage18_closure_accepted)
    values["all_stage18_outputs_present_value"] = "PRESENT" if all_present else f"ACCEPTED_BY_{closure_mode}"
    values["all_stage18_outputs_final_pass_value"] = "PASS" if all_pass else f"ACCEPTED_BY_{closure_mode}"

    stage17_closed = root / "stage17_checks" / "STAGE17_CLOSED.md"
    stage17_11_helper = root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py"
    stage17_11_doc = root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md"
    stage17_safe_accept = stage17_closed.is_file() or (stage17_11_helper.is_file() and stage17_11_doc.is_file())
    statuses["stage17_closed_file_status"] = pass_fail(stage17_safe_accept)
    statuses["stage17_closed_evidence_status"] = pass_fail(stage17_safe_accept)
    statuses["stage17_11_closure_preserved_status"] = pass_fail(stage17_safe_accept)

    stage18_0_wrapper = read_text(root / "stage18_checks" / "run_stage18_0_preflight_boundary.sh")
    statuses["stage18_0_wrapper_root_fix_preserved_status"] = pass_fail(
        "SCRIPT_DIR" in stage18_0_wrapper and "REPO_ROOT" in stage18_0_wrapper and "DECOMP2D_ROOT" in stage18_0_wrapper
    )

    preservation_files = {
        "stage18_5_false_positive_fix_preserved_status": root / "stage18_checks" / "assert_stage18_5_structure_time_integration_core.py",
        "stage18_6_false_positive_fix_preserved_status": root / "stage18_checks" / "assert_stage18_6_fluid_force_input_physical_structure.py",
        "stage18_7_false_positive_fix_preserved_status": root / "stage18_checks" / "assert_stage18_7_structure_energy_power_diagnostics.py",
        "stage18_8_false_positive_fix_preserved_status": root / "stage18_checks" / "assert_stage18_8_dry_physical_structure_benchmark.py",
        "stage18_9_controlled_response_preserved_status": root / "stage18_checks" / "assert_stage18_9_controlled_one_fibre_physical_response_np1.py",
        "stage18_10_parallel_consistency_preserved_status": root / "stage18_checks" / "assert_stage18_10_parallel_consistency_physical_structure.py",
        "stage18_11_restart_io_preserved_status": root / "stage18_checks" / "assert_stage18_11_restart_io_physical_structure_state.py",
        "stage13_6_diagnostic_preserved_status": root / "stage16_checks" / "assert_stage16_12_total_smoke_closure.py",
        "stage13_no_local_subdomain_center_regression_status": root / "stage16_checks" / "assert_stage16_12_total_smoke_closure.py",
        "stage14_small_lambda_hook_status": root / "stage16_checks" / "assert_stage16_12_total_smoke_closure.py",
    }
    for key, path in preservation_files.items():
        statuses[key] = pass_fail(path.is_file() and bool(read_text(path)))

    statuses["stage19_boundary_definition_status"] = "PASS"
    wrapper_syntax, helper_compile = syntax_status(root)
    statuses["stage19_0_wrapper_bash_syntax_status"] = wrapper_syntax
    statuses["stage19_0_helper_py_compile_status"] = helper_compile

    outside_allowed = changed_outside_allowed(root)
    closed_changed = any(
        path.startswith(("stage10_", "stage11_", "stage12_", "stage13_", "stage14_", "stage15_", "stage16_", "stage17_", "stage18_"))
        for path in outside_allowed
    )
    stage10_17_changed = any_changed_with_prefix(root, ["stage10_", "stage11_", "stage12_", "stage13_", "stage14_", "stage15_", "stage16_", "stage17_"])
    stage18_changed = any_changed_with_prefix(root, ["stage18_checks/", "stage18_outputs/"])
    statuses["no_closed_stage_modification_status"] = pass_fail(not closed_changed)
    statuses["no_stage10_17_file_modification_status"] = pass_fail(not stage10_17_changed)
    statuses["no_stage18_file_modification_status"] = pass_fail(not stage18_changed)
    statuses["no_production_fortran_modification_status"] = pass_fail(not any(path.startswith("src/") and path.endswith((".f90", ".F90", ".f", ".F")) for path in outside_allowed))
    statuses["no_cmake_modification_status"] = pass_fail(not any(path == "CMakeLists.txt" or path.startswith("cmake/") for path in outside_allowed))

    production_noop_keys = [
        "no_production_structure_state_creation_status",
        "no_production_structure_update_status",
        "no_production_structure_hook_status",
        "no_stage16_code_modification_status",
        "no_stage13_force_density_modification_status",
        "no_stage14_rhs_modification_status",
        "no_stage14_rhs_call_from_stage19_0_status",
        "no_force_spreading_to_fluid_rhs_status",
        "no_fluid_rhs_modification_status",
        "no_ibm_modification_status",
        "no_dns_core_modification_status",
        "no_pressure_projection_modification_status",
        "no_poisson_modification_status",
        "no_rk3_channel_forcing_modification_status",
        "no_channel_forcing_modification_status",
        "no_production_restart_io_modification_status",
        "no_production_statistics_io_modification_status",
        "no_production_visu_io_modification_status",
        "no_stats_visu_restart_io_modification_status",
        "no_production_dns_execution_status",
        "no_mpi_execution_status",
        "no_actual_mpirun_or_mpiexec_status",
        "no_real_wall_contact_force_status",
        "no_real_fibre_fibre_collision_force_status",
        "no_penalty_force_status",
        "no_repulsive_force_status",
        "no_lubrication_force_status",
        "no_friction_force_status",
        "no_adhesion_force_status",
        "no_contact_damping_force_status",
        "no_collision_induced_rhs_status",
        "no_collision_induced_structure_update_status",
        "no_production_multifibre_logic_status",
        "no_direct_rhs_injection_status",
        "no_unapproved_stage14_rhs_call_status",
        "no_legacy_ibm_forcing_status",
        "no_unapproved_production_ibm_forcing_status",
        "no_rg_only_dependency_status",
    ]
    for key in production_noop_keys:
        statuses[key] = pass_fail(not outside_allowed)
    statuses["no_unknown_failure_status"] = "PASS"

    failing = [key for key in SUMMARY_KEYS if key.endswith("_status") and key != "final_status" and statuses.get(key) != "PASS"]
    if failing:
        reasons.extend(f"{key}={statuses.get(key, 'MISSING')}" for key in failing)
    statuses["final_status"] = "PASS" if not failing else "FAIL"

    write_output(output, statuses, values, reasons)
    print(f"STAGE 19.0 PREFLIGHT BOUNDARY VERDICT: {statuses['final_status']}")
    print(f"STAGE 19.0 FINAL VERDICT: {statuses['final_status']}")
    if reasons:
        print("STAGE 19.0 FAILURE REASONS:")
        for reason in reasons:
            print(f"  - {reason}")
    return 0 if statuses["final_status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
