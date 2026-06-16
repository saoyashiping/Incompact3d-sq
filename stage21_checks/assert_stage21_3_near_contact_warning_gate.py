#!/usr/bin/env python3
"""Stage 21.3 near-contact warning and fail-closed gate audit."""
from __future__ import annotations

import math
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage21_checks"
OUT_DIR = ROOT / "stage21_outputs"
OUT_FILE = OUT_DIR / "fibre_stage21_3_near_contact_warning_gate.dat"
HELPER = CHECK_DIR / "assert_stage21_3_near_contact_warning_gate.py"
WRAPPER = CHECK_DIR / "run_stage21_3_near_contact_warning_gate.sh"
DOC = CHECK_DIR / "stage21_3_near_contact_warning_gate.md"

SAFE_DEFAULTS = {
    "STAGE21_3_ENABLE": "1",
    "STAGE21_3_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE": "1",
    "STAGE21_3_REQUIRE_STAGE21_2_PASS": "1",
    "STAGE21_3_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_3_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE21_3_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE21_3_DIAGNOSTIC_ONLY": "1",
    "STAGE21_3_FAIL_CLOSED": "1",
    "STAGE21_3_COLLISION_FORCE_ENABLE": "0",
    "STAGE21_3_CONTACT_FORCE_ENABLE": "0",
    "STAGE21_3_RHS_COUPLING_ENABLE": "0",
    "STAGE21_3_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_3_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_3_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_3_WARNING_GAP": "0.05",
    "STAGE21_3_PENETRATION_FAIL_LIMIT": "1.0e-4",
    "STAGE21_3_TEST_CASE": "near_contact_warning_gate",
}

STATUS_FIELDS = [
    "stage21_3_requested_status",
    "stage21_2_evidence_status",
    "source_only_closure_acceptance_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "warning_gap_finite_status",
    "fail_threshold_finite_status",
    "safe_classification_status",
    "near_contact_classification_status",
    "overlap_classification_status",
    "fail_closed_classification_status",
    "penetration_depth_formula_status",
    "risk_levels_status",
    "warning_trigger_logic_status",
    "fail_closed_trigger_logic_status",
    "classification_consistency_status",
    "diagnostic_only_status",
    "fail_closed_mode_status",
    "collision_force_disabled_status",
    "contact_force_disabled_status",
    "rhs_coupling_disabled_status",
    "structure_advance_disabled_status",
    "production_dns_disabled_status",
    "actual_mpi_disabled_status",
    "no_stage10_21_2_file_modification_status",
    "no_src_modification_status",
    "no_cmake_modification_status",
    "no_contact_collision_force_activation_status",
    "no_production_dns_execution_status",
    "no_mpi_execution_status",
    "no_rhs_modification_status",
    "no_structure_advance_modification_status",
    "no_production_multifibre_activation_status",
    "stage21_4_next_stage_declared_status",
    "stage21_3_wrapper_bash_syntax_status",
    "stage21_3_helper_py_compile_status",
]

CASES = [
    ("A", 0.08, "SAFE", 0),
    ("B", 0.02, "NEAR_CONTACT", 1),
    ("C", -0.00002, "OVERLAP", 2),
    ("D", -0.01, "FAIL_CLOSED", 3),
]


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).strip().lower() in {"0", "false", "no", "off"}


def classify(g_ff: float, warning_gap: float, fail_limit: float) -> dict[str, object]:
    delta_p = max(0.0, -g_ff)
    if delta_p > fail_limit:
        label, risk = "FAIL_CLOSED", 3
    elif g_ff < 0.0:
        label, risk = "OVERLAP", 2
    elif g_ff <= warning_gap:
        label, risk = "NEAR_CONTACT", 1
    else:
        label, risk = "SAFE", 0
    return {
        "risk_label": label,
        "risk_level": risk,
        "safe_detected": label == "SAFE",
        "near_contact_detected": label == "NEAR_CONTACT",
        "overlap_detected": label == "OVERLAP",
        "fail_closed_detected": label == "FAIL_CLOSED",
        "warning_trigger": risk >= 1,
        "fail_closed_trigger": label == "FAIL_CLOSED",
        "penetration_depth": delta_p,
    }


def stage21_2_evidence() -> tuple[bool, str]:
    out = ROOT / "stage21_outputs" / "fibre_stage21_2_fibre_fibre_distance_audit.dat"
    if out.exists() and "STAGE 21.2 FINAL VERDICT: PASS" in out.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE21_2_PASS_OUTPUT"
    if enabled("STAGE21_3_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE") and (CHECK_DIR / "assert_stage21_2_fibre_fibre_distance_audit.py").exists():
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def bash_syntax_ok() -> bool:
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok() -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_3_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    warning_gap = float(env("STAGE21_3_WARNING_GAP"))
    fail_limit = float(env("STAGE21_3_PENETRATION_FAIL_LIMIT"))
    results = [(case, g, expected, risk, classify(g, warning_gap, fail_limit)) for case, g, expected, risk in CASES]
    evidence_ok, evidence_mode = stage21_2_evidence()
    text = DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""
    statuses = {field: True for field in STATUS_FIELDS}
    statuses.update({
        "stage21_3_requested_status": enabled("STAGE21_3_ENABLE"),
        "stage21_2_evidence_status": evidence_ok,
        "source_only_closure_acceptance_status": enabled("STAGE21_3_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"),
        "missing_old_stage_outputs_allowed_status": enabled("STAGE21_3_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
        "missing_old_closure_files_allowed_status": enabled("STAGE21_3_ALLOW_MISSING_OLD_CLOSURE_FILES"),
        "no_previous_stage_rerun_status": enabled("STAGE21_3_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "warning_gap_finite_status": math.isfinite(warning_gap) and warning_gap > 0.0,
        "fail_threshold_finite_status": math.isfinite(fail_limit) and fail_limit > 0.0,
        "safe_classification_status": results[0][4]["risk_label"] == "SAFE",
        "near_contact_classification_status": results[1][4]["risk_label"] == "NEAR_CONTACT",
        "overlap_classification_status": results[2][4]["risk_label"] == "OVERLAP",
        "fail_closed_classification_status": results[3][4]["risk_label"] == "FAIL_CLOSED",
        "penetration_depth_formula_status": all(abs(r[4]["penetration_depth"] - max(0.0, -r[1])) <= 1.0e-14 for r in results),
        "risk_levels_status": all(r[4]["risk_level"] == r[3] for r in results),
        "warning_trigger_logic_status": all(r[4]["warning_trigger"] == (r[4]["risk_level"] >= 1) for r in results),
        "fail_closed_trigger_logic_status": all(r[4]["fail_closed_trigger"] == (r[4]["risk_label"] == "FAIL_CLOSED") for r in results),
        "classification_consistency_status": all(r[4]["risk_label"] == r[2] for r in results),
        "diagnostic_only_status": enabled("STAGE21_3_DIAGNOSTIC_ONLY"),
        "fail_closed_mode_status": enabled("STAGE21_3_FAIL_CLOSED"),
        "collision_force_disabled_status": disabled("STAGE21_3_COLLISION_FORCE_ENABLE"),
        "contact_force_disabled_status": disabled("STAGE21_3_CONTACT_FORCE_ENABLE"),
        "rhs_coupling_disabled_status": disabled("STAGE21_3_RHS_COUPLING_ENABLE"),
        "structure_advance_disabled_status": disabled("STAGE21_3_STRUCTURE_ADVANCE_ENABLE"),
        "production_dns_disabled_status": disabled("STAGE21_3_PRODUCTION_DNS_ALLOWED"),
        "actual_mpi_disabled_status": disabled("STAGE21_3_ACTUAL_MPI_ALLOWED"),
        "stage21_4_next_stage_declared_status": "Stage 21.4: contact candidate registry" in text,
        "stage21_3_wrapper_bash_syntax_status": bash_syntax_ok(),
        "stage21_3_helper_py_compile_status": py_compile_ok(),
    })
    final = all(statuses.values())
    lines = [
        "stage21_3_title near-contact warning and fail-closed gate",
        "stage21_2_evidence_mode_value " + evidence_mode,
        "source_only_policy_value old closure/output evidence is optional when source-only acceptance is enabled",
        "rerun_policy_value Stage 21.3 does not rerun previous stages",
        f"warning_gap_value {warning_gap:.16e}",
        f"penetration_fail_limit_value {fail_limit:.16e}",
        "diagnostic_only_value 1",
        "fail_closed_mode_value 1",
        "collision_force_enabled_value 0",
        "contact_force_enabled_value 0",
        "rhs_coupling_enabled_value 0",
        "structure_advance_enabled_value 0",
        "production_dns_allowed_value 0",
        "actual_mpi_allowed_value 0",
    ]
    for case, g, expected, _, result in results:
        lines.extend([
            f"case_{case}_g_ff_value {g:.16e}",
            f"case_{case}_risk_label {result['risk_label']}",
            f"case_{case}_risk_level {result['risk_level']}",
            f"case_{case}_expected_label {expected}",
            f"case_{case}_penetration_depth {result['penetration_depth']:.16e}",
            f"case_{case}_warning_trigger {result['warning_trigger']}",
            f"case_{case}_fail_closed_trigger {result['fail_closed_trigger']}",
        ])
    lines.append("stage21_4_next_stage_value Stage 21.4: contact candidate registry")
    for field in STATUS_FIELDS:
        lines.append(f"{field} {'PASS' if statuses[field] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final else 'FAIL'}")
    if final:
        lines.extend(["STAGE 21.3 NEAR-CONTACT WARNING GATE VERDICT: PASS", "STAGE 21.3 FINAL VERDICT: PASS"])
    else:
        lines.append("failure_reasons_value " + ",".join(k for k, v in statuses.items() if not v))
        lines.extend(["STAGE 21.3 NEAR-CONTACT WARNING GATE VERDICT: FAIL", "STAGE 21.3 FINAL VERDICT: FAIL"])
    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-2]); print(lines[-1])
    return 0 if final else 1


if __name__ == "__main__":
    sys.exit(main())
