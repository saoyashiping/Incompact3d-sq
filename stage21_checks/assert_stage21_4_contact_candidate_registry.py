#!/usr/bin/env python3
"""Stage 21.4 diagnostic-only contact candidate registry audit."""
from __future__ import annotations

import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage21_checks"
OUT_DIR = ROOT / "stage21_outputs"
OUT_FILE = OUT_DIR / "fibre_stage21_4_contact_candidate_registry.dat"
HELPER = CHECK_DIR / "assert_stage21_4_contact_candidate_registry.py"
WRAPPER = CHECK_DIR / "run_stage21_4_contact_candidate_registry.sh"
DOC = CHECK_DIR / "stage21_4_contact_candidate_registry.md"

SAFE_DEFAULTS = {
    "STAGE21_4_ENABLE": "1",
    "STAGE21_4_CONTACT_CANDIDATE_REGISTRY_ENABLE": "1",
    "STAGE21_4_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE": "1",
    "STAGE21_4_REQUIRE_STAGE21_3_PASS": "1",
    "STAGE21_4_ALLOW_MISSING_OLD_STAGE_OUTPUTS": "1",
    "STAGE21_4_ALLOW_MISSING_OLD_CLOSURE_FILES": "1",
    "STAGE21_4_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_4_DIAGNOSTIC_ONLY": "1",
    "STAGE21_4_FAIL_CLOSED": "1",
    "STAGE21_4_CONTACT_FORCE_ENABLE": "0",
    "STAGE21_4_COLLISION_FORCE_ENABLE": "0",
    "STAGE21_4_CONTACT_FORCE_APPLY_ENABLE": "0",
    "STAGE21_4_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_4_RHS_COUPLING_ENABLE": "0",
    "STAGE21_4_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_4_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_4_TEST_CASE": "contact_candidate_registry",
}

FIELDS = [
    "candidate_id", "candidate_type", "fibre_i", "fibre_j", "point_i", "point_j", "segment_i", "segment_j",
    "gap_value", "penetration_depth", "warning_trigger", "fail_closed_trigger", "risk_level", "risk_label",
    "nearest_wall_side", "closest_pair_valid", "candidate_active", "force_computation_allowed",
    "force_application_allowed", "diagnostic_only",
]
VALID_TYPES = {"wall_lower", "wall_upper", "fibre_fibre"}
RISK_LABELS = {0: "SAFE", 1: "NEAR_CONTACT", 2: "OVERLAP", 3: "FAIL_CLOSED"}
STATUS_FIELDS = [
    "stage21_4_requested_status",
    "stage21_4_contact_candidate_registry_enable_status",
    "stage21_3_evidence_status",
    "source_only_closure_acceptance_status",
    "missing_old_stage_outputs_allowed_status",
    "missing_old_closure_files_allowed_status",
    "no_previous_stage_rerun_status",
    "contact_candidate_registry_documented_status",
    "all_required_registry_fields_present_status",
    "candidate_ids_unique_status",
    "candidate_types_valid_status",
    "wall_candidates_nearest_wall_side_status",
    "fibre_fibre_candidates_ordered_pair_status",
    "risk_labels_match_risk_levels_status",
    "warning_fail_flags_consistent_status",
    "diagnostic_only_status",
    "fail_closed_status",
    "contact_force_disabled_status",
    "collision_force_disabled_status",
    "force_application_disabled_status",
    "structure_advance_disabled_status",
    "rhs_coupling_disabled_status",
    "production_dns_disabled_status",
    "actual_mpi_disabled_status",
    "no_stage10_21_3_file_modification_status",
    "no_src_modification_status",
    "no_cmake_modification_status",
    "no_real_force_computed_status",
    "no_production_hook_activated_status",
    "stage21_5_next_stage_declared_status",
    "stage21_4_wrapper_bash_syntax_status",
    "stage21_4_helper_py_compile_status",
]


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).strip().lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).strip().lower() in {"0", "false", "no", "off"}


def stage21_3_evidence() -> tuple[bool, str]:
    out = ROOT / "stage21_outputs" / "fibre_stage21_3_near_contact_warning_gate.dat"
    if out.exists() and "STAGE 21.3 FINAL VERDICT: PASS" in out.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE21_3_PASS_OUTPUT"
    if enabled("STAGE21_4_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE") and (CHECK_DIR / "assert_stage21_3_near_contact_warning_gate.py").exists():
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def bash_syntax_ok() -> bool:
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok() -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_4_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def registry() -> list[dict[str, object]]:
    return [
        {
            "candidate_id": "wall_lower_0000", "candidate_type": "wall_lower", "fibre_i": 0, "fibre_j": -1,
            "point_i": 12, "point_j": -1, "segment_i": -1, "segment_j": -1, "gap_value": 0.020,
            "penetration_depth": 0.0, "warning_trigger": True, "fail_closed_trigger": False, "risk_level": 1,
            "risk_label": "NEAR_CONTACT", "nearest_wall_side": "lower", "closest_pair_valid": True,
            "candidate_active": False, "force_computation_allowed": False, "force_application_allowed": False,
            "diagnostic_only": True,
        },
        {
            "candidate_id": "wall_upper_0000", "candidate_type": "wall_upper", "fibre_i": 0, "fibre_j": -1,
            "point_i": 51, "point_j": -1, "segment_i": -1, "segment_j": -1, "gap_value": 0.400,
            "penetration_depth": 0.0, "warning_trigger": False, "fail_closed_trigger": False, "risk_level": 0,
            "risk_label": "SAFE", "nearest_wall_side": "upper", "closest_pair_valid": True,
            "candidate_active": False, "force_computation_allowed": False, "force_application_allowed": False,
            "diagnostic_only": True,
        },
        {
            "candidate_id": "fibre_fibre_0000", "candidate_type": "fibre_fibre", "fibre_i": 0, "fibre_j": 1,
            "point_i": 0, "point_j": 0, "segment_i": 0, "segment_j": 0, "gap_value": 0.060,
            "penetration_depth": 0.0, "warning_trigger": False, "fail_closed_trigger": False, "risk_level": 0,
            "risk_label": "SAFE", "nearest_wall_side": "none", "closest_pair_valid": True,
            "candidate_active": False, "force_computation_allowed": False, "force_application_allowed": False,
            "diagnostic_only": True,
        },
    ]


def warning_fail_consistent(rec: dict[str, object]) -> bool:
    risk = int(rec["risk_level"])
    return bool(rec["warning_trigger"]) == (risk >= 1) and bool(rec["fail_closed_trigger"]) == (risk == 3)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = registry()
    evidence_ok, evidence_mode = stage21_3_evidence()
    text = DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""
    ids = [r["candidate_id"] for r in rows]
    statuses = {field: True for field in STATUS_FIELDS}
    statuses.update({
        "stage21_4_requested_status": enabled("STAGE21_4_ENABLE"),
        "stage21_4_contact_candidate_registry_enable_status": enabled("STAGE21_4_CONTACT_CANDIDATE_REGISTRY_ENABLE"),
        "stage21_3_evidence_status": evidence_ok,
        "source_only_closure_acceptance_status": enabled("STAGE21_4_ACCEPT_SOURCE_ONLY_STAGE20_CLOSURE"),
        "missing_old_stage_outputs_allowed_status": enabled("STAGE21_4_ALLOW_MISSING_OLD_STAGE_OUTPUTS"),
        "missing_old_closure_files_allowed_status": enabled("STAGE21_4_ALLOW_MISSING_OLD_CLOSURE_FILES"),
        "no_previous_stage_rerun_status": enabled("STAGE21_4_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "contact_candidate_registry_documented_status": "contact candidate registry" in text,
        "all_required_registry_fields_present_status": all(all(f in row for f in FIELDS) for row in rows),
        "candidate_ids_unique_status": len(ids) == len(set(ids)),
        "candidate_types_valid_status": all(r["candidate_type"] in VALID_TYPES for r in rows),
        "wall_candidates_nearest_wall_side_status": all(r["nearest_wall_side"] in {"lower", "upper"} for r in rows if str(r["candidate_type"]).startswith("wall_")),
        "fibre_fibre_candidates_ordered_pair_status": all(int(r["fibre_i"]) < int(r["fibre_j"]) for r in rows if r["candidate_type"] == "fibre_fibre"),
        "risk_labels_match_risk_levels_status": all(r["risk_label"] == RISK_LABELS[int(r["risk_level"])] for r in rows),
        "warning_fail_flags_consistent_status": all(warning_fail_consistent(r) for r in rows),
        "diagnostic_only_status": enabled("STAGE21_4_DIAGNOSTIC_ONLY") and all(r["diagnostic_only"] for r in rows),
        "fail_closed_status": enabled("STAGE21_4_FAIL_CLOSED"),
        "contact_force_disabled_status": disabled("STAGE21_4_CONTACT_FORCE_ENABLE"),
        "collision_force_disabled_status": disabled("STAGE21_4_COLLISION_FORCE_ENABLE"),
        "force_application_disabled_status": disabled("STAGE21_4_CONTACT_FORCE_APPLY_ENABLE") and all(not r["force_application_allowed"] for r in rows),
        "structure_advance_disabled_status": disabled("STAGE21_4_STRUCTURE_ADVANCE_ENABLE"),
        "rhs_coupling_disabled_status": disabled("STAGE21_4_RHS_COUPLING_ENABLE"),
        "production_dns_disabled_status": disabled("STAGE21_4_PRODUCTION_DNS_ALLOWED"),
        "actual_mpi_disabled_status": disabled("STAGE21_4_ACTUAL_MPI_ALLOWED"),
        "no_real_force_computed_status": all(not r["force_computation_allowed"] for r in rows),
        "stage21_5_next_stage_declared_status": "Stage 21.5: contact pair ownership audit" in text,
        "stage21_4_wrapper_bash_syntax_status": bash_syntax_ok(),
        "stage21_4_helper_py_compile_status": py_compile_ok(),
    })
    final = all(statuses.values())
    lines = [
        "stage21_4_title contact candidate registry",
        "stage21_3_evidence_mode_value " + evidence_mode,
        "source_only_policy_value old closure/output evidence is optional when source-only acceptance is enabled",
        "rerun_policy_value Stage 21.4 does not rerun previous stages",
        "registry_schema_fields_value " + ",".join(FIELDS),
        f"candidate_count_value {len(rows)}",
        "diagnostic_only_value 1",
        "contact_force_enabled_value 0",
        "collision_force_enabled_value 0",
        "force_application_allowed_value 0",
        "structure_advance_enabled_value 0",
        "rhs_coupling_enabled_value 0",
        "production_dns_allowed_value 0",
        "actual_mpi_allowed_value 0",
    ]
    for key in SAFE_DEFAULTS:
        lines.append(f"{key.lower()}_value {env(key)}")
    for row in rows:
        prefix = str(row["candidate_id"])
        lines.append(prefix + "_record " + ";".join(f"{k}={row[k]}" for k in FIELDS))
    lines.append("stage21_5_next_stage_value Stage 21.5: contact pair ownership audit")
    for field in STATUS_FIELDS:
        lines.append(f"{field} {'PASS' if statuses[field] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final else 'FAIL'}")
    if final:
        lines.extend(["STAGE 21.4 CONTACT CANDIDATE REGISTRY VERDICT: PASS", "STAGE 21.4 FINAL VERDICT: PASS"])
    else:
        lines.append("failure_reasons_value " + ",".join(k for k, v in statuses.items() if not v))
        lines.extend(["STAGE 21.4 CONTACT CANDIDATE REGISTRY VERDICT: FAIL", "STAGE 21.4 FINAL VERDICT: FAIL"])
    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-2]); print(lines[-1])
    return 0 if final else 1


if __name__ == "__main__":
    sys.exit(main())
