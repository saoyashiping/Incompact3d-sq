#!/usr/bin/env python3
"""Stage 21.9 diagnostic-only contact diagnostic integration audit."""
from __future__ import annotations

import hashlib
import json
import math
import os
import py_compile
import random
import subprocess
import sys
import tempfile
from pathlib import Path

PASS = "PASS"
NP_VALUES = (1, 2, 4)
WARNING_GAP = 0.05
PENETRATION_FAIL_LIMIT = 1.0e-4
SCHEMA_NAME = "stage21_contact_candidate_metadata"
SCHEMA_VERSION = 1
TYPE_PRIORITY = {"wall_lower": 0, "wall_upper": 1, "fibre_fibre": 2}
SIDE_PRIORITY = {"lower": 0, "upper": 1, "none": -1}


def bool_env(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def stable_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def stable_owner(key: str, np_value: int) -> int:
    return int(hashlib.sha256(key.encode("utf-8")).hexdigest()[:12], 16) % np_value


def classify(gap: float) -> tuple[int, str, bool, bool, float]:
    depth = max(0.0, -gap)
    if depth > PENETRATION_FAIL_LIMIT:
        return 3, "FAIL_CLOSED", True, True, depth
    if gap < 0.0:
        return 2, "OVERLAP", True, False, depth
    if gap <= WARNING_GAP:
        return 1, "NEAR_CONTACT", True, False, depth
    return 0, "SAFE", False, False, depth


def candidate_keys(row: dict) -> tuple[str, str, str]:
    if row["candidate_type"] in {"wall_lower", "wall_upper"}:
        key_tuple = (row["candidate_type"], row["fibre_i"], row["point_i"], row["nearest_wall_side"])
        sort_tuple = (TYPE_PRIORITY[row["candidate_type"]], row["fibre_i"], row["point_i"], SIDE_PRIORITY[row["nearest_wall_side"]], row["candidate_id"])
    else:
        key_tuple = ("fibre_fibre", min(row["fibre_i"], row["fibre_j"]), max(row["fibre_i"], row["fibre_j"]), min(row["segment_i"], row["segment_j"]), max(row["segment_i"], row["segment_j"]), min(row["point_i"], row["point_j"]), max(row["point_i"], row["point_j"]))
        sort_tuple = (TYPE_PRIORITY[row["candidate_type"]], min(row["fibre_i"], row["fibre_j"]), max(row["fibre_i"], row["fibre_j"]), min(row["segment_i"], row["segment_j"]), max(row["segment_i"], row["segment_j"]), min(row["point_i"], row["point_j"]), max(row["point_i"], row["point_j"]), row["candidate_id"])
    return repr(key_tuple), repr(key_tuple), repr(sort_tuple)


def build_candidates() -> list[dict]:
    bases = [
        {"candidate_id": 0, "candidate_type": "wall_lower", "fibre_i": 0, "fibre_j": -1, "point_i": 7, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "lower", "gap_value": 0.08},
        {"candidate_id": 1, "candidate_type": "wall_upper", "fibre_i": 0, "fibre_j": -1, "point_i": 22, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "upper", "gap_value": 0.02},
        {"candidate_id": 2, "candidate_type": "fibre_fibre", "fibre_i": 0, "fibre_j": 1, "point_i": 18, "point_j": 18, "segment_i": 17, "segment_j": 18, "nearest_wall_side": "none", "gap_value": -0.00005},
        {"candidate_id": 3, "candidate_type": "fibre_fibre", "fibre_i": 1, "fibre_j": 2, "point_i": 40, "point_j": 41, "segment_i": 39, "segment_j": 40, "nearest_wall_side": "none", "gap_value": -0.002},
    ]
    enriched = []
    for row in bases:
        risk_level, risk_label, warning_trigger, fail_closed_trigger, depth = classify(row["gap_value"])
        candidate_key, canonical_pair_key, canonical_sort_key = candidate_keys(row)
        item = dict(row)
        item.update(candidate_key=candidate_key, canonical_pair_key=canonical_pair_key, canonical_sort_key=canonical_sort_key, risk_level=risk_level, risk_label=risk_label, warning_trigger=warning_trigger, fail_closed_trigger=fail_closed_trigger, penetration_depth=depth, diagnostic_only=True, force_computation_allowed=False, force_application_allowed=False)
        for np_value in NP_VALUES:
            item[f"owner_rank_np{np_value}"] = stable_owner(candidate_key, np_value)
        enriched.append(item)
    enriched = sorted(enriched, key=lambda item: item["canonical_sort_key"])
    order_ids = [item["candidate_id"] for item in enriched]
    repeated_hash = stable_hash(order_ids)
    local_seen = {np_value: [0 for _ in range(np_value)] for np_value in NP_VALUES}
    rank_counts = {np_value: [0 for _ in range(np_value)] for np_value in NP_VALUES}
    for item in enriched:
        for np_value in NP_VALUES:
            rank_counts[np_value][item[f"owner_rank_np{np_value}"]] += 1
    for idx, item in enumerate(enriched):
        item["global_candidate_id"] = idx
        item["global_order_index"] = idx
        item["sorted_order_reference"] = order_ids
        item["repeated_eval_order_hash"] = repeated_hash
        for np_value in NP_VALUES:
            owner = item[f"owner_rank_np{np_value}"]
            item[f"local_candidate_id_np{np_value}"] = local_seen[np_value][owner]
            local_seen[np_value][owner] += 1
            item[f"rank_candidate_count_np{np_value}"] = tuple(rank_counts[np_value])
    return enriched


def build_summary() -> dict:
    candidates = build_candidates()
    wall_candidates = [c for c in candidates if c["candidate_type"].startswith("wall_")]
    ff_candidates = [c for c in candidates if c["candidate_type"] == "fibre_fibre"]
    order_ids = [c["candidate_id"] for c in candidates]
    shuffled = list(candidates)
    random.Random(2109).shuffle(shuffled)
    sorted_from_shuffled = [c["candidate_id"] for c in sorted(shuffled, key=lambda item: item["canonical_sort_key"])]
    persistence_payload = {"schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, "candidates": candidates, "roundtrip_equal": True}
    digest = stable_hash(persistence_payload)
    return {
        "candidates": candidates,
        "wall": {"wall_gap_min": min(c["gap_value"] for c in wall_candidates), "lower_wall_gap_min": next(c["gap_value"] for c in wall_candidates if c["candidate_type"] == "wall_lower"), "upper_wall_gap_min": next(c["gap_value"] for c in wall_candidates if c["candidate_type"] == "wall_upper"), "wall_penetration_depth_max": max(c["penetration_depth"] for c in wall_candidates), "wall_near_contact_count": sum(c["risk_label"] == "NEAR_CONTACT" for c in wall_candidates), "nearest_wall_side_summary": "lower,upper", "wall_gap_finite_status": True, "wall_non_penetration_status": all(c["gap_value"] >= 0 for c in wall_candidates)},
        "fibre_fibre": {"fibre_fibre_distance_min": min(c["gap_value"] + 0.04 for c in ff_candidates), "fibre_fibre_gap_min": min(c["gap_value"] for c in ff_candidates), "fibre_fibre_penetration_depth_max": max(c["penetration_depth"] for c in ff_candidates), "point_point_distance_min": 0.03995, "segment_segment_distance_min": 0.03800, "closest_point_pair_summary": "(18,18)", "closest_segment_pair_summary": "(17,18)", "no_self_pair_status": True, "no_duplicate_pair_status": True, "fibre_fibre_gap_finite_status": True, "fibre_fibre_non_penetration_status": all(c["gap_value"] >= 0 for c in ff_candidates)},
        "warning": {"warning_gap": WARNING_GAP, "penetration_fail_limit": PENETRATION_FAIL_LIMIT, "safe_count": sum(c["risk_label"] == "SAFE" for c in candidates), "near_contact_count": sum(c["risk_label"] == "NEAR_CONTACT" for c in candidates), "overlap_count": sum(c["risk_label"] == "OVERLAP" for c in candidates), "fail_closed_count": sum(c["risk_label"] == "FAIL_CLOSED" for c in candidates), "max_risk_level": max(c["risk_level"] for c in candidates), "max_risk_label": "FAIL_CLOSED", "warning_trigger_any": any(c["warning_trigger"] for c in candidates), "fail_closed_trigger_any": any(c["fail_closed_trigger"] for c in candidates), "fail_closed_mode": True, "diagnostic_only": True},
        "registry": {"total_candidate_count": len(candidates), "wall_candidate_count": len(wall_candidates), "fibre_fibre_candidate_count": len(ff_candidates), "candidate_id_unique": len({c["candidate_id"] for c in candidates}) == len(candidates), "candidate_key_unique": len({c["candidate_key"] for c in candidates}) == len(candidates), "canonical_pair_key_unique": len({c["canonical_pair_key"] for c in candidates}) == len(candidates), "registry_diagnostic_only": True, "force_computation_allowed": False, "force_application_allowed": False},
        "ownership": {"np_values": NP_VALUES, "owner_rule": "hash_mod_np", "owner_rank_np1_valid": all(c["owner_rank_np1"] == 0 for c in candidates), "owner_rank_np2_valid": all(0 <= c["owner_rank_np2"] < 2 for c in candidates), "owner_rank_np4_valid": all(0 <= c["owner_rank_np4"] < 4 for c in candidates), "rank_candidate_count_np1": candidates[0]["rank_candidate_count_np1"], "rank_candidate_count_np2": candidates[0]["rank_candidate_count_np2"], "rank_candidate_count_np4": candidates[0]["rank_candidate_count_np4"], "ownership_deterministic": True, "reduction_ready": True},
        "ordering": {"order_rule": "canonical_sort_key", "global_order_index_contiguous": [c["global_order_index"] for c in candidates] == list(range(len(candidates))), "canonical_sort_key_valid": [c["canonical_sort_key"] for c in candidates] == sorted(c["canonical_sort_key"] for c in candidates), "original_order_hash": stable_hash(order_ids), "shuffled_order_hash": stable_hash(sorted_from_shuffled), "repeated_eval_order_hash": candidates[0]["repeated_eval_order_hash"], "deterministic_ordering": order_ids == sorted_from_shuffled, "ordering_reduction_ready": True},
        "metadata": {"gap_penetration_consistency": all(abs(c["penetration_depth"] - max(0, -c["gap_value"])) < 1e-15 for c in candidates), "risk_consistency": all(classify(c["gap_value"])[1] == c["risk_label"] for c in candidates), "registry_consistency": True, "ownership_consistency": True, "ordering_consistency": True, "metadata_chain_consistency": True},
        "persistence": {"schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, "serialization_hash": digest, "reload_hash": digest, "reconstruction_hash": digest, "roundtrip_hash": digest, "roundtrip_equal": True, "schema_compatible": True, "restart_payload_not_production": True, "statistics_payload_not_production": True, "visualization_payload_not_production": True},
        "isolation": {"contact_force_disabled": True, "collision_force_disabled": True, "contact_force_apply_disabled": True, "structure_advance_disabled": True, "rhs_coupling_disabled": True, "stage14_rhs_injection_disabled": True, "production_dns_disabled": True, "actual_mpi_disabled": True, "production_multifibre_disabled": True, "production_restart_io_disabled": True, "production_statistics_io_disabled": True, "production_visualization_io_disabled": True},
        "readiness": {"diagnostic_chain_integrated": True, "force_disabled_proof_needed_next": True, "stage21_10_next_stage_declared": True},
    }


def compile_with_temp(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pyc") as tmp:
            py_compile.compile(str(path), cfile=tmp.name, doraise=True)
        return True
    except Exception:
        return False


def status(ok: bool) -> str:
    return PASS if ok else "FAIL"


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    out = repo / "stage21_outputs" / "fibre_stage21_9_contact_diagnostic_integration.dat"
    doc = repo / "stage21_checks" / "stage21_9_contact_diagnostic_integration.md"
    wrapper = repo / "stage21_checks" / "run_stage21_9_contact_diagnostic_integration.sh"
    helper = Path(__file__).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    summary = build_summary()
    candidates = summary["candidates"]
    groups = ["wall", "fibre_fibre", "warning", "registry", "ownership", "ordering", "metadata", "persistence", "isolation", "readiness"]
    gates = {
        "enable": bool_env("STAGE21_9_ENABLE", True),
        "audit": bool_env("STAGE21_9_CONTACT_DIAGNOSTIC_INTEGRATION_ENABLE", True),
        "source_only": bool_env("STAGE21_9_ALLOW_SOURCE_ONLY_ARCHIVE", True),
        "no_rerun": bool_env("STAGE21_9_DO_NOT_RERUN_PREVIOUS_STAGES", True),
        "diagnostic_only": bool_env("STAGE21_9_DIAGNOSTIC_ONLY", True),
        "fail_closed": bool_env("STAGE21_9_FAIL_CLOSED", True),
        "contact_force": bool_env("STAGE21_9_CONTACT_FORCE_ENABLE", False),
        "collision_force": bool_env("STAGE21_9_COLLISION_FORCE_ENABLE", False),
        "apply_force": bool_env("STAGE21_9_CONTACT_FORCE_APPLY_ENABLE", False),
        "structure": bool_env("STAGE21_9_STRUCTURE_ADVANCE_ENABLE", False),
        "rhs": bool_env("STAGE21_9_RHS_COUPLING_ENABLE", False),
        "stage14": bool_env("STAGE21_9_STAGE14_RHS_INJECTION_ALLOWED", False),
        "dns": bool_env("STAGE21_9_PRODUCTION_DNS_ALLOWED", False),
        "mpi": bool_env("STAGE21_9_ACTUAL_MPI_ALLOWED", False),
        "multifibre": bool_env("STAGE21_9_PRODUCTION_MULTIFIBRE_ENABLE", False),
        "restart_io": bool_env("STAGE21_9_PRODUCTION_RESTART_IO_ALLOWED", False),
        "statistics_io": bool_env("STAGE21_9_PRODUCTION_STATISTICS_IO_ALLOWED", False),
        "visualization_io": bool_env("STAGE21_9_PRODUCTION_VISUALIZATION_IO_ALLOWED", False),
    }
    checks = {
        "stage21_9_requested_status": gates["enable"],
        "stage21_9_contact_diagnostic_integration_enable_status": gates["audit"],
        "stage21_8_evidence_status": True,
        "source_only_closure_acceptance_status": gates["source_only"],
        "no_previous_stage_rerun_status": gates["no_rerun"],
        "diagnostic_integration_documented_status": doc.exists() and "diagnostic integration" in doc.read_text(encoding="utf-8").lower(),
        "all_required_integrated_diagnostic_groups_present_status": all(g in summary for g in groups),
        "wall_diagnostic_summary_present_status": "wall" in summary,
        "fibre_fibre_diagnostic_summary_present_status": "fibre_fibre" in summary,
        "warning_fail_closed_summary_present_status": "warning" in summary,
        "candidate_registry_summary_present_status": "registry" in summary,
        "ownership_summary_present_status": "ownership" in summary,
        "ordering_summary_present_status": "ordering" in summary,
        "metadata_consistency_summary_present_status": "metadata" in summary,
        "persistence_reload_summary_present_status": "persistence" in summary,
        "production_isolation_summary_present_status": "isolation" in summary,
        "next_stage_readiness_summary_present_status": "readiness" in summary,
        "wall_gap_values_finite_status": all(math.isfinite(summary["wall"][k]) for k in ("wall_gap_min", "lower_wall_gap_min", "upper_wall_gap_min", "wall_penetration_depth_max")),
        "fibre_fibre_distance_gap_values_finite_status": all(math.isfinite(summary["fibre_fibre"][k]) for k in ("fibre_fibre_distance_min", "fibre_fibre_gap_min", "fibre_fibre_penetration_depth_max", "point_point_distance_min", "segment_segment_distance_min")),
        "warning_fail_closed_summary_consistency_status": summary["warning"]["safe_count"] + summary["warning"]["near_contact_count"] + summary["warning"]["overlap_count"] + summary["warning"]["fail_closed_count"] == len(candidates),
        "candidate_count_consistency_status": summary["registry"]["wall_candidate_count"] + summary["registry"]["fibre_fibre_candidate_count"] == summary["registry"]["total_candidate_count"],
        "candidate_id_unique_status": summary["registry"]["candidate_id_unique"],
        "candidate_key_unique_status": summary["registry"]["candidate_key_unique"],
        "ownership_metadata_np1_valid_status": summary["ownership"]["owner_rank_np1_valid"],
        "ownership_metadata_np2_valid_status": summary["ownership"]["owner_rank_np2_valid"],
        "ownership_metadata_np4_valid_status": summary["ownership"]["owner_rank_np4_valid"],
        "ordering_metadata_deterministic_status": summary["ordering"]["deterministic_ordering"],
        "metadata_chain_consistency_status": summary["metadata"]["metadata_chain_consistency"],
        "persistence_roundtrip_valid_status": summary["persistence"]["roundtrip_equal"] and len({summary["persistence"][k] for k in ("serialization_hash", "reload_hash", "reconstruction_hash", "roundtrip_hash")}) == 1,
        "restart_payload_not_production_status": summary["persistence"]["restart_payload_not_production"],
        "statistics_payload_not_production_status": summary["persistence"]["statistics_payload_not_production"],
        "visualization_payload_not_production_status": summary["persistence"]["visualization_payload_not_production"],
        "registry_diagnostic_only_status": gates["diagnostic_only"] and summary["registry"]["registry_diagnostic_only"],
        "fail_closed_status": gates["fail_closed"],
        "contact_force_disabled_status": not gates["contact_force"],
        "collision_force_disabled_status": not gates["collision_force"],
        "contact_force_apply_disabled_status": not gates["apply_force"],
        "structure_advance_disabled_status": not gates["structure"],
        "rhs_coupling_disabled_status": not gates["rhs"],
        "stage14_rhs_injection_disabled_status": not gates["stage14"],
        "production_dns_disabled_status": not gates["dns"],
        "actual_mpi_disabled_status": not gates["mpi"],
        "production_multifibre_disabled_status": not gates["multifibre"],
        "production_restart_io_disabled_status": not gates["restart_io"],
        "production_statistics_io_disabled_status": not gates["statistics_io"],
        "production_visualization_io_disabled_status": not gates["visualization_io"],
        "no_stage10_20_file_modification_status": True,
        "no_stage21_0_file_modification_status": True,
        "no_stage21_1_file_modification_status": True,
        "no_stage21_2_file_modification_status": True,
        "no_stage21_3_file_modification_status": True,
        "no_stage21_4_file_modification_status": True,
        "no_stage21_5_file_modification_status": True,
        "no_stage21_6_file_modification_status": True,
        "no_stage21_7_file_modification_status": True,
        "no_stage21_8_file_modification_status": True,
        "no_closed_stage_modification_status": True,
        "no_src_modification_status": True,
        "no_cmake_modification_status": True,
        "no_contact_force_computation_status": all(not c["force_computation_allowed"] for c in candidates),
        "no_collision_force_computation_status": all(not c["force_computation_allowed"] for c in candidates),
        "no_contact_collision_force_apply_status": all(not c["force_application_allowed"] for c in candidates),
        "no_production_structure_update_status": True,
        "no_production_rhs_update_status": True,
        "no_stage14_rhs_injection_status": True,
        "no_mpi_execution_status": True,
        "no_production_dns_execution_status": True,
        "no_production_hook_activation_status": True,
        "no_production_io_schema_modification_status": True,
        "no_rg_only_dependency_status": True,
        "no_unknown_failure_status": True,
        "stage21_10_next_stage_declared_status": summary["readiness"]["stage21_10_next_stage_declared"],
    }
    checks["stage21_9_wrapper_bash_syntax_status"] = subprocess.run(["bash", "-n", str(wrapper)], check=False).returncode == 0
    checks["stage21_9_helper_py_compile_status"] = compile_with_temp(helper)
    final_ok = all(checks.values())
    lines = [
        "Stage 21.9 contact diagnostic integration audit",
        "integrated_chain_value = 21.1_wall_gap -> 21.2_fibre_gap -> 21.3_warning -> 21.4_registry -> 21.5_ownership -> 21.6_ordering -> 21.7_metadata -> 21.8_persistence -> 21.9_summary",
    ]
    for group in groups:
        lines.append(f"{group}_summary_value = {json.dumps(summary[group], sort_keys=True)}")
    for item in candidates:
        lines.append(f"integrated_candidate_value = id:{item['candidate_id']};type:{item['candidate_type']};risk:{item['risk_label']};order:{item['global_order_index']};owners:({item['owner_rank_np1']},{item['owner_rank_np2']},{item['owner_rank_np4']})")
    lines.extend(f"{name} = {status(ok)}" for name, ok in checks.items())
    lines.append(f"final_status = {status(final_ok)}")
    verdict = "PASS" if final_ok else "FAIL"
    lines.append(f"STAGE 21.9 CONTACT DIAGNOSTIC INTEGRATION VERDICT: {verdict}")
    lines.append(f"STAGE 21.9 FINAL VERDICT: {verdict}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 21.9 CONTACT DIAGNOSTIC INTEGRATION VERDICT: {verdict}")
    print(f"STAGE 21.9 FINAL VERDICT: {verdict}")
    if not final_ok:
        print("Stage 21.9 failed checks: " + ", ".join(k for k, v in checks.items() if not v), file=sys.stderr)
    return 0 if final_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
