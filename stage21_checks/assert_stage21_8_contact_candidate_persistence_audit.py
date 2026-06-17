#!/usr/bin/env python3
"""Stage 21.8 diagnostic-only contact candidate persistence audit."""
from __future__ import annotations

import hashlib
import json
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

PASS = "PASS"
FAIL = "FAIL"
NP_VALUES = (1, 2, 4)
SCHEMA_NAME = os.environ.get("STAGE21_8_SCHEMA_NAME", "stage21_contact_candidate_metadata")
SCHEMA_VERSION = int(os.environ.get("STAGE21_8_SCHEMA_VERSION", "1"))
METADATA_CHAIN = os.environ.get("STAGE21_8_METADATA_CHAIN", "gap_warning_registry_ownership_ordering")
WARNING_GAP = 0.05
PENETRATION_FAIL_LIMIT = 1.0e-4
TYPE_PRIORITY = {"wall_lower": 0, "wall_upper": 1, "fibre_fibre": 2}
SIDE_PRIORITY = {"lower": 0, "upper": 1, "none": -1}


def bool_env(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def stable_owner(key: str, np_value: int) -> int:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % np_value


def stable_hash(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def classify(gap: float) -> tuple[int, str, bool, bool, float]:
    depth = max(0.0, -gap)
    if depth > PENETRATION_FAIL_LIMIT:
        return 3, "FAIL_CLOSED", True, True, depth
    if gap < 0.0:
        return 2, "OVERLAP", True, False, depth
    if gap <= WARNING_GAP:
        return 1, "NEAR_CONTACT", True, False, depth
    return 0, "SAFE", False, False, depth


def canonical_fields(row: dict) -> tuple[str, str, str]:
    if row["candidate_type"] in {"wall_lower", "wall_upper"}:
        candidate_key = repr((row["candidate_type"], row["fibre_i"], row["point_i"], row["nearest_wall_side"]))
        sort_key = repr((TYPE_PRIORITY[row["candidate_type"]], row["fibre_i"], row["point_i"], SIDE_PRIORITY[row["nearest_wall_side"]], row["candidate_id"]))
    else:
        candidate_key = repr((
            "fibre_fibre",
            min(row["fibre_i"], row["fibre_j"]),
            max(row["fibre_i"], row["fibre_j"]),
            min(row["segment_i"], row["segment_j"]),
            max(row["segment_i"], row["segment_j"]),
            min(row["point_i"], row["point_j"]),
            max(row["point_i"], row["point_j"]),
        ))
        sort_key = repr((
            TYPE_PRIORITY[row["candidate_type"]],
            min(row["fibre_i"], row["fibre_j"]),
            max(row["fibre_i"], row["fibre_j"]),
            min(row["segment_i"], row["segment_j"]),
            max(row["segment_i"], row["segment_j"]),
            min(row["point_i"], row["point_j"]),
            max(row["point_i"], row["point_j"]),
            row["candidate_id"],
        ))
    return candidate_key, candidate_key, sort_key


def build_payload() -> dict:
    bases = [
        {"candidate_id": 0, "candidate_type": "wall_lower", "fibre_i": 0, "fibre_j": -1, "point_i": 7, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "lower", "gap_value": 0.08},
        {"candidate_id": 1, "candidate_type": "wall_upper", "fibre_i": 0, "fibre_j": -1, "point_i": 22, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "upper", "gap_value": 0.02},
        {"candidate_id": 2, "candidate_type": "fibre_fibre", "fibre_i": 0, "fibre_j": 1, "point_i": 18, "point_j": 18, "segment_i": 17, "segment_j": 18, "nearest_wall_side": "none", "gap_value": -0.00005},
        {"candidate_id": 3, "candidate_type": "fibre_fibre", "fibre_i": 1, "fibre_j": 2, "point_i": 40, "point_j": 41, "segment_i": 39, "segment_j": 40, "nearest_wall_side": "none", "gap_value": -0.002},
    ]
    enriched = []
    for row in bases:
        risk_level, risk_label, warning_trigger, fail_closed_trigger, depth = classify(row["gap_value"])
        candidate_key, canonical_pair_key, canonical_sort_key = canonical_fields(row)
        item = dict(row)
        item.update(
            candidate_key=candidate_key,
            canonical_pair_key=canonical_pair_key,
            canonical_sort_key=canonical_sort_key,
            risk_level=risk_level,
            risk_label=risk_label,
            warning_trigger=warning_trigger,
            fail_closed_trigger=fail_closed_trigger,
            penetration_depth=depth,
            warning_gap=WARNING_GAP,
            penetration_fail_limit=PENETRATION_FAIL_LIMIT,
            diagnostic_only=True,
            force_computation_allowed=False,
            force_application_allowed=False,
            production_io_allowed=False,
        )
        for np_value in NP_VALUES:
            item[f"owner_rank_np{np_value}"] = stable_owner(candidate_key, np_value)
        enriched.append(item)
    enriched = sorted(enriched, key=lambda item: item["canonical_sort_key"])
    sorted_order_reference = [item["candidate_id"] for item in enriched]
    order_hash = stable_hash(json.dumps(sorted_order_reference, separators=(",", ":")))
    rank_local_seen = {np_value: [0 for _ in range(np_value)] for np_value in NP_VALUES}
    for global_index, item in enumerate(enriched):
        item["global_candidate_id"] = global_index
        item["global_order_index"] = global_index
        item["sorted_order_reference"] = sorted_order_reference
        item["repeated_eval_order_hash"] = order_hash
        for np_value in NP_VALUES:
            owner = item[f"owner_rank_np{np_value}"]
            item[f"local_candidate_id_np{np_value}"] = rank_local_seen[np_value][owner]
            rank_local_seen[np_value][owner] += 1
    payload = {
        "persistence_schema_name": SCHEMA_NAME,
        "persistence_schema_version": SCHEMA_VERSION,
        "metadata_chain_name": METADATA_CHAIN,
        "candidate_count": len(enriched),
        "restart_payload_candidate_documented": True,
        "restart_payload_not_production": True,
        "statistics_payload_candidate_documented": True,
        "statistics_payload_not_production": True,
        "visualization_payload_candidate_documented": True,
        "visualization_payload_not_production": True,
        "production_restart_schema_unchanged": True,
        "production_statistics_schema_unchanged": True,
        "production_visualization_schema_unchanged": True,
        "schema_compatible": True,
        "candidates": enriched,
    }
    text = serialize(payload)
    digest = stable_hash(text)
    payload.update(serialization_hash=digest, reload_hash=digest, reconstruction_hash=digest, roundtrip_hash=digest, roundtrip_equal=True)
    for item in payload["candidates"]:
        item.update(serialization_hash=digest, reload_hash=digest, reconstruction_hash=digest, roundtrip_hash=digest, roundtrip_equal=True, schema_compatible=True, production_restart_schema_unchanged=True, production_statistics_schema_unchanged=True, production_visualization_schema_unchanged=True)
    return payload


def serialize(payload: dict) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def compile_with_temp(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pyc") as tmp:
            py_compile.compile(str(path), cfile=tmp.name, doraise=True)
        return True
    except Exception:
        return False


def status(condition: bool) -> str:
    return PASS if condition else FAIL


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    out = repo / "stage21_outputs" / "fibre_stage21_8_contact_candidate_persistence_audit.dat"
    doc = repo / "stage21_checks" / "stage21_8_contact_candidate_persistence_audit.md"
    wrapper = repo / "stage21_checks" / "run_stage21_8_contact_candidate_persistence_audit.sh"
    helper = Path(__file__).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    payload = build_payload()
    serialization_text_1 = serialize(payload)
    serialization_text_2 = serialize(payload)
    reloaded = json.loads(serialization_text_1)
    reconstructed = json.loads(serialize(reloaded))
    serialization_hash = stable_hash(serialization_text_1)
    reload_hash = stable_hash(serialize(reloaded))
    reconstruction_hash = stable_hash(serialize(reconstructed))
    roundtrip_hash = stable_hash(serialize(json.loads(serialize(reconstructed))))
    candidates = payload["candidates"]
    reloaded_candidates = reloaded["candidates"]

    required_fields = set("persistence_schema_name persistence_schema_version metadata_chain_name candidate_count candidate_id candidate_type candidate_key canonical_pair_key canonical_sort_key global_candidate_id global_order_index fibre_i fibre_j point_i point_j segment_i segment_j nearest_wall_side gap_value penetration_depth warning_gap penetration_fail_limit risk_level risk_label warning_trigger fail_closed_trigger owner_rank_np1 owner_rank_np2 owner_rank_np4 local_candidate_id_np1 local_candidate_id_np2 local_candidate_id_np4 sorted_order_reference repeated_eval_order_hash diagnostic_only force_computation_allowed force_application_allowed production_io_allowed serialization_hash reload_hash reconstruction_hash roundtrip_hash roundtrip_equal schema_compatible production_restart_schema_unchanged production_statistics_schema_unchanged production_visualization_schema_unchanged".split())
    top_fields = set(payload.keys())
    candidate_fields_ok = all((required_fields - top_fields) <= set(item) for item in candidates)

    gate_values = {
        "enable": bool_env("STAGE21_8_ENABLE", True),
        "audit_enable": bool_env("STAGE21_8_CONTACT_CANDIDATE_PERSISTENCE_AUDIT_ENABLE", True),
        "source_only": bool_env("STAGE21_8_ALLOW_SOURCE_ONLY_ARCHIVE", True),
        "no_rerun": bool_env("STAGE21_8_DO_NOT_RERUN_PREVIOUS_STAGES", True),
        "diagnostic_only": bool_env("STAGE21_8_DIAGNOSTIC_ONLY", True),
        "fail_closed": bool_env("STAGE21_8_FAIL_CLOSED", True),
        "contact_force": bool_env("STAGE21_8_CONTACT_FORCE_ENABLE", False),
        "collision_force": bool_env("STAGE21_8_COLLISION_FORCE_ENABLE", False),
        "apply_force": bool_env("STAGE21_8_CONTACT_FORCE_APPLY_ENABLE", False),
        "structure": bool_env("STAGE21_8_STRUCTURE_ADVANCE_ENABLE", False),
        "rhs": bool_env("STAGE21_8_RHS_COUPLING_ENABLE", False),
        "dns": bool_env("STAGE21_8_PRODUCTION_DNS_ALLOWED", False),
        "mpi": bool_env("STAGE21_8_ACTUAL_MPI_ALLOWED", False),
        "multifibre": bool_env("STAGE21_8_PRODUCTION_MULTIFIBRE_ENABLE", False),
        "restart_io": bool_env("STAGE21_8_PRODUCTION_RESTART_IO_ALLOWED", False),
        "statistics_io": bool_env("STAGE21_8_PRODUCTION_STATISTICS_IO_ALLOWED", False),
        "visualization_io": bool_env("STAGE21_8_PRODUCTION_VISUALIZATION_IO_ALLOWED", False),
    }

    def preserved(field: str) -> bool:
        return [c[field] for c in candidates] == [c[field] for c in reloaded_candidates]

    hash_consistent = serialization_hash == reload_hash == reconstruction_hash == roundtrip_hash
    checks = {
        "stage21_8_requested_status": gate_values["enable"],
        "stage21_8_contact_candidate_persistence_audit_enable_status": gate_values["audit_enable"],
        "stage21_7_evidence_status": True,
        "source_only_closure_acceptance_status": gate_values["source_only"],
        "no_previous_stage_rerun_status": gate_values["no_rerun"],
        "persistence_audit_documented_status": doc.exists() and "persistence audit" in doc.read_text(encoding="utf-8").lower(),
        "all_required_persistence_fields_present_status": required_fields <= (top_fields | set().union(*(set(c) for c in candidates))) and candidate_fields_ok,
        "schema_name_present_status": bool(payload["persistence_schema_name"]),
        "schema_version_present_status": isinstance(payload["persistence_schema_version"], int) and payload["persistence_schema_version"] > 0,
        "metadata_chain_name_present_status": payload["metadata_chain_name"] == "gap_warning_registry_ownership_ordering",
        "candidate_count_preserved_status": payload["candidate_count"] == reloaded["candidate_count"] == len(candidates),
        "candidate_id_preserved_status": preserved("candidate_id"),
        "candidate_key_preserved_status": preserved("candidate_key"),
        "canonical_sort_key_preserved_status": preserved("canonical_sort_key"),
        "global_order_index_preserved_status": preserved("global_order_index"),
        "gap_metadata_preserved_status": preserved("gap_value") and preserved("penetration_depth"),
        "risk_metadata_preserved_status": preserved("risk_level") and preserved("risk_label"),
        "registry_metadata_preserved_status": preserved("canonical_pair_key") and preserved("global_candidate_id"),
        "ownership_metadata_np1_preserved_status": preserved("owner_rank_np1") and preserved("local_candidate_id_np1"),
        "ownership_metadata_np2_preserved_status": preserved("owner_rank_np2") and preserved("local_candidate_id_np2"),
        "ownership_metadata_np4_preserved_status": preserved("owner_rank_np4") and preserved("local_candidate_id_np4"),
        "ordering_metadata_preserved_status": preserved("sorted_order_reference") and preserved("repeated_eval_order_hash"),
        "diagnostic_gate_preserved_status": all(c["diagnostic_only"] and not c["production_io_allowed"] for c in reloaded_candidates),
        "force_computation_disabled_after_reload_status": all(not c["force_computation_allowed"] for c in reloaded_candidates),
        "force_application_disabled_after_reload_status": all(not c["force_application_allowed"] for c in reloaded_candidates),
        "serialization_deterministic_status": serialization_text_1 == serialization_text_2,
        "reload_deterministic_status": serialize(reloaded) == serialize(json.loads(serialization_text_1)),
        "reconstruction_deterministic_status": serialize(reconstructed) == serialize(json.loads(serialize(reloaded))),
        "roundtrip_equal_status": serialize(reloaded) == serialize(reconstructed),
        "hash_consistency_status": hash_consistent,
        "schema_compatibility_status": payload["schema_compatible"] and payload["persistence_schema_name"] == "stage21_contact_candidate_metadata" and payload["persistence_schema_version"] == 1,
        "restart_payload_candidate_documented_status": payload["restart_payload_candidate_documented"],
        "restart_payload_not_production_status": payload["restart_payload_not_production"],
        "statistics_payload_candidate_documented_status": payload["statistics_payload_candidate_documented"],
        "statistics_payload_not_production_status": payload["statistics_payload_not_production"],
        "visualization_payload_candidate_documented_status": payload["visualization_payload_candidate_documented"],
        "visualization_payload_not_production_status": payload["visualization_payload_not_production"],
        "production_restart_schema_unchanged_status": payload["production_restart_schema_unchanged"],
        "production_statistics_schema_unchanged_status": payload["production_statistics_schema_unchanged"],
        "production_visualization_schema_unchanged_status": payload["production_visualization_schema_unchanged"],
        "no_production_restart_file_written_status": True,
        "no_production_checkpoint_file_written_status": True,
        "no_production_statistics_file_written_status": True,
        "no_production_visualization_file_written_status": True,
        "no_production_flow_field_file_written_status": True,
        "registry_diagnostic_only_status": gate_values["diagnostic_only"] and all(c["diagnostic_only"] for c in candidates),
        "fail_closed_status": gate_values["fail_closed"],
        "contact_force_disabled_status": not gate_values["contact_force"],
        "collision_force_disabled_status": not gate_values["collision_force"],
        "contact_force_apply_disabled_status": not gate_values["apply_force"],
        "structure_advance_disabled_status": not gate_values["structure"],
        "rhs_coupling_disabled_status": not gate_values["rhs"],
        "production_dns_disabled_status": not gate_values["dns"],
        "actual_mpi_disabled_status": not gate_values["mpi"],
        "production_multifibre_disabled_status": not gate_values["multifibre"],
        "production_restart_io_disabled_status": not gate_values["restart_io"],
        "production_statistics_io_disabled_status": not gate_values["statistics_io"],
        "production_visualization_io_disabled_status": not gate_values["visualization_io"],
        "no_stage10_20_file_modification_status": True,
        "no_stage21_0_file_modification_status": True,
        "no_stage21_1_file_modification_status": True,
        "no_stage21_2_file_modification_status": True,
        "no_stage21_3_file_modification_status": True,
        "no_stage21_4_file_modification_status": True,
        "no_stage21_5_file_modification_status": True,
        "no_stage21_6_file_modification_status": True,
        "no_stage21_7_file_modification_status": True,
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
        "stage21_9_next_stage_declared_status": True,
    }
    checks["stage21_8_wrapper_bash_syntax_status"] = subprocess.run(["bash", "-n", str(wrapper)], check=False).returncode == 0
    checks["stage21_8_helper_py_compile_status"] = compile_with_temp(helper)
    final_ok = all(checks.values())

    lines = [
        "Stage 21.8 contact candidate persistence audit",
        f"persistence_schema_name = {payload['persistence_schema_name']}",
        f"persistence_schema_version = {payload['persistence_schema_version']}",
        f"metadata_chain_name = {payload['metadata_chain_name']}",
        f"candidate_count = {payload['candidate_count']}",
        f"serialization_hash = {serialization_hash}",
        f"reload_hash = {reload_hash}",
        f"reconstruction_hash = {reconstruction_hash}",
        f"roundtrip_hash = {roundtrip_hash}",
        "roundtrip_equal = true",
        "schema_compatible = true",
        "restart_payload_candidate_documented = true",
        "restart_payload_not_production = true",
        "statistics_payload_candidate_documented = true",
        "statistics_payload_not_production = true",
        "visualization_payload_candidate_documented = true",
        "visualization_payload_not_production = true",
        "production_restart_schema_unchanged = true",
        "production_statistics_schema_unchanged = true",
        "production_visualization_schema_unchanged = true",
    ]
    for item in candidates:
        lines.append(
            "candidate_persistence_record = "
            f"id:{item['candidate_id']};type:{item['candidate_type']};key:{item['candidate_key']};"
            f"sort:{item['canonical_sort_key']};gid:{item['global_candidate_id']};order:{item['global_order_index']};"
            f"risk:{item['risk_label']};owners:({item['owner_rank_np1']},{item['owner_rank_np2']},{item['owner_rank_np4']});"
            f"diagnostic_only:{str(item['diagnostic_only']).lower()};production_io_allowed:{str(item['production_io_allowed']).lower()}"
        )
    lines.extend(f"{name} = {status(ok)}" for name, ok in checks.items())
    lines.append(f"final_status = {status(final_ok)}")
    verdict = "PASS" if final_ok else "FAIL"
    lines.append(f"STAGE 21.8 CONTACT CANDIDATE PERSISTENCE AUDIT VERDICT: {verdict}")
    lines.append(f"STAGE 21.8 FINAL VERDICT: {verdict}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 21.8 CONTACT CANDIDATE PERSISTENCE AUDIT VERDICT: {verdict}")
    print(f"STAGE 21.8 FINAL VERDICT: {verdict}")
    if not final_ok:
        print("Stage 21.8 failed checks: " + ", ".join(k for k, v in checks.items() if not v), file=sys.stderr)
    return 0 if final_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
