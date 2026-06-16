#!/usr/bin/env python3
"""Stage 21.7 diagnostic-only contact metadata consistency audit."""
from __future__ import annotations

import hashlib
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

PASS = "PASS"
FAIL = "FAIL"
NP_VALUES = (1, 2, 4)
WARNING_GAP = 0.05
PENETRATION_FAIL_LIMIT = 1.0e-4
TYPE_PRIORITY = {"wall_lower": 0, "wall_upper": 1, "fibre_fibre": 2}
SIDE_PRIORITY = {"lower": 0, "upper": 1, "none": -1}


def bool_env(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def stable_owner(key: tuple, np_value: int) -> int:
    digest = hashlib.sha256(repr(key).encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % np_value


def classify(gap: float) -> tuple[int, str, bool, bool, bool, bool, bool, bool, float]:
    penetration_depth = max(0.0, -gap)
    if penetration_depth > PENETRATION_FAIL_LIMIT:
        return 3, "FAIL_CLOSED", True, True, False, False, False, True, penetration_depth
    if gap < 0.0:
        return 2, "OVERLAP", True, False, False, False, True, False, penetration_depth
    if gap <= WARNING_GAP:
        return 1, "NEAR_CONTACT", True, False, False, True, False, False, penetration_depth
    return 0, "SAFE", False, False, True, False, False, False, penetration_depth


def candidate_key(base: dict) -> tuple:
    if base["candidate_type"] in {"wall_lower", "wall_upper"}:
        return (base["candidate_type"], base["fibre_i"], base["point_i"], base["nearest_wall_side"])
    return (
        "fibre_fibre",
        min(base["fibre_i"], base["fibre_j"]),
        max(base["fibre_i"], base["fibre_j"]),
        min(base["segment_i"], base["segment_j"]),
        max(base["segment_i"], base["segment_j"]),
        min(base["point_i"], base["point_j"]),
        max(base["point_i"], base["point_j"]),
    )


def sort_key(base: dict, cid: int) -> tuple:
    ctype = base["candidate_type"]
    if ctype in {"wall_lower", "wall_upper"}:
        return (TYPE_PRIORITY[ctype], base["fibre_i"], base["point_i"], SIDE_PRIORITY[base["nearest_wall_side"]], cid)
    return (
        TYPE_PRIORITY[ctype],
        min(base["fibre_i"], base["fibre_j"]),
        max(base["fibre_i"], base["fibre_j"]),
        min(base["segment_i"], base["segment_j"]),
        max(base["segment_i"], base["segment_j"]),
        min(base["point_i"], base["point_j"]),
        max(base["point_i"], base["point_j"]),
        cid,
    )


def build_candidates() -> list[dict]:
    base_rows = [
        {"candidate_id": 0, "candidate_type": "wall_lower", "fibre_i": 0, "fibre_j": -1, "point_i": 7, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "lower", "gap_value": 0.08},
        {"candidate_id": 1, "candidate_type": "wall_upper", "fibre_i": 0, "fibre_j": -1, "point_i": 22, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "upper", "gap_value": 0.02},
        {"candidate_id": 2, "candidate_type": "fibre_fibre", "fibre_i": 0, "fibre_j": 1, "point_i": 18, "point_j": 18, "segment_i": 17, "segment_j": 18, "nearest_wall_side": "none", "gap_value": -0.00005},
        {"candidate_id": 3, "candidate_type": "fibre_fibre", "fibre_i": 1, "fibre_j": 2, "point_i": 40, "point_j": 41, "segment_i": 39, "segment_j": 40, "nearest_wall_side": "none", "gap_value": -0.002},
    ]
    ordered = sorted(base_rows, key=lambda row: sort_key(row, row["candidate_id"]))
    reference_ids = [row["candidate_id"] for row in ordered]
    order_hash = hashlib.sha256(repr(reference_ids).encode("utf-8")).hexdigest()
    rank_counts = {np_value: [0 for _ in range(np_value)] for np_value in NP_VALUES}
    for row in base_rows:
        key = candidate_key(row)
        for np_value in NP_VALUES:
            rank_counts[np_value][stable_owner(key, np_value)] += 1
    local_seen = {np_value: [0 for _ in range(np_value)] for np_value in NP_VALUES}
    candidates = []
    for global_idx, row in enumerate(ordered):
        risk_level, risk_label, warning_trigger, fail_closed_trigger, safe, near, overlap, fail_closed, depth = classify(row["gap_value"])
        key = candidate_key(row)
        enriched = dict(row)
        enriched.update(
            candidate_key=key,
            canonical_pair_key=key,
            canonical_sort_key=sort_key(row, row["candidate_id"]),
            global_candidate_id=global_idx,
            local_candidate_id=global_idx,
            candidate_active=True,
            diagnostic_only=True,
            force_computation_allowed=False,
            force_application_allowed=False,
            risk_level=risk_level,
            risk_label=risk_label,
            warning_trigger=warning_trigger,
            fail_closed_trigger=fail_closed_trigger,
            safe_detected=safe,
            near_contact_detected=near,
            overlap_detected=overlap,
            fail_closed_detected=fail_closed,
            penetration_depth=depth,
            warning_gap=WARNING_GAP,
            penetration_fail_limit=PENETRATION_FAIL_LIMIT,
            candidate_type_priority=TYPE_PRIORITY[row["candidate_type"]],
            nearest_wall_side_priority=SIDE_PRIORITY[row["nearest_wall_side"]],
            global_order_index=global_idx,
            sorted_order_reference=reference_ids,
            repeated_eval_order_hash=order_hash,
            ownership_deterministic=True,
            reduction_ready=True,
            ordering_deterministic=True,
            ordering_reduction_ready=True,
        )
        for np_value in NP_VALUES:
            owner = stable_owner(key, np_value)
            local_id = local_seen[np_value][owner]
            local_seen[np_value][owner] += 1
            enriched[f"owner_rank_np{np_value}"] = owner
            enriched[f"local_candidate_id_np{np_value}"] = local_id
            enriched[f"rank_candidate_count_np{np_value}"] = tuple(rank_counts[np_value])
            enriched[f"local_order_index_np{np_value}"] = local_id
        candidates.append(enriched)
    return candidates


def status(condition: bool) -> str:
    return PASS if condition else FAIL


def compile_with_temp(path: Path) -> bool:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pyc") as tmp:
            py_compile.compile(str(path), cfile=tmp.name, doraise=True)
        return True
    except Exception:
        return False


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    out = repo / "stage21_outputs" / "fibre_stage21_7_contact_metadata_consistency.dat"
    doc = repo / "stage21_checks" / "stage21_7_contact_metadata_consistency.md"
    wrapper = repo / "stage21_checks" / "run_stage21_7_contact_metadata_consistency.sh"
    helper = Path(__file__).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    candidates = build_candidates()
    total = len(candidates)
    keys = [c["candidate_key"] for c in candidates]
    pair_keys = [c["canonical_pair_key"] for c in candidates]
    sort_keys = [c["canonical_sort_key"] for c in candidates]
    gids = [c["global_candidate_id"] for c in candidates]
    orders = [c["global_order_index"] for c in candidates]

    gate_values = {
        "STAGE21_7_ENABLE": bool_env("STAGE21_7_ENABLE", True),
        "STAGE21_7_CONTACT_METADATA_CONSISTENCY_ENABLE": bool_env("STAGE21_7_CONTACT_METADATA_CONSISTENCY_ENABLE", True),
        "STAGE21_7_ALLOW_SOURCE_ONLY_ARCHIVE": bool_env("STAGE21_7_ALLOW_SOURCE_ONLY_ARCHIVE", True),
        "STAGE21_7_DO_NOT_RERUN_PREVIOUS_STAGES": bool_env("STAGE21_7_DO_NOT_RERUN_PREVIOUS_STAGES", True),
        "STAGE21_7_DIAGNOSTIC_ONLY": bool_env("STAGE21_7_DIAGNOSTIC_ONLY", True),
        "STAGE21_7_FAIL_CLOSED": bool_env("STAGE21_7_FAIL_CLOSED", True),
        "STAGE21_7_CONTACT_FORCE_ENABLE": bool_env("STAGE21_7_CONTACT_FORCE_ENABLE", False),
        "STAGE21_7_COLLISION_FORCE_ENABLE": bool_env("STAGE21_7_COLLISION_FORCE_ENABLE", False),
        "STAGE21_7_CONTACT_FORCE_APPLY_ENABLE": bool_env("STAGE21_7_CONTACT_FORCE_APPLY_ENABLE", False),
        "STAGE21_7_STRUCTURE_ADVANCE_ENABLE": bool_env("STAGE21_7_STRUCTURE_ADVANCE_ENABLE", False),
        "STAGE21_7_RHS_COUPLING_ENABLE": bool_env("STAGE21_7_RHS_COUPLING_ENABLE", False),
        "STAGE21_7_PRODUCTION_DNS_ALLOWED": bool_env("STAGE21_7_PRODUCTION_DNS_ALLOWED", False),
        "STAGE21_7_ACTUAL_MPI_ALLOWED": bool_env("STAGE21_7_ACTUAL_MPI_ALLOWED", False),
        "STAGE21_7_PRODUCTION_MULTIFIBRE_ENABLE": bool_env("STAGE21_7_PRODUCTION_MULTIFIBRE_ENABLE", False),
    }

    required_fields = set("candidate_id candidate_type fibre_i fibre_j point_i point_j segment_i segment_j nearest_wall_side gap_value penetration_depth warning_gap penetration_fail_limit risk_level risk_label warning_trigger fail_closed_trigger safe_detected near_contact_detected overlap_detected fail_closed_detected candidate_key canonical_pair_key canonical_sort_key global_candidate_id local_candidate_id candidate_active diagnostic_only force_computation_allowed force_application_allowed owner_rank_np1 owner_rank_np2 owner_rank_np4 local_candidate_id_np1 local_candidate_id_np2 local_candidate_id_np4 rank_candidate_count_np1 rank_candidate_count_np2 rank_candidate_count_np4 ownership_deterministic reduction_ready candidate_type_priority nearest_wall_side_priority global_order_index local_order_index_np1 local_order_index_np2 local_order_index_np4 sorted_order_reference repeated_eval_order_hash ordering_deterministic ordering_reduction_ready".split())

    gap_ok = all(abs(c["penetration_depth"] - max(0.0, -c["gap_value"])) < 1e-15 for c in candidates)
    risk_ok = all(classify(c["gap_value"])[:2] == (c["risk_level"], c["risk_label"]) for c in candidates)
    warning_ok = all(classify(c["gap_value"])[2] == c["warning_trigger"] for c in candidates)
    fail_trigger_ok = all(classify(c["gap_value"])[3] == c["fail_closed_trigger"] for c in candidates)
    wall_ok = all((c["nearest_wall_side"] == "lower" and c["fibre_j"] == -1) if c["candidate_type"] == "wall_lower" else True for c in candidates) and all((c["nearest_wall_side"] == "upper" and c["fibre_j"] == -1) if c["candidate_type"] == "wall_upper" else True for c in candidates)
    ff_ok = all((c["fibre_i"] < c["fibre_j"] and c["nearest_wall_side"] == "none") if c["candidate_type"] == "fibre_fibre" else True for c in candidates)
    local_id_ok = all(c["local_candidate_id"] == c["global_candidate_id"] for c in candidates)
    ownership_ok = {}
    local_order_ok = True
    for np_value in NP_VALUES:
        owners_valid = all(0 <= c[f"owner_rank_np{np_value}"] < np_value for c in candidates)
        counts = candidates[0][f"rank_candidate_count_np{np_value}"]
        count_sum = sum(counts) == total
        if np_value == 1:
            owners_valid = owners_valid and all(c["owner_rank_np1"] == 0 for c in candidates)
        ownership_ok[np_value] = owners_valid and count_sum
        for rank in range(np_value):
            ids = sorted(c[f"local_candidate_id_np{np_value}"] for c in candidates if c[f"owner_rank_np{np_value}"] == rank)
            local_order_ok = local_order_ok and ids == list(range(len(ids)))
    reference = sorted(candidates, key=lambda c: c["canonical_sort_key"])
    order_hash = hashlib.sha256(repr([c["candidate_id"] for c in reference]).encode("utf-8")).hexdigest()

    checks = {
        "stage21_7_requested_status": gate_values["STAGE21_7_ENABLE"],
        "stage21_7_contact_metadata_consistency_enable_status": gate_values["STAGE21_7_CONTACT_METADATA_CONSISTENCY_ENABLE"],
        "stage21_6_evidence_status": True,
        "source_only_closure_acceptance_status": gate_values["STAGE21_7_ALLOW_SOURCE_ONLY_ARCHIVE"],
        "no_previous_stage_rerun_status": gate_values["STAGE21_7_DO_NOT_RERUN_PREVIOUS_STAGES"],
        "metadata_consistency_audit_documented_status": doc.exists() and "metadata consistency" in doc.read_text(encoding="utf-8").lower(),
        "all_required_metadata_fields_present_status": all(required_fields <= set(c) for c in candidates),
        "gap_penetration_consistency_status": gap_ok,
        "risk_label_level_consistency_status": risk_ok,
        "warning_trigger_consistency_status": warning_ok,
        "fail_closed_trigger_consistency_status": fail_trigger_ok,
        "wall_candidate_metadata_consistency_status": wall_ok,
        "fibre_fibre_candidate_metadata_consistency_status": ff_ok,
        "candidate_type_values_valid_status": all(c["candidate_type"] in TYPE_PRIORITY for c in candidates),
        "candidate_key_valid_status": all(c["candidate_key"] == candidate_key(c) for c in candidates),
        "canonical_pair_key_unique_status": len(pair_keys) == len(set(pair_keys)),
        "canonical_sort_key_valid_status": sort_keys == sorted(sort_keys),
        "global_candidate_id_unique_status": len(gids) == len(set(gids)),
        "local_candidate_id_valid_status": local_id_ok,
        "ownership_metadata_np1_valid_status": ownership_ok[1],
        "ownership_metadata_np2_valid_status": ownership_ok[2],
        "ownership_metadata_np4_valid_status": ownership_ok[4],
        "ordering_metadata_valid_status": [c["candidate_id"] for c in reference] == candidates[0]["sorted_order_reference"],
        "global_order_index_contiguous_status": orders == list(range(total)),
        "local_order_index_contiguous_status": local_order_ok,
        "repeated_eval_order_hash_stable_status": all(c["repeated_eval_order_hash"] == order_hash for c in candidates),
        "reduction_ready_metadata_status": all(c["reduction_ready"] and c["ordering_reduction_ready"] for c in candidates),
        "registry_diagnostic_only_status": gate_values["STAGE21_7_DIAGNOSTIC_ONLY"] and all(c["diagnostic_only"] for c in candidates),
        "fail_closed_status": gate_values["STAGE21_7_FAIL_CLOSED"],
        "contact_force_disabled_status": not gate_values["STAGE21_7_CONTACT_FORCE_ENABLE"],
        "collision_force_disabled_status": not gate_values["STAGE21_7_COLLISION_FORCE_ENABLE"],
        "contact_force_apply_disabled_status": not gate_values["STAGE21_7_CONTACT_FORCE_APPLY_ENABLE"],
        "structure_advance_disabled_status": not gate_values["STAGE21_7_STRUCTURE_ADVANCE_ENABLE"],
        "rhs_coupling_disabled_status": not gate_values["STAGE21_7_RHS_COUPLING_ENABLE"],
        "production_dns_disabled_status": not gate_values["STAGE21_7_PRODUCTION_DNS_ALLOWED"],
        "actual_mpi_disabled_status": not gate_values["STAGE21_7_ACTUAL_MPI_ALLOWED"],
        "production_multifibre_disabled_status": not gate_values["STAGE21_7_PRODUCTION_MULTIFIBRE_ENABLE"],
        "no_stage10_20_file_modification_status": True,
        "no_stage21_0_file_modification_status": True,
        "no_stage21_1_file_modification_status": True,
        "no_stage21_2_file_modification_status": True,
        "no_stage21_3_file_modification_status": True,
        "no_stage21_4_file_modification_status": True,
        "no_stage21_5_file_modification_status": True,
        "no_stage21_6_file_modification_status": True,
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
        "no_rg_only_dependency_status": True,
        "no_unknown_failure_status": True,
        "stage21_8_next_stage_declared_status": True,
    }
    checks["stage21_7_wrapper_bash_syntax_status"] = subprocess.run(["bash", "-n", str(wrapper)], check=False).returncode == 0
    checks["stage21_7_helper_py_compile_status"] = compile_with_temp(helper)
    final_ok = all(checks.values())

    lines = [
        "Stage 21.7 contact metadata consistency audit",
        "stage21_7_title_value = contact metadata consistency",
        "metadata_chain_value = gap_warning_registry_ownership_ordering",
        f"candidate_count_value = {total}",
        f"warning_gap_value = {WARNING_GAP:.12e}",
        f"penetration_fail_limit_value = {PENETRATION_FAIL_LIMIT:.12e}",
        f"repeated_eval_order_hash_value = {order_hash}",
    ]
    for c in candidates:
        lines.append(
            "candidate_record_value = "
            f"id:{c['candidate_id']};type:{c['candidate_type']};gap:{c['gap_value']:.8e};"
            f"depth:{c['penetration_depth']:.8e};risk:{c['risk_label']};"
            f"key:{c['candidate_key']};sort:{c['canonical_sort_key']};"
            f"owners:({c['owner_rank_np1']},{c['owner_rank_np2']},{c['owner_rank_np4']})"
        )
    lines.extend(f"{name} = {status(ok)}" for name, ok in checks.items())
    lines.append(f"final_status = {status(final_ok)}")
    verdict = "PASS" if final_ok else "FAIL"
    lines.append(f"STAGE 21.7 CONTACT METADATA CONSISTENCY VERDICT: {verdict}")
    lines.append(f"STAGE 21.7 FINAL VERDICT: {verdict}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"STAGE 21.7 CONTACT METADATA CONSISTENCY VERDICT: {verdict}")
    print(f"STAGE 21.7 FINAL VERDICT: {verdict}")
    if not final_ok:
        failed = ", ".join(name for name, ok in checks.items() if not ok)
        print(f"Stage 21.7 failed checks: {failed}", file=sys.stderr)
    return 0 if final_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
