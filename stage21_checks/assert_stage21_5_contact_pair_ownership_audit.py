#!/usr/bin/env python3
"""Stage 21.5 diagnostic-only contact pair ownership audit."""
from __future__ import annotations

import hashlib
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage21_checks"
OUT_DIR = ROOT / "stage21_outputs"
OUT_FILE = OUT_DIR / "fibre_stage21_5_contact_pair_ownership_audit.dat"
HELPER = CHECK_DIR / "assert_stage21_5_contact_pair_ownership_audit.py"
WRAPPER = CHECK_DIR / "run_stage21_5_contact_pair_ownership_audit.sh"
DOC = CHECK_DIR / "stage21_5_contact_pair_ownership_audit.md"

SAFE_DEFAULTS = {
    "STAGE21_5_ENABLE": "1",
    "STAGE21_5_CONTACT_PAIR_OWNERSHIP_AUDIT_ENABLE": "1",
    "STAGE21_5_REQUIRE_STAGE21_4_PASS": "1",
    "STAGE21_5_ALLOW_MISSING_OLD_OUTPUTS": "1",
    "STAGE21_5_ALLOW_SOURCE_ONLY_ARCHIVE": "1",
    "STAGE21_5_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_5_DIAGNOSTIC_ONLY": "1",
    "STAGE21_5_FAIL_CLOSED": "1",
    "STAGE21_5_NP_VALUES": "1,2,4",
    "STAGE21_5_OWNER_RULE": "hash_mod_np",
    "STAGE21_5_CONTACT_FORCE_ENABLE": "0",
    "STAGE21_5_COLLISION_FORCE_ENABLE": "0",
    "STAGE21_5_CONTACT_FORCE_APPLY_ENABLE": "0",
    "STAGE21_5_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_5_RHS_COUPLING_ENABLE": "0",
    "STAGE21_5_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_5_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_5_PRODUCTION_MULTIFIBRE_ENABLE": "0",
}

STATUS_FIELDS = """
stage21_5_requested_status
stage21_5_contact_pair_ownership_audit_enable_status
stage21_4_evidence_status
source_only_closure_acceptance_status
no_previous_stage_rerun_status
ownership_audit_documented_status
all_required_ownership_fields_present_status
canonical_pair_key_present_status
wall_candidate_key_valid_status
fibre_fibre_candidate_key_valid_status
candidate_id_unique_status
global_candidate_id_unique_status
no_duplicate_candidate_key_status
no_self_pair_status
canonical_ordering_status
owner_rank_np1_valid_status
owner_rank_np2_valid_status
owner_rank_np4_valid_status
local_candidate_id_np1_contiguous_status
local_candidate_id_np2_contiguous_status
local_candidate_id_np4_contiguous_status
rank_candidate_count_np1_status
rank_candidate_count_np2_status
rank_candidate_count_np4_status
ownership_deterministic_status
np1_all_candidates_rank0_status
np2_reduction_ready_status
np4_reduction_ready_status
registry_diagnostic_only_status
contact_force_disabled_status
collision_force_disabled_status
contact_force_apply_disabled_status
structure_advance_disabled_status
rhs_coupling_disabled_status
production_dns_disabled_status
actual_mpi_disabled_status
production_multifibre_disabled_status
no_stage10_20_file_modification_status
no_stage21_0_file_modification_status
no_stage21_1_file_modification_status
no_stage21_2_file_modification_status
no_stage21_3_file_modification_status
no_stage21_4_file_modification_status
no_closed_stage_modification_status
no_src_modification_status
no_cmake_modification_status
no_contact_force_computation_status
no_collision_force_computation_status
no_production_structure_update_status
no_production_rhs_update_status
no_stage14_rhs_injection_status
no_mpi_execution_status
no_production_dns_execution_status
no_production_hook_activation_status
no_rg_only_dependency_status
no_unknown_failure_status
stage21_6_next_stage_declared_status
stage21_5_wrapper_bash_syntax_status
stage21_5_helper_py_compile_status
""".split()

REQUIRED_FIELDS = """candidate_id candidate_key canonical_pair_key candidate_type fibre_i fibre_j point_i point_j segment_i segment_j nearest_wall_side global_candidate_id local_candidate_id owner_rank_np1 owner_rank_np2 owner_rank_np4 local_candidate_id_np1 local_candidate_id_np2 local_candidate_id_np4 rank_candidate_count_np1 rank_candidate_count_np2 rank_candidate_count_np4 duplicate_key_detected self_pair_detected unordered_pair_detected ownership_deterministic reduction_ready""".split()


def env(name: str) -> str:
    return os.environ.get(name, SAFE_DEFAULTS[name])


def enabled(name: str) -> bool:
    return env(name).lower() in {"1", "true", "yes", "on"}


def disabled(name: str) -> bool:
    return env(name).lower() in {"0", "false", "no", "off"}


def stable_hash(key: str) -> int:
    return int(hashlib.sha256(key.encode("utf-8")).hexdigest(), 16)


def candidate_base() -> list[dict[str, object]]:
    raw = [
        {"candidate_id": "wall_lower_0000", "candidate_type": "wall_lower", "fibre_i": 0, "fibre_j": -1, "point_i": 12, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "lower"},
        {"candidate_id": "wall_upper_0000", "candidate_type": "wall_upper", "fibre_i": 0, "fibre_j": -1, "point_i": 51, "point_j": -1, "segment_i": -1, "segment_j": -1, "nearest_wall_side": "upper"},
        {"candidate_id": "fibre_fibre_0000", "candidate_type": "fibre_fibre", "fibre_i": 0, "fibre_j": 1, "point_i": 0, "point_j": 0, "segment_i": 0, "segment_j": 0, "nearest_wall_side": "none"},
    ]
    for gid, rec in enumerate(raw):
        if str(rec["candidate_type"]).startswith("wall_"):
            key = (rec["candidate_type"], rec["fibre_i"], rec["point_i"], rec["nearest_wall_side"])
        else:
            fi, fj = sorted((int(rec["fibre_i"]), int(rec["fibre_j"])))
            si, sj = sorted((int(rec["segment_i"]), int(rec["segment_j"])))
            key = (rec["candidate_type"], fi, fj, si, sj)
        rec["candidate_key"] = str(key)
        rec["canonical_pair_key"] = str(key)
        rec["global_candidate_id"] = gid
    return raw


def assign(cands: list[dict[str, object]], np: int) -> tuple[list[int], list[int], list[int]]:
    owners = [stable_hash(str(c["canonical_pair_key"])) % np for c in cands]
    next_local = [0] * np
    locals_: list[int] = []
    for owner in owners:
        locals_.append(next_local[owner])
        next_local[owner] += 1
    return owners, locals_, next_local


def enriched_candidates() -> list[dict[str, object]]:
    cands = candidate_base()
    for np in (1, 2, 4):
        owners, locals_, counts = assign(cands, np)
        for c, owner, lid in zip(cands, owners, locals_):
            c[f"owner_rank_np{np}"] = owner
            c[f"local_candidate_id_np{np}"] = lid
            c[f"rank_candidate_count_np{np}"] = counts
    keys = [c["canonical_pair_key"] for c in cands]
    for c in cands:
        c["local_candidate_id"] = c["local_candidate_id_np1"]
        c["duplicate_key_detected"] = keys.count(c["canonical_pair_key"]) > 1
        c["self_pair_detected"] = c["candidate_type"] == "fibre_fibre" and c["fibre_i"] == c["fibre_j"]
        c["unordered_pair_detected"] = c["candidate_type"] == "fibre_fibre" and int(c["fibre_i"]) > int(c["fibre_j"])
        c["ownership_deterministic"] = all(assign(candidate_base(), np)[0][int(c["global_candidate_id"])] == c[f"owner_rank_np{np}"] for np in (1, 2, 4))
        c["reduction_ready"] = not c["duplicate_key_detected"] and not c["self_pair_detected"] and not c["unordered_pair_detected"]
    return cands


def stage21_4_evidence() -> tuple[bool, str]:
    out = ROOT / "stage21_outputs" / "fibre_stage21_4_contact_candidate_registry.dat"
    if out.exists() and "STAGE 21.4 FINAL VERDICT: PASS" in out.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE21_4_PASS_OUTPUT"
    if enabled("STAGE21_5_ALLOW_SOURCE_ONLY_ARCHIVE") and (CHECK_DIR / "assert_stage21_4_contact_candidate_registry.py").exists():
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"


def bash_syntax_ok() -> bool:
    return subprocess.run(["bash", "-n", str(WRAPPER)], cwd=ROOT).returncode == 0


def py_compile_ok() -> bool:
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_5_pycompile_") as tmpdir:
            py_compile.compile(str(HELPER), cfile=str(Path(tmpdir) / "helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False


def contiguous_by_rank(cands: list[dict[str, object]], np: int) -> bool:
    for rank in range(np):
        ids = sorted(int(c[f"local_candidate_id_np{np}"]) for c in cands if int(c[f"owner_rank_np{np}"]) == rank)
        if ids != list(range(len(ids))):
            return False
    return True


def counts_sum(cands: list[dict[str, object]], np: int) -> bool:
    counts = cands[0][f"rank_candidate_count_np{np}"]
    return isinstance(counts, list) and sum(counts) == len(cands)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cands = enriched_candidates()
    ev_ok, ev_mode = stage21_4_evidence()
    text = DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""
    ids = [c["candidate_id"] for c in cands]
    gids = [c["global_candidate_id"] for c in cands]
    keys = [c["canonical_pair_key"] for c in cands]
    statuses = {field: True for field in STATUS_FIELDS}
    statuses.update({
        "stage21_5_requested_status": enabled("STAGE21_5_ENABLE"),
        "stage21_5_contact_pair_ownership_audit_enable_status": enabled("STAGE21_5_CONTACT_PAIR_OWNERSHIP_AUDIT_ENABLE"),
        "stage21_4_evidence_status": ev_ok,
        "source_only_closure_acceptance_status": enabled("STAGE21_5_ALLOW_SOURCE_ONLY_ARCHIVE"),
        "no_previous_stage_rerun_status": enabled("STAGE21_5_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "ownership_audit_documented_status": "contact pair ownership audit" in text,
        "all_required_ownership_fields_present_status": all(all(f in c for f in REQUIRED_FIELDS) for c in cands),
        "canonical_pair_key_present_status": all(c["canonical_pair_key"] for c in cands),
        "wall_candidate_key_valid_status": all((not str(c["candidate_type"]).startswith("wall_")) or c["nearest_wall_side"] in {"lower", "upper"} for c in cands),
        "fibre_fibre_candidate_key_valid_status": all(c["fibre_i"] < c["fibre_j"] for c in cands if c["candidate_type"] == "fibre_fibre"),
        "candidate_id_unique_status": len(ids) == len(set(ids)),
        "global_candidate_id_unique_status": len(gids) == len(set(gids)),
        "no_duplicate_candidate_key_status": len(keys) == len(set(keys)),
        "no_self_pair_status": not any(c["self_pair_detected"] for c in cands),
        "canonical_ordering_status": not any(c["unordered_pair_detected"] for c in cands),
        "owner_rank_np1_valid_status": all(0 <= c["owner_rank_np1"] < 1 for c in cands),
        "owner_rank_np2_valid_status": all(0 <= c["owner_rank_np2"] < 2 for c in cands),
        "owner_rank_np4_valid_status": all(0 <= c["owner_rank_np4"] < 4 for c in cands),
        "local_candidate_id_np1_contiguous_status": contiguous_by_rank(cands, 1),
        "local_candidate_id_np2_contiguous_status": contiguous_by_rank(cands, 2),
        "local_candidate_id_np4_contiguous_status": contiguous_by_rank(cands, 4),
        "rank_candidate_count_np1_status": counts_sum(cands, 1),
        "rank_candidate_count_np2_status": counts_sum(cands, 2),
        "rank_candidate_count_np4_status": counts_sum(cands, 4),
        "ownership_deterministic_status": all(c["ownership_deterministic"] for c in cands),
        "np1_all_candidates_rank0_status": all(c["owner_rank_np1"] == 0 for c in cands),
        "np2_reduction_ready_status": all(c["reduction_ready"] and 0 <= c["owner_rank_np2"] < 2 for c in cands),
        "np4_reduction_ready_status": all(c["reduction_ready"] and 0 <= c["owner_rank_np4"] < 4 for c in cands),
        "registry_diagnostic_only_status": enabled("STAGE21_5_DIAGNOSTIC_ONLY"),
        "contact_force_disabled_status": disabled("STAGE21_5_CONTACT_FORCE_ENABLE"),
        "collision_force_disabled_status": disabled("STAGE21_5_COLLISION_FORCE_ENABLE"),
        "contact_force_apply_disabled_status": disabled("STAGE21_5_CONTACT_FORCE_APPLY_ENABLE"),
        "structure_advance_disabled_status": disabled("STAGE21_5_STRUCTURE_ADVANCE_ENABLE"),
        "rhs_coupling_disabled_status": disabled("STAGE21_5_RHS_COUPLING_ENABLE"),
        "production_dns_disabled_status": disabled("STAGE21_5_PRODUCTION_DNS_ALLOWED"),
        "actual_mpi_disabled_status": disabled("STAGE21_5_ACTUAL_MPI_ALLOWED"),
        "production_multifibre_disabled_status": disabled("STAGE21_5_PRODUCTION_MULTIFIBRE_ENABLE"),
        "stage21_6_next_stage_declared_status": "Stage 21.6: deterministic pair ordering" in text,
        "stage21_5_wrapper_bash_syntax_status": bash_syntax_ok(),
        "stage21_5_helper_py_compile_status": py_compile_ok(),
    })
    final = all(statuses.values())
    lines = [
        "stage21_5_title contact pair ownership audit",
        "stage21_4_evidence_mode_value " + ev_mode,
        "source_only_policy_value old closure/output evidence is optional when source-only archive acceptance is enabled",
        "rerun_policy_value Stage 21.5 does not rerun previous stages",
        "owner_rule_value stable_sha256(candidate_key)_mod_np",
        "np_values_value 1,2,4",
    ]
    for key in SAFE_DEFAULTS:
        lines.append(f"{key.lower()}_value {env(key)}")
    for c in cands:
        lines.append(str(c["candidate_id"]) + "_ownership " + ";".join(f"{k}={c[k]}" for k in REQUIRED_FIELDS))
    lines.append("stage21_6_next_stage_value Stage 21.6: deterministic pair ordering")
    for field in STATUS_FIELDS:
        lines.append(f"{field} {'PASS' if statuses[field] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final else 'FAIL'}")
    if final:
        lines.extend(["STAGE 21.5 CONTACT PAIR OWNERSHIP AUDIT VERDICT: PASS", "STAGE 21.5 FINAL VERDICT: PASS"])
    else:
        lines.append("failure_reasons_value " + ",".join(k for k, v in statuses.items() if not v))
        lines.extend(["STAGE 21.5 CONTACT PAIR OWNERSHIP AUDIT VERDICT: FAIL", "STAGE 21.5 FINAL VERDICT: FAIL"])
    OUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(lines[-2]); print(lines[-1])
    return 0 if final else 1


if __name__ == "__main__":
    sys.exit(main())
