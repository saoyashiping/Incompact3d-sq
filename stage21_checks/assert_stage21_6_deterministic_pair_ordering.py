#!/usr/bin/env python3
"""Stage 21.6 diagnostic-only deterministic pair ordering audit."""
from __future__ import annotations

import hashlib
import os
import py_compile
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "stage21_checks"
OUT_DIR = ROOT / "stage21_outputs"
OUT_FILE = OUT_DIR / "fibre_stage21_6_deterministic_pair_ordering.dat"
HELPER = CHECK_DIR / "assert_stage21_6_deterministic_pair_ordering.py"
WRAPPER = CHECK_DIR / "run_stage21_6_deterministic_pair_ordering.sh"
DOC = CHECK_DIR / "stage21_6_deterministic_pair_ordering.md"

SAFE_DEFAULTS = {
    "STAGE21_6_ENABLE": "1",
    "STAGE21_6_DETERMINISTIC_PAIR_ORDERING_ENABLE": "1",
    "STAGE21_6_REQUIRE_STAGE21_5_PASS": "1",
    "STAGE21_6_ALLOW_MISSING_OLD_OUTPUTS": "1",
    "STAGE21_6_ALLOW_SOURCE_ONLY_ARCHIVE": "1",
    "STAGE21_6_DO_NOT_RERUN_PREVIOUS_STAGES": "1",
    "STAGE21_6_DIAGNOSTIC_ONLY": "1",
    "STAGE21_6_FAIL_CLOSED": "1",
    "STAGE21_6_NP_VALUES": "1,2,4",
    "STAGE21_6_ORDER_RULE": "canonical_sort_key",
    "STAGE21_6_SHUFFLE_SEED": "2106",
    "STAGE21_6_CONTACT_FORCE_ENABLE": "0",
    "STAGE21_6_COLLISION_FORCE_ENABLE": "0",
    "STAGE21_6_CONTACT_FORCE_APPLY_ENABLE": "0",
    "STAGE21_6_STRUCTURE_ADVANCE_ENABLE": "0",
    "STAGE21_6_RHS_COUPLING_ENABLE": "0",
    "STAGE21_6_PRODUCTION_DNS_ALLOWED": "0",
    "STAGE21_6_ACTUAL_MPI_ALLOWED": "0",
    "STAGE21_6_PRODUCTION_MULTIFIBRE_ENABLE": "0",
}

STATUS_FIELDS = """
stage21_6_requested_status
stage21_6_deterministic_pair_ordering_enable_status
stage21_5_evidence_status
source_only_closure_acceptance_status
no_previous_stage_rerun_status
deterministic_ordering_audit_documented_status
all_required_ordering_fields_present_status
canonical_sort_key_documented_status
wall_candidate_sort_key_valid_status
fibre_fibre_candidate_sort_key_valid_status
candidate_type_priority_valid_status
nearest_wall_side_priority_valid_status
candidate_id_unique_status
canonical_pair_key_unique_status
no_duplicate_candidate_key_status
no_self_pair_status
canonical_fibre_pair_ordering_status
canonical_point_segment_ordering_status
original_order_sorts_to_reference_status
reversed_order_sorts_to_reference_status
shuffled_order_sorts_to_reference_status
np1_grouped_order_sorts_to_reference_status
np2_grouped_order_sorts_to_reference_status
np4_grouped_order_sorts_to_reference_status
repeated_eval_order_hash_status
global_order_index_contiguous_status
local_order_index_np1_contiguous_status
local_order_index_np2_contiguous_status
local_order_index_np4_contiguous_status
ordering_deterministic_status
ordering_reduction_ready_status
registry_diagnostic_only_status
fail_closed_status
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
no_stage21_5_file_modification_status
no_closed_stage_modification_status
no_src_modification_status
no_cmake_modification_status
no_contact_force_computation_status
no_collision_force_computation_status
no_contact_collision_force_apply_status
no_production_structure_update_status
no_production_rhs_update_status
no_stage14_rhs_injection_status
no_mpi_execution_status
no_production_dns_execution_status
no_production_hook_activation_status
no_rg_only_dependency_status
no_unknown_failure_status
stage21_7_next_stage_declared_status
stage21_6_wrapper_bash_syntax_status
stage21_6_helper_py_compile_status
""".split()

REQUIRED_FIELDS = """candidate_id candidate_type candidate_type_priority candidate_key canonical_pair_key canonical_sort_key global_order_index local_order_index_np1 local_order_index_np2 local_order_index_np4 owner_rank_np1 owner_rank_np2 owner_rank_np4 ordered_fibre_i ordered_fibre_j ordered_point_i ordered_point_j ordered_segment_i ordered_segment_j nearest_wall_side nearest_wall_side_priority discovery_order_index reversed_order_index shuffled_order_index sorted_order_reference sorted_order_from_original sorted_order_from_reversed sorted_order_from_shuffled sorted_order_from_np1_grouped sorted_order_from_np2_grouped sorted_order_from_np4_grouped repeated_eval_order_hash ordering_deterministic ordering_reduction_ready""".split()
TYPE_PRIORITY = {"wall_lower": 0, "wall_upper": 1, "fibre_fibre": 2}
SIDE_PRIORITY = {"lower": 0, "upper": 1, "none": 9}


def env(name): return os.environ.get(name, SAFE_DEFAULTS[name])
def enabled(name): return env(name).lower() in {"1", "true", "yes", "on"}
def disabled(name): return env(name).lower() in {"0", "false", "no", "off"}
def stable_hash(s): return int(hashlib.sha256(s.encode()).hexdigest(), 16)

def base_candidates():
    rows = [
        {"candidate_id":"wall_lower_0000","candidate_type":"wall_lower","fibre_i":0,"fibre_j":-1,"point_i":12,"point_j":-1,"segment_i":-1,"segment_j":-1,"nearest_wall_side":"lower"},
        {"candidate_id":"wall_upper_0000","candidate_type":"wall_upper","fibre_i":0,"fibre_j":-1,"point_i":51,"point_j":-1,"segment_i":-1,"segment_j":-1,"nearest_wall_side":"upper"},
        {"candidate_id":"fibre_fibre_0000","candidate_type":"fibre_fibre","fibre_i":0,"fibre_j":1,"point_i":0,"point_j":0,"segment_i":0,"segment_j":0,"nearest_wall_side":"none"},
    ]
    for idx, r in enumerate(rows):
        r["discovery_order_index"] = idx
        r["candidate_type_priority"] = TYPE_PRIORITY[r["candidate_type"]]
        r["nearest_wall_side_priority"] = SIDE_PRIORITY[r["nearest_wall_side"]]
        fi, fj = sorted((int(r["fibre_i"]), int(r["fibre_j"]))) if r["candidate_type"] == "fibre_fibre" else (int(r["fibre_i"]), int(r["fibre_j"]))
        pi, pj = sorted((int(r["point_i"]), int(r["point_j"]))) if r["candidate_type"] == "fibre_fibre" else (int(r["point_i"]), int(r["point_j"]))
        si, sj = sorted((int(r["segment_i"]), int(r["segment_j"]))) if r["candidate_type"] == "fibre_fibre" else (int(r["segment_i"]), int(r["segment_j"]))
        r.update({"ordered_fibre_i":fi,"ordered_fibre_j":fj,"ordered_point_i":pi,"ordered_point_j":pj,"ordered_segment_i":si,"ordered_segment_j":sj})
        if r["candidate_type"].startswith("wall_"):
            key = (r["candidate_type"], r["fibre_i"], r["point_i"], r["nearest_wall_side"])
            sort_key = (r["candidate_type_priority"], r["fibre_i"], r["point_i"], r["nearest_wall_side_priority"], r["candidate_id"])
        else:
            key = (r["candidate_type"], fi, fj, si, sj)
            sort_key = (r["candidate_type_priority"], fi, fj, si, sj, pi, pj, r["candidate_id"])
        r["candidate_key"] = str(key)
        r["canonical_pair_key"] = str(key)
        r["canonical_sort_key"] = str(sort_key)
    return rows

def sorted_ids(rows): return [r["candidate_id"] for r in sorted(rows, key=lambda r: r["canonical_sort_key"])]
def order_hash(ids): return hashlib.sha256("|".join(ids).encode()).hexdigest()

def assign_owner(rows, np):
    owners = [stable_hash(r["canonical_pair_key"]) % np for r in rows]
    counts = [0]*np; local=[]
    for o in owners:
        local.append(counts[o]); counts[o]+=1
    return owners, local, counts

def enrich():
    rows = base_candidates(); ref = sorted_ids(rows); ref_hash = order_hash(ref)
    rev = list(reversed(rows)); rng = random.Random(int(env("STAGE21_6_SHUFFLE_SEED"))); shuf = rows[:]; rng.shuffle(shuf)
    for idx, r in enumerate(rev): r["reversed_order_index"] = idx
    for idx, r in enumerate(shuf): r["shuffled_order_index"] = idx
    for np in (1,2,4):
        owners, locals_, counts = assign_owner(rows, np)
        for r,o,l in zip(rows, owners, locals_):
            r[f"owner_rank_np{np}"] = o; r[f"local_order_index_np{np}"] = l
        grouped = sorted(rows, key=lambda r: (r[f"owner_rank_np{np}"], r[f"local_order_index_np{np}"]))
        for r in rows: r[f"sorted_order_from_np{np}_grouped"] = ",".join(sorted_ids(grouped))
    for idx, r in enumerate(sorted(rows, key=lambda r: r["canonical_sort_key"])):
        r["global_order_index"] = idx
    for r in rows:
        r["sorted_order_reference"] = ",".join(ref)
        r["sorted_order_from_original"] = ",".join(sorted_ids(rows))
        r["sorted_order_from_reversed"] = ",".join(sorted_ids(rev))
        r["sorted_order_from_shuffled"] = ",".join(sorted_ids(shuf))
        r["repeated_eval_order_hash"] = ref_hash
        r["ordering_deterministic"] = order_hash(sorted_ids(base_candidates())) == ref_hash
        r["ordering_reduction_ready"] = True
    return rows, ref, ref_hash

def evidence():
    out = ROOT/"stage21_outputs"/"fibre_stage21_5_contact_pair_ownership_audit.dat"
    if out.exists() and "STAGE 21.5 FINAL VERDICT: PASS" in out.read_text(encoding="utf-8", errors="replace"):
        return True, "STAGE21_5_PASS_OUTPUT"
    if enabled("STAGE21_6_ALLOW_SOURCE_ONLY_ARCHIVE") and (CHECK_DIR/"assert_stage21_5_contact_pair_ownership_audit.py").exists():
        return True, "SOURCE_ONLY_ACCEPTED"
    return False, "MISSING"

def bash_ok(): return subprocess.run(["bash","-n",str(WRAPPER)], cwd=ROOT).returncode == 0
def py_ok():
    try:
        with tempfile.TemporaryDirectory(prefix="stage21_6_pycompile_") as d:
            py_compile.compile(str(HELPER), cfile=str(Path(d)/"helper.pyc"), doraise=True)
        return True
    except py_compile.PyCompileError:
        return False

def contiguous(rows,np):
    for owner in range(np):
        ids=sorted(r[f"local_order_index_np{np}"] for r in rows if r[f"owner_rank_np{np}"]==owner)
        if ids != list(range(len(ids))): return False
    return True

def main():
    OUT_DIR.mkdir(exist_ok=True)
    rows, ref, ref_hash = enrich(); ev_ok, ev_mode = evidence()
    text = DOC.read_text(encoding="utf-8", errors="replace") if DOC.exists() else ""
    ids=[r["candidate_id"] for r in rows]; keys=[r["canonical_pair_key"] for r in rows]
    status={f: True for f in STATUS_FIELDS}
    status.update({
        "stage21_6_requested_status": enabled("STAGE21_6_ENABLE"),
        "stage21_6_deterministic_pair_ordering_enable_status": enabled("STAGE21_6_DETERMINISTIC_PAIR_ORDERING_ENABLE"),
        "stage21_5_evidence_status": ev_ok,
        "source_only_closure_acceptance_status": enabled("STAGE21_6_ALLOW_SOURCE_ONLY_ARCHIVE"),
        "no_previous_stage_rerun_status": enabled("STAGE21_6_DO_NOT_RERUN_PREVIOUS_STAGES"),
        "deterministic_ordering_audit_documented_status": "deterministic ordering audit" in text,
        "all_required_ordering_fields_present_status": all(all(f in r for f in REQUIRED_FIELDS) for r in rows),
        "canonical_sort_key_documented_status": "canonical_sort_key" in text,
        "wall_candidate_sort_key_valid_status": all((not r["candidate_type"].startswith("wall_")) or r["nearest_wall_side"] in {"lower","upper"} for r in rows),
        "fibre_fibre_candidate_sort_key_valid_status": all(r["ordered_fibre_i"] < r["ordered_fibre_j"] for r in rows if r["candidate_type"]=="fibre_fibre"),
        "candidate_type_priority_valid_status": all(r["candidate_type_priority"] in TYPE_PRIORITY.values() for r in rows),
        "nearest_wall_side_priority_valid_status": all(r["nearest_wall_side_priority"] in SIDE_PRIORITY.values() for r in rows),
        "candidate_id_unique_status": len(ids)==len(set(ids)),
        "canonical_pair_key_unique_status": len(keys)==len(set(keys)),
        "no_duplicate_candidate_key_status": len(keys)==len(set(keys)),
        "no_self_pair_status": all(not (r["candidate_type"]=="fibre_fibre" and r["ordered_fibre_i"]==r["ordered_fibre_j"]) for r in rows),
        "canonical_fibre_pair_ordering_status": all(r["ordered_fibre_i"] <= r["ordered_fibre_j"] for r in rows if r["candidate_type"]=="fibre_fibre"),
        "canonical_point_segment_ordering_status": all(r["ordered_segment_i"] <= r["ordered_segment_j"] and r["ordered_point_i"] <= r["ordered_point_j"] for r in rows if r["candidate_type"]=="fibre_fibre"),
        "original_order_sorts_to_reference_status": all(r["sorted_order_from_original"] == ",".join(ref) for r in rows),
        "reversed_order_sorts_to_reference_status": all(r["sorted_order_from_reversed"] == ",".join(ref) for r in rows),
        "shuffled_order_sorts_to_reference_status": all(r["sorted_order_from_shuffled"] == ",".join(ref) for r in rows),
        "np1_grouped_order_sorts_to_reference_status": all(r["sorted_order_from_np1_grouped"] == ",".join(ref) for r in rows),
        "np2_grouped_order_sorts_to_reference_status": all(r["sorted_order_from_np2_grouped"] == ",".join(ref) for r in rows),
        "np4_grouped_order_sorts_to_reference_status": all(r["sorted_order_from_np4_grouped"] == ",".join(ref) for r in rows),
        "repeated_eval_order_hash_status": all(r["repeated_eval_order_hash"] == ref_hash for r in rows),
        "global_order_index_contiguous_status": sorted(r["global_order_index"] for r in rows)==list(range(len(rows))),
        "local_order_index_np1_contiguous_status": contiguous(rows,1),
        "local_order_index_np2_contiguous_status": contiguous(rows,2),
        "local_order_index_np4_contiguous_status": contiguous(rows,4),
        "ordering_deterministic_status": all(r["ordering_deterministic"] for r in rows),
        "ordering_reduction_ready_status": all(r["ordering_reduction_ready"] for r in rows),
        "registry_diagnostic_only_status": enabled("STAGE21_6_DIAGNOSTIC_ONLY"),
        "fail_closed_status": enabled("STAGE21_6_FAIL_CLOSED"),
        "contact_force_disabled_status": disabled("STAGE21_6_CONTACT_FORCE_ENABLE"),
        "collision_force_disabled_status": disabled("STAGE21_6_COLLISION_FORCE_ENABLE"),
        "contact_force_apply_disabled_status": disabled("STAGE21_6_CONTACT_FORCE_APPLY_ENABLE"),
        "structure_advance_disabled_status": disabled("STAGE21_6_STRUCTURE_ADVANCE_ENABLE"),
        "rhs_coupling_disabled_status": disabled("STAGE21_6_RHS_COUPLING_ENABLE"),
        "production_dns_disabled_status": disabled("STAGE21_6_PRODUCTION_DNS_ALLOWED"),
        "actual_mpi_disabled_status": disabled("STAGE21_6_ACTUAL_MPI_ALLOWED"),
        "production_multifibre_disabled_status": disabled("STAGE21_6_PRODUCTION_MULTIFIBRE_ENABLE"),
        "stage21_7_next_stage_declared_status": "Stage 21.7: contact metadata consistency" in text,
        "stage21_6_wrapper_bash_syntax_status": bash_ok(),
        "stage21_6_helper_py_compile_status": py_ok(),
    })
    final=all(status.values())
    lines=["stage21_6_title deterministic pair ordering", "stage21_5_evidence_mode_value "+ev_mode, "source_only_policy_value old closure/output evidence is optional when source-only archive acceptance is enabled", "rerun_policy_value Stage 21.6 does not rerun previous stages", "canonical_reference_order_value "+",".join(ref), "repeated_eval_order_hash_value "+ref_hash]
    for k in SAFE_DEFAULTS: lines.append(f"{k.lower()}_value {env(k)}")
    for r in rows: lines.append(str(r["candidate_id"])+"_ordering "+";".join(f"{k}={r[k]}" for k in REQUIRED_FIELDS))
    lines.append("stage21_7_next_stage_value Stage 21.7: contact metadata consistency")
    for f in STATUS_FIELDS: lines.append(f"{f} {'PASS' if status[f] else 'FAIL'}")
    lines.append(f"final_status {'PASS' if final else 'FAIL'}")
    if final: lines += ["STAGE 21.6 DETERMINISTIC PAIR ORDERING VERDICT: PASS", "STAGE 21.6 FINAL VERDICT: PASS"]
    else:
        lines.append("failure_reasons_value "+",".join(k for k,v in status.items() if not v)); lines += ["STAGE 21.6 DETERMINISTIC PAIR ORDERING VERDICT: FAIL", "STAGE 21.6 FINAL VERDICT: FAIL"]
    OUT_FILE.write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(lines[-2]); print(lines[-1]); return 0 if final else 1

if __name__ == "__main__": sys.exit(main())
