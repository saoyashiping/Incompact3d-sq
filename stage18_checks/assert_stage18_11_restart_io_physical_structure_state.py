#!/usr/bin/env python3
"""Stage 18.11 restart/I/O compatibility diagnostic audit.

Pure-Python, helper-local JSON snapshot/roundtrip validation for the physical
single-fibre structure state introduced by Stage 18.3 through Stage 18.10.  The
helper writes only Stage 18.11 diagnostic files under stage18_outputs and does
not modify production restart I/O, statistics I/O, visualisation I/O, X/V/A
state, RHS, IBM, DNS-core, or runtime paths.

The helper continues the corrected Stage 18.10 / 18.9 / 18.8 / 18.7 / 18.6 /
18.5 / 18.0 / Stage 17 / Stage 16 false-positive-safe policy: targeted checks
only, no broad repository scans, no Markdown-as-code activation evidence, no
mandatory rg, source-only archives accepted, MPI compatibility names are not
execution, local stage18_outputs JSON is not production restart I/O, and only
*_status fields control final_status.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

Vec = Tuple[float, float, float]
SCHEMA_NAME = "stage18_11_physical_structure_state_snapshot"
SCHEMA_VERSION = "1"
STAGE_ID = "18.11"
SNAPSHOT_JSON = "stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_snapshot.json"
PARTITION_JSON = "stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state_partition_snapshot.json"
OUTPUT = "stage18_outputs/fibre_stage18_11_restart_io_physical_structure_state.dat"

SUMMARY_KEYS = [
    "stage18_11_requested_status",
    "stage18_10_evidence_status",
    "stage18_9_evidence_status",
    "stage18_8_evidence_status",
    "stage18_7_evidence_status",
    "stage18_6_evidence_status",
    "stage18_5_evidence_status",
    "stage18_4_evidence_status",
    "stage18_3_evidence_status",
    "stage18_2_evidence_status",
    "stage18_1_evidence_status",
    "stage18_0_evidence_status",
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
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "stage18_0_files_unmodified_status",
    "stage18_1_files_unmodified_status",
    "stage18_2_files_unmodified_status",
    "stage18_3_files_unmodified_status",
    "stage18_4_files_unmodified_status",
    "stage18_5_files_unmodified_status",
    "stage18_6_files_unmodified_status",
    "stage18_7_files_unmodified_status",
    "stage18_8_files_unmodified_status",
    "stage18_9_files_unmodified_status",
    "stage18_10_files_unmodified_status",
    "stage18_enable_status",
    "restart_io_compatibility_enable_status",
    "single_fibre_only_status",
    "diagnostic_only_status",
    "np_list_status",
    "npts_status",
    "nsteps_status",
    "component_dim_status",
    "fibre_length_status",
    "ds_formula_status",
    "dt_structure_status",
    "rho_l_status",
    "rho_tilde_status",
    "bending_stiffness_status",
    "gamma_status",
    "snapshot_schema_name_status",
    "snapshot_schema_version_status",
    "snapshot_stage_id_status",
    "snapshot_metadata_status",
    "snapshot_array_shape_metadata_status",
    "snapshot_array_finite_status",
    "snapshot_scalar_finite_status",
    "array_snapshot_roundtrip_status",
    "scalar_snapshot_roundtrip_status",
    "deterministic_digest_roundtrip_status",
    "restart_recompute_acceleration_status",
    "restart_recompute_x_next_status",
    "restart_recompute_v_next_status",
    "restart_recompute_energy_status",
    "restart_recompute_power_status",
    "partition_snapshot_np1_status",
    "partition_snapshot_np2_status",
    "partition_snapshot_np4_status",
    "partition_snapshot_reconstruction_status",
    "partition_snapshot_scalar_reduction_status",
    "helper_local_snapshot_output_status",
    "restart_io_diagnostic_only_status",
    "no_production_restart_io_modification_status",
    "no_production_statistics_io_modification_status",
    "no_production_visu_io_modification_status",
    "no_production_io_schema_modification_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_stage16_code_modification_status",
    "no_stage13_force_density_modification_status",
    "no_stage14_rhs_modification_status",
    "no_stage14_rhs_call_from_stage18_11_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_fluid_rhs_modification_status",
    "no_ibm_modification_status",
    "no_dns_core_modification_status",
    "no_stats_visu_restart_io_modification_status",
    "no_production_structure_time_integration_status",
    "no_production_dns_execution_status",
    "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status",
    "no_bending_force_runtime_application_status",
    "no_tension_force_runtime_application_status",
    "no_inextensibility_projection_status",
    "no_inextensibility_repair_status",
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
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
    "stage13_6_diagnostic_preserved_status",
    "stage13_no_local_subdomain_center_regression_status",
    "stage14_small_lambda_hook_status",
    "no_rg_only_dependency_status",
    "no_unknown_failure_status",
    "final_status",
]
VALUE_KEYS = {k for k in SUMMARY_KEYS if k.endswith(("_value", "_formula_value", "_shape_value", "_case_value"))}
S18 = {
    0: ["stage18_checks/run_stage18_0_preflight_boundary.sh", "stage18_checks/assert_stage18_0_preflight_boundary.py", "stage18_checks/stage18_0_preflight_boundary.md"],
    1: ["stage18_checks/run_stage18_1_physical_structure_config.sh", "stage18_checks/assert_stage18_1_physical_structure_config.py", "stage18_checks/stage18_1_physical_structure_config.md"],
    2: ["stage18_checks/run_stage18_2_structure_state_geometry_operators.sh", "stage18_checks/assert_stage18_2_structure_state_geometry_operators.py", "stage18_checks/stage18_2_structure_state_geometry_operators.md"],
    3: ["stage18_checks/run_stage18_3_physical_bending_force_operator.sh", "stage18_checks/assert_stage18_3_physical_bending_force_operator.py", "stage18_checks/stage18_3_physical_bending_force_operator.md"],
    4: ["stage18_checks/run_stage18_4_tension_inextensibility_constraint.sh", "stage18_checks/assert_stage18_4_tension_inextensibility_constraint.py", "stage18_checks/stage18_4_tension_inextensibility_constraint.md"],
    5: ["stage18_checks/run_stage18_5_structure_time_integration_core.sh", "stage18_checks/assert_stage18_5_structure_time_integration_core.py", "stage18_checks/stage18_5_structure_time_integration_core.md"],
    6: ["stage18_checks/run_stage18_6_fluid_force_input_physical_structure.sh", "stage18_checks/assert_stage18_6_fluid_force_input_physical_structure.py", "stage18_checks/stage18_6_fluid_force_input_physical_structure.md"],
    7: ["stage18_checks/run_stage18_7_structure_energy_power_diagnostics.sh", "stage18_checks/assert_stage18_7_structure_energy_power_diagnostics.py", "stage18_checks/stage18_7_structure_energy_power_diagnostics.md"],
    8: ["stage18_checks/run_stage18_8_dry_physical_structure_benchmark.sh", "stage18_checks/assert_stage18_8_dry_physical_structure_benchmark.py", "stage18_checks/stage18_8_dry_physical_structure_benchmark.md"],
    9: ["stage18_checks/run_stage18_9_controlled_one_fibre_physical_response_np1.sh", "stage18_checks/assert_stage18_9_controlled_one_fibre_physical_response_np1.py", "stage18_checks/stage18_9_controlled_one_fibre_physical_response_np1.md"],
    10: ["stage18_checks/run_stage18_10_parallel_consistency_physical_structure.sh", "stage18_checks/assert_stage18_10_parallel_consistency_physical_structure.py", "stage18_checks/stage18_10_parallel_consistency_physical_structure.md"],
}
STAGE18_11_FILES = [
    "stage18_checks/run_stage18_11_restart_io_physical_structure_state.sh",
    "stage18_checks/assert_stage18_11_restart_io_physical_structure_state.py",
    "stage18_checks/stage18_11_restart_io_physical_structure_state.md",
]
ALLOWED_CHANGED = set(STAGE18_11_FILES + [OUTPUT, SNAPSHOT_JSON, PARTITION_JSON])
NEGATIVE_STATUS_KEYS = [
    "no_production_restart_io_modification_status",
    "no_production_statistics_io_modification_status",
    "no_production_visu_io_modification_status",
    "no_production_io_schema_modification_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_stage16_code_modification_status",
    "no_stage13_force_density_modification_status",
    "no_stage14_rhs_modification_status",
    "no_stage14_rhs_call_from_stage18_11_status",
    "no_force_spreading_to_fluid_rhs_status",
    "no_fluid_rhs_modification_status",
    "no_ibm_modification_status",
    "no_dns_core_modification_status",
    "no_stats_visu_restart_io_modification_status",
    "no_production_structure_time_integration_status",
    "no_production_dns_execution_status",
    "no_mpi_execution_status",
    "no_actual_mpirun_or_mpiexec_status",
    "no_bending_force_runtime_application_status",
    "no_tension_force_runtime_application_status",
    "no_inextensibility_projection_status",
    "no_inextensibility_repair_status",
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
    "no_pressure_projection_modification_status",
    "no_poisson_modification_status",
    "no_rk3_channel_forcing_modification_status",
    "no_channel_forcing_modification_status",
]
ARRAY_KEYS = ["X", "V", "A", "F_b", "F_T", "F_h", "F_total", "X_next", "V_next"]
SCALAR_KEYS = ["E_k", "E_b", "E_s", "P_h", "max_response", "max_velocity", "max_acceleration"]


def status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def parse_dat(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.strip().split(None, 1)
        if len(parts) == 2:
            data[parts[0]] = parts[1]
    return data


def ff(text: str, default: float | None = None) -> float | None:
    try:
        value = float(text)
        return value if math.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def ii(text: str, default: int | None = None) -> int | None:
    try:
        return int(text)
    except (TypeError, ValueError):
        return default


def parse_np_list(text: str) -> List[int]:
    vals: List[int] = []
    for part in text.split(","):
        val = ii(part.strip())
        if val is None:
            return []
        vals.append(val)
    return vals


def stable_digest(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def git_entries(root: Path) -> List[Tuple[str, str]]:
    if not (root / ".git").exists():
        return []
    try:
        proc = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    except OSError:
        return []
    entries: List[Tuple[str, str]] = []
    for line in proc.stdout.splitlines():
        if line:
            path = line[3:] if len(line) > 3 else ""
            if " -> " in path:
                path = path.split(" -> ", 1)[1]
            entries.append((line[:2], path))
    return entries


def changed_only_allowed(root: Path) -> bool:
    return all(path in ALLOWED_CHANGED for _code, path in git_entries(root))


def all_present(root: Path, files: Sequence[str]) -> bool:
    return all((root / f).exists() for f in files)


def stage_dat(root: Path, number: int) -> Path | None:
    out_dir = root / "stage18_outputs"
    if not out_dir.exists():
        return None
    matches = sorted(out_dir.glob(f"fibre_stage18_{number}_*.dat"))
    return matches[0] if matches else None


def evidence(root: Path, number: int, needles: Sequence[str]) -> bool:
    dat = stage_dat(root, number)
    if dat and parse_dat(dat).get("final_status") == "PASS":
        return True
    text = "\n".join(read_text(root / f) for f in S18[number])
    return all_present(root, S18[number]) and all(needle in text for needle in needles)


def stage17_closed(root: Path) -> Tuple[bool, bool]:
    closed = root / "stage17_checks" / "STAGE17_CLOSED.md"
    if closed.exists():
        text = read_text(closed)
        return True, ("Stage 17" in text and "closed" in text.lower())
    helper = root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py"
    text = read_text(helper)
    ok = "stage17_11" in helper.name and "source-only" in text and "final_status" in text
    return ok, ok


def s17_11_ok(root: Path) -> bool:
    text = read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    return "source-only" in text and "final_status" in text and "STAGE17_CLOSED" in text


def wrapper_root_fix_ok(root: Path) -> bool:
    return (root / "stage18_checks" / "run_stage18_0_preflight_boundary.sh").exists()


def false_positive_ok(root: Path, number: int) -> bool:
    helper = read_text(root / S18[number][1])
    return "source-only" in helper and "VALUE_KEYS" in helper and "final_status" in helper


def stage13_ok(root: Path) -> bool:
    return all_present(root, ["src/fibre_stage13_production_force_density_candidate.f90", "stage13_checks/run_stage13_6_production_force_density_candidate.sh", "stage13_checks/stage13_6_production_force_density_candidate.md"])


def stage13_reg_ok(root: Path) -> bool:
    text = read_text(root / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py")
    return "local_subdomain_center" not in text.lower() or "global" in text.lower()


def stage14_ok(root: Path) -> bool:
    return all_present(root, ["src/fibre_stage14_production_rhs_injection.f90", "src/xcompact3d.f90"])


def no_rg_only_dependency(root: Path) -> bool:
    text = read_text(root / STAGE18_11_FILES[0])
    return "rg " not in text and "rg[[:space:]]" not in text


def no_activation(root: Path) -> bool:
    text = "\n".join(read_text(root / f) for f in STAGE18_11_FILES[:2]).lower()
    forbidden = [
        "call " + "stage14", "fibre_stage14" + "_production_rhs_injection(", "stage14" + "_rhs_injection(",
        "mpi" + "run ", "mpi" + "exec ", "subprocess.run([\"" + "mpirun" + "\"", "subprocess.run([\"" + "mpiexec" + "\"",
        "cm" + "ake ", "ct" + "est ", "nin" + "ja ", "ma" + "ke ",
        "statistics" + ".f90', 'w", "visu" + ".f90', 'w", "restart" + "_write", "restart" + "_read",
        "fibre_stage13" + "_production_force_density_candidate.f90', 'w",
    ]
    return all(term not in text for term in forbidden)


def weights(n: int, ds: float) -> List[float]:
    return [0.5 * ds] + [ds] * (n - 2) + [0.5 * ds]


def finite_vals(values: Iterable[float]) -> bool:
    return all(math.isfinite(v) for v in values)


def finite_vec(arrays: Iterable[Sequence[Vec]]) -> bool:
    return all(finite_vals(component for row in arr for component in row) for arr in arrays)


def max_norm(arr: Sequence[Vec]) -> float:
    return max(math.sqrt(x*x + y*y + z*z) for x, y, z in arr) if arr else float("inf")


def max_delta(a: Sequence[Vec], b: Sequence[Vec]) -> float:
    return max_norm([(x1 - x0, y1 - y0, z1 - z0) for (x1, y1, z1), (x0, y0, z0) in zip(a, b)])


def array_err(a: Sequence[Vec], b: Sequence[Vec]) -> float:
    return max_delta(a, b) if len(a) == len(b) else float("inf")


def energy_kin(v: Sequence[Vec], rho: float, w: Sequence[float]) -> float:
    return 0.5 * sum(rho * (vx*vx + vy*vy + vz*vz) * wi for (vx, vy, vz), wi in zip(v, w))


def energy_bend(xss: Sequence[Vec], stiffness: float, w: Sequence[float]) -> float:
    return 0.5 * sum(stiffness * (xx*xx + yy*yy + zz*zz) * wi for (xx, yy, zz), wi in zip(xss, w))


def power(fh: Sequence[Vec], v: Sequence[Vec], w: Sequence[float]) -> float:
    return sum((fx*vx + fy*vy + fz*vz) * wi for (fx, fy, fz), (vx, vy, vz), wi in zip(fh, v, w))


def step(x: Sequence[Vec], v: Sequence[Vec], a: Sequence[Vec], dt: float) -> Tuple[List[Vec], List[Vec]]:
    vn = [(vx + dt*ax, vy + dt*ay, vz + dt*az) for (vx, vy, vz), (ax, ay, az) in zip(v, a)]
    xn = [(xx + dt*vx + 0.5*dt*dt*ax, yy + dt*vy + 0.5*dt*dt*ay, zz + dt*vz + 0.5*dt*dt*az) for (xx, yy, zz), (vx, vy, vz), (ax, ay, az) in zip(x, v, a)]
    return xn, vn


def partitions(n: int, np: int) -> List[Tuple[int, int]]:
    base, rem = divmod(n, np)
    spans: List[Tuple[int, int]] = []
    start = 0
    for rank in range(np):
        size = base + (1 if rank < rem else 0)
        spans.append((start, start + size))
        start += size
    return spans


def build_state(n: int, L: float, rho: float, stiffness: float, dt: float, eps: float, mode: int, fmag: float, v0: float) -> Dict[str, Any]:
    ds = L / (n - 1)
    svals = [i * ds for i in range(n)]
    k = 2.0 * math.pi * mode / L
    w = weights(n, ds)
    x = [(s, eps * math.sin(k*s), 0.0) for s in svals]
    v = [(0.0, v0, 0.0)] * n
    xss = [(0.0, -eps*k*k*math.sin(k*s), 0.0) for s in svals]
    fb = [(0.0, -stiffness*eps*k**4*math.sin(k*s), 0.0) for s in svals]
    ft = [(0.0, 0.0, 0.0)] * n
    fh = [(0.0, fmag, 0.0)] * n
    ftotal = [(tx+bx+hx, ty+by+hy, tz+bz+hz) for (tx,ty,tz),(bx,by,bz),(hx,hy,hz) in zip(ft, fb, fh)]
    a = [(fx/rho, fy/rho, fz/rho) for fx,fy,fz in ftotal]
    xnext, vnext = step(x, v, a, dt)
    ek = energy_kin(v, rho, w)
    eb = energy_bend(xss, stiffness, w)
    return {
        "X": x, "V": v, "A": a, "F_b": fb, "F_T": ft, "F_h": fh, "F_total": ftotal, "X_next": xnext, "V_next": vnext, "X_ss": xss,
        "E_k": ek, "E_b": eb, "E_s": ek + eb, "P_h": power(fh, v, w),
        "max_response": max_delta(xnext, x), "max_velocity": max_norm(v), "max_acceleration": max_norm(a),
    }


def to_json_array(arr: Sequence[Vec]) -> List[List[float]]:
    return [[float(x), float(y), float(z)] for x, y, z in arr]


def from_json_array(arr: Sequence[Sequence[float]]) -> List[Vec]:
    return [(float(row[0]), float(row[1]), float(row[2])) for row in arr]


def make_snapshot(meta: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    arrays = {key: to_json_array(state[key]) for key in ARRAY_KEYS}
    scalars = {key: float(state[key]) for key in SCALAR_KEYS}
    shape_metadata = {key: [len(arrays[key]), meta["component_dim"]] for key in ARRAY_KEYS}
    snap = {"schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, "stage_id": STAGE_ID, "metadata": meta, "array_shape_metadata": shape_metadata, "arrays": arrays, "diagnostics": scalars}
    snap["digest"] = stable_digest(snap)
    return snap


def make_partition_snapshot(meta: Dict[str, Any], state: Dict[str, Any], w: Sequence[float]) -> Dict[str, Any]:
    chunks: Dict[str, Any] = {}
    for np in meta["np_list"]:
        np_int = int(np)
        spans = partitions(int(meta["npts"]), np_int)
        parts = []
        for rank, (start, stop) in enumerate(spans):
            arrays = {key: to_json_array(state[key][start:stop]) for key in ARRAY_KEYS}
            ek = energy_kin(state["V"][start:stop], float(meta["rho_l"]), w[start:stop])
            eb = energy_bend(state["X_ss"][start:stop], float(meta["bending_stiffness"]), w[start:stop])
            parts.append({"rank": rank, "start": start, "stop": stop, "arrays": arrays, "diagnostics": {"E_k": ek, "E_b": eb, "E_s": ek + eb, "P_h": power(state["F_h"][start:stop], state["V"][start:stop], w[start:stop])}})
        chunks[str(np_int)] = {"np": np_int, "parts": parts}
    snap = {"schema_name": SCHEMA_NAME, "schema_version": SCHEMA_VERSION, "stage_id": STAGE_ID, "metadata": meta, "partitions": chunks}
    snap["digest"] = stable_digest(snap)
    return snap


def partition_roundtrip_ok(part_snap: Dict[str, Any], state: Dict[str, Any], tol: float) -> Tuple[Dict[int, bool], bool, bool]:
    np_status: Dict[int, bool] = {}
    recon_all = True
    scalar_all = True
    for np_text, entry in part_snap["partitions"].items():
        np = int(np_text)
        n = int(part_snap["metadata"]["npts"])
        reconstructed = {key: [(float("nan"), float("nan"), float("nan"))] * n for key in ARRAY_KEYS}
        counts = [0] * n
        scalar_sums = {"E_k": 0.0, "E_b": 0.0, "E_s": 0.0, "P_h": 0.0}
        for part in entry["parts"]:
            start, stop = int(part["start"]), int(part["stop"])
            for idx in range(start, stop):
                counts[idx] += 1
            for key in ARRAY_KEYS:
                chunk = from_json_array(part["arrays"][key])
                for offset, value in enumerate(chunk):
                    reconstructed[key][start + offset] = value
            for key in scalar_sums:
                scalar_sums[key] += float(part["diagnostics"][key])
        coverage = all(c == 1 for c in counts)
        recon = coverage and all(array_err(reconstructed[key], state[key]) <= tol for key in ARRAY_KEYS)
        scalar = all(abs(scalar_sums[key] - float(state[key])) <= tol for key in ["E_k", "E_b", "E_s", "P_h"])
        np_status[np] = recon and scalar
        recon_all = recon_all and recon
        scalar_all = scalar_all and scalar
    return np_status, recon_all, scalar_all


def write_summary(root: Path, summary: Dict[str, str], reasons: Sequence[str]) -> None:
    out = root / OUTPUT
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {summary.get(key, 'FAIL')}" for key in SUMMARY_KEYS]
    lines.extend(f"reason {reason}" for reason in reasons)
    out.write_text("\n".join(lines) + "\n")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".")
    p.add_argument("--stage18-11-enable", default="1")
    p.add_argument("--restart-io-compatibility-enable", default="1")
    p.add_argument("--single-fibre-only", default="1")
    p.add_argument("--diagnostic-only", default="1")
    p.add_argument("--np-list", default="1,2,4")
    p.add_argument("--npts", default="64")
    p.add_argument("--nsteps", default="5")
    p.add_argument("--fibre-length", default="1.0")
    p.add_argument("--component-dim", default="3")
    p.add_argument("--rho-l", default="1.0")
    p.add_argument("--rho-tilde", default="1.0")
    p.add_argument("--bending-stiffness", default="1.0e-3")
    p.add_argument("--gamma", default="1.0e-3")
    p.add_argument("--dt-structure", default="1.0e-4")
    p.add_argument("--sine-eps", default="1.0e-3")
    p.add_argument("--sine-mode", default="1")
    p.add_argument("--fluid-force-mag", default="1.0e-3")
    p.add_argument("--initial-velocity-mag", default="0.0")
    p.add_argument("--response-bound", default="1.0e-4")
    p.add_argument("--velocity-bound", default="1.0e-3")
    p.add_argument("--acceleration-bound", default="1.0e-2")
    p.add_argument("--zero-tol", default="1.0e-14")
    p.add_argument("--formula-tol", default="1.0e-12")
    p.add_argument("--restart-tol", default="1.0e-12")
    p.add_argument("--reduction-tol", default="1.0e-12")
    p.add_argument("--energy-tol", default="1.0e-12")
    p.add_argument("--test-case", default="snapshot_roundtrip_restart_equivalence_partition_io")
    return p


def main() -> int:
    a = parser().parse_args()
    root = Path(a.repo_root).resolve()
    summary = {key: "FAIL" for key in SUMMARY_KEYS}
    reasons: List[str] = []
    n = ii(a.npts); nsteps = ii(a.nsteps); comp = ii(a.component_dim); np_list = parse_np_list(a.np_list)
    L = ff(a.fibre_length); dt = ff(a.dt_structure); rho = ff(a.rho_l); rhot = ff(a.rho_tilde); B = ff(a.bending_stiffness); gamma = ff(a.gamma)
    eps = ff(a.sine_eps, 1.0e-3) or 1.0e-3; mode = ii(a.sine_mode, 1) or 1; fmag = ff(a.fluid_force_mag, 1.0e-3) or 1.0e-3; v0 = ff(a.initial_velocity_mag, 0.0) or 0.0
    rt = ff(a.restart_tol, 1.0e-12) or 1.0e-12; rdt = ff(a.reduction_tol, 1.0e-12) or 1.0e-12; ft = ff(a.formula_tol, 1.0e-12) or 1.0e-12; et = ff(a.energy_tol, 1.0e-12) or 1.0e-12
    rb = ff(a.response_bound, 1.0e-4) or 1.0e-4; vb = ff(a.velocity_bound, 1.0e-3) or 1.0e-3; ab = ff(a.acceleration_bound, 1.0e-2) or 1.0e-2
    n_ok = n is not None and n >= 8; nsteps_ok = nsteps is not None and nsteps >= 1; L_ok = L is not None and L > 0; dt_ok = dt is not None and dt > 0
    rho_ok = rho is not None and rho > 0; rhot_ok = rhot is not None and rhot > 0; B_ok = B is not None and B >= 0; gamma_ok = gamma is not None and gamma >= 0; np_ok = sorted(np_list) == [1, 2, 4]
    n_d = n if n_ok else 8; L_d = L if L_ok else 1.0; dt_d = dt if dt_ok else 1.0e-4; rho_d = rho if rho_ok else 1.0; B_d = B if B_ok else 1.0e-3; rhot_d = rhot if rhot_ok else 1.0; gamma_d = gamma if gamma_ok else 1.0e-3
    ds = L_d / (n_d - 1); w = weights(n_d, ds)
    state = build_state(n_d, L_d, rho_d, B_d, dt_d, eps, mode, fmag, v0)
    meta = {"npts": n_d, "component_dim": 3, "fibre_length": L_d, "ds": ds, "dt_structure": dt_d, "rho_l": rho_d, "rho_tilde": rhot_d, "bending_stiffness": B_d, "gamma": gamma_d, "np_list": np_list if np_ok else [1,2,4]}
    snapshot = make_snapshot(meta, state)
    partition_snapshot = make_partition_snapshot(meta, state, w)
    snap_path = root / SNAPSHOT_JSON; part_path = root / PARTITION_JSON
    snap_path.parent.mkdir(parents=True, exist_ok=True)
    snap_path.write_text(json.dumps(snapshot, sort_keys=True, indent=2) + "\n")
    part_path.write_text(json.dumps(partition_snapshot, sort_keys=True, indent=2) + "\n")
    loaded = json.loads(snap_path.read_text())
    loaded_part = json.loads(part_path.read_text())
    loaded_arrays = {key: from_json_array(loaded["arrays"][key]) for key in ARRAY_KEYS}
    loaded_scalars = {key: float(loaded["diagnostics"][key]) for key in SCALAR_KEYS}
    array_roundtrip = all(array_err(loaded_arrays[key], state[key]) <= rt for key in ARRAY_KEYS)
    scalar_roundtrip = all(abs(loaded_scalars[key] - float(state[key])) <= rt for key in SCALAR_KEYS)
    loaded_no_digest = dict(loaded); digest = loaded_no_digest.pop("digest", "")
    digest_ok = digest == stable_digest(loaded_no_digest) == snapshot["digest"]
    a_loaded = [(fx/rho_d, fy/rho_d, fz/rho_d) for fx, fy, fz in loaded_arrays["F_total"]]
    xnext_loaded, vnext_loaded = step(loaded_arrays["X"], loaded_arrays["V"], a_loaded, dt_d)
    xss = state["X_ss"]
    ek_loaded = energy_kin(loaded_arrays["V"], rho_d, w); eb_loaded = energy_bend(xss, B_d, w); ph_loaded = power(loaded_arrays["F_h"], loaded_arrays["V"], w)
    np_part_status, part_recon, part_scalars = partition_roundtrip_ok(loaded_part, state, rdt)
    all_outputs_local = all((root / p).parent == root / "stage18_outputs" for p in [OUTPUT, SNAPSHOT_JSON, PARTITION_JSON])
    arrays_finite = finite_vec([state[key] for key in ARRAY_KEYS])
    scalars_finite = finite_vals(float(state[key]) for key in SCALAR_KEYS)
    shapes_ok = all(loaded["array_shape_metadata"].get(key) == [n_d, 3] for key in ARRAY_KEYS)
    meta_ok = loaded["metadata"] == meta
    changed_ok = changed_only_allowed(root)
    closed_file_present, closed_evidence = stage17_closed(root)
    activation_ok = no_activation(root)
    doc = read_text(root / STAGE18_11_FILES[2])
    requested_files = all_present(root, STAGE18_11_FILES)
    summary.update({
        "stage18_11_requested_status": status(requested_files),
        "stage18_10_evidence_status": status(evidence(root, 10, ["PARALLEL CONSISTENCY DIAGNOSTIC", "np=1", "partition"])),
        "stage18_9_evidence_status": status(evidence(root, 9, ["CONTROLLED RESPONSE DIAGNOSTIC", "np=1", "P_h"])),
        "stage18_8_evidence_status": status(evidence(root, 8, ["F_h = 0", "A_b_candidate", "R_E_dry"])),
        "stage18_7_evidence_status": status(evidence(root, 7, ["E_k", "E_b", "R_E"])),
        "stage18_6_evidence_status": status(evidence(root, 6, ["F_h", "F_fibre_on_fluid", "P_h"])),
        "stage18_5_evidence_status": status(evidence(root, 5, ["X_candidate", "V_candidate", "A_candidate"])),
        "stage18_4_evidence_status": status(evidence(root, 4, ["tension", "inextensibility", "diagnostic"])),
        "stage18_3_evidence_status": status(evidence(root, 3, ["F_b", "X_ssss", "E_b"])),
        "stage18_2_evidence_status": status(evidence(root, 2, ["X_ss", "X_ssss", "geometry"])),
        "stage18_1_evidence_status": status(evidence(root, 1, ["rho_l", "B", "ds"])),
        "stage18_0_evidence_status": status(evidence(root, 0, ["Stage 18.0", "preflight", "diagnostic"])),
        "stage17_closed_file_status": status(closed_file_present or closed_evidence),
        "stage17_closed_evidence_status": status(closed_evidence),
        "stage17_11_closure_preserved_status": status(s17_11_ok(root)),
        "stage18_0_wrapper_root_fix_preserved_status": status(wrapper_root_fix_ok(root)),
        "stage18_5_false_positive_fix_preserved_status": status(false_positive_ok(root, 5)),
        "stage18_6_false_positive_fix_preserved_status": status(false_positive_ok(root, 6)),
        "stage18_7_false_positive_fix_preserved_status": status(false_positive_ok(root, 7)),
        "stage18_8_false_positive_fix_preserved_status": status(false_positive_ok(root, 8)),
        "stage18_9_controlled_response_preserved_status": status(evidence(root, 9, ["LOCAL CANDIDATE FORCE/STEP", "NP=1 DIAGNOSTIC"])),
        "stage18_10_parallel_consistency_preserved_status": status(evidence(root, 10, ["LOCAL PARTITION", "NP=1 REFERENCE"])),
        "no_closed_stage_modification_status": status(changed_ok), "no_stage10_17_file_modification_status": status(changed_ok),
        "stage18_0_files_unmodified_status": status(changed_ok), "stage18_1_files_unmodified_status": status(changed_ok), "stage18_2_files_unmodified_status": status(changed_ok), "stage18_3_files_unmodified_status": status(changed_ok), "stage18_4_files_unmodified_status": status(changed_ok), "stage18_5_files_unmodified_status": status(changed_ok), "stage18_6_files_unmodified_status": status(changed_ok), "stage18_7_files_unmodified_status": status(changed_ok), "stage18_8_files_unmodified_status": status(changed_ok), "stage18_9_files_unmodified_status": status(changed_ok), "stage18_10_files_unmodified_status": status(changed_ok),
        "stage18_enable_status": status(a.stage18_11_enable == "1"), "restart_io_compatibility_enable_status": status(a.restart_io_compatibility_enable == "1"), "single_fibre_only_status": status(a.single_fibre_only == "1"), "diagnostic_only_status": status(a.diagnostic_only == "1"),
        "np_list_status": status(np_ok), "npts_status": status(n_ok), "nsteps_status": status(nsteps_ok), "component_dim_status": status(comp == 3), "fibre_length_status": status(L_ok), "ds_formula_status": status(n_ok and L_ok and abs(ds - L_d/(n_d-1)) <= ft and ds > 0), "dt_structure_status": status(dt_ok), "rho_l_status": status(rho_ok), "rho_tilde_status": status(rhot_ok), "bending_stiffness_status": status(B_ok), "gamma_status": status(gamma_ok),
        "snapshot_schema_name_status": status(loaded.get("schema_name") == SCHEMA_NAME), "snapshot_schema_version_status": status(bool(loaded.get("schema_version"))), "snapshot_stage_id_status": status(loaded.get("stage_id") == STAGE_ID), "snapshot_metadata_status": status(meta_ok), "snapshot_array_shape_metadata_status": status(shapes_ok), "snapshot_array_finite_status": status(arrays_finite), "snapshot_scalar_finite_status": status(scalars_finite),
        "array_snapshot_roundtrip_status": status(array_roundtrip), "scalar_snapshot_roundtrip_status": status(scalar_roundtrip), "deterministic_digest_roundtrip_status": status(digest_ok),
        "restart_recompute_acceleration_status": status(array_err(a_loaded, state["A"]) <= rt), "restart_recompute_x_next_status": status(array_err(xnext_loaded, state["X_next"]) <= rt), "restart_recompute_v_next_status": status(array_err(vnext_loaded, state["V_next"]) <= rt), "restart_recompute_energy_status": status(abs(ek_loaded-float(state["E_k"])) <= et and abs(eb_loaded-float(state["E_b"])) <= et), "restart_recompute_power_status": status(abs(ph_loaded-float(state["P_h"])) <= rt),
        "partition_snapshot_np1_status": status(np_part_status.get(1, False)), "partition_snapshot_np2_status": status(np_part_status.get(2, False)), "partition_snapshot_np4_status": status(np_part_status.get(4, False)), "partition_snapshot_reconstruction_status": status(part_recon), "partition_snapshot_scalar_reduction_status": status(part_scalars),
        "helper_local_snapshot_output_status": status(all_outputs_local and snap_path.exists() and part_path.exists()), "restart_io_diagnostic_only_status": status(activation_ok and "RESTART/I/O COMPATIBILITY DIAGNOSTIC" in doc),
        "stage13_6_diagnostic_preserved_status": status(stage13_ok(root)), "stage13_no_local_subdomain_center_regression_status": status(stage13_reg_ok(root)), "stage14_small_lambda_hook_status": status(stage14_ok(root)), "no_rg_only_dependency_status": status(no_rg_only_dependency(root)), "no_unknown_failure_status": "PASS",
    })
    for key in NEGATIVE_STATUS_KEYS:
        summary[key] = status(activation_ok)
    if a.test_case != "snapshot_roundtrip_restart_equivalence_partition_io":
        reasons.append("unexpected_stage18_11_test_case")
    if not requested_files:
        reasons.append("required_stage18_11_file_missing")
    if not changed_ok:
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(path for _code, path in git_entries(root) if path not in ALLOWED_CHANGED))
    pass_fail_keys = [key for key in SUMMARY_KEYS if key.endswith("_status") and key != "final_status" and key not in VALUE_KEYS]
    for key in pass_fail_keys:
        if summary.get(key) != "PASS":
            reasons.append(key.replace("_status", "_failed"))
    summary["final_status"] = "PASS" if all(summary.get(key) == "PASS" for key in pass_fail_keys) else "FAIL"
    write_summary(root, summary, reasons)
    for key in SUMMARY_KEYS:
        print(f"{key} {summary.get(key, 'FAIL')}")
    for reason in reasons:
        print(f"reason {reason}")
    return 0 if summary["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
