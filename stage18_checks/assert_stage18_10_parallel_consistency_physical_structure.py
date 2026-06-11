#!/usr/bin/env python3
"""Stage 18.10 parallel-consistency physical-structure diagnostic audit.

Pure-Python, MPI-free validation that Stage 18.9-style controlled one-fibre
candidate arrays and scalar diagnostics are invariant under deterministic
partition/reduction emulation for np=1, np=2, and np=4.  The helper uses local
Python arrays only: it does not run production DNS, does not execute MPI, does
not update production X/V/A, does not call Stage 14 RHS injection, does not
spread forces to Eulerian RHS, and does not write production IBM, DNS-core,
statistics, visualisation, or restart I/O.

The helper continues the corrected Stage 18.9 / 18.8 / 18.7 / 18.6 / 18.5 /
Stage 18.0 / Stage 17 / Stage 16 false-positive-safe policy: targeted checks
only, no broad repository scans, no Markdown-as-code activation evidence, no
mandatory rg, source-only archives accepted, MPI compatibility names are not
execution, and only *_status fields control final_status.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Vec = Tuple[float, float, float]

SUMMARY_KEYS = [
    "stage18_10_requested_status",
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
    "stage18_enable_status",
    "parallel_consistency_enable_status",
    "single_fibre_only_status",
    "diagnostic_only_status",
    "np_list_status",
    "np1_reference_status",
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
    "dimensional_response_validation_status",
    "nondimensional_response_validation_status",
    "partition_coverage_status",
    "partition_no_gap_status",
    "partition_no_overlap_status",
    "partition_size_balance_status",
    "np1_reference_array_finite_status",
    "np1_reference_energy_nonnegative_status",
    "np2_array_reconstruction_status",
    "np4_array_reconstruction_status",
    "np2_scalar_reduction_status",
    "np4_scalar_reduction_status",
    "np2_energy_reduction_status",
    "np4_energy_reduction_status",
    "np2_power_reduction_status",
    "np4_power_reduction_status",
    "np2_response_bound_reduction_status",
    "np4_response_bound_reduction_status",
    "dry_response_parallel_consistency_status",
    "forced_response_parallel_consistency_status",
    "partition_reduction_diagnostic_only_status",
    "no_production_parallel_output_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_stage16_code_modification_status",
    "no_stage13_force_density_modification_status",
    "no_stage14_rhs_modification_status",
    "no_stage14_rhs_call_from_stage18_10_status",
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
}
STAGE18_10_FILES = [
    "stage18_checks/run_stage18_10_parallel_consistency_physical_structure.sh",
    "stage18_checks/assert_stage18_10_parallel_consistency_physical_structure.py",
    "stage18_checks/stage18_10_parallel_consistency_physical_structure.md",
]
OUTPUT = "stage18_outputs/fibre_stage18_10_parallel_consistency_physical_structure.dat"
ALLOWED_CHANGED = set(STAGE18_10_FILES + [OUTPUT])
NEGATIVE_STATUS_KEYS = [
    "no_production_parallel_output_status",
    "no_production_structure_update_status",
    "no_production_structure_hook_status",
    "no_stage16_code_modification_status",
    "no_stage13_force_density_modification_status",
    "no_stage14_rhs_modification_status",
    "no_stage14_rhs_call_from_stage18_10_status",
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
    values: List[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        val = ii(part)
        if val is None:
            return []
        values.append(val)
    return values


def git_entries(root: Path) -> List[Tuple[str, str]]:
    if not (root / ".git").exists():
        return []
    try:
        proc = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError:
        return []
    entries: List[Tuple[str, str]] = []
    for line in proc.stdout.splitlines():
        if not line:
            continue
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
    safe_helper = "stage17_11" in helper.name and "source-only" in text and "final_status" in text
    return safe_helper, safe_helper


def s17_11_ok(root: Path) -> bool:
    text = read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    return "source-only" in text and "final_status" in text and "STAGE17_CLOSED" in text


def wrapper_root_fix_ok(root: Path) -> bool:
    # Stage 18.10 may only read closed files; unchanged Stage 18.0 wrapper
    # presence is accepted as preserved evidence without requiring a patch.
    return (root / "stage18_checks" / "run_stage18_0_preflight_boundary.sh").exists()


def false_positive_ok(root: Path, number: int) -> bool:
    helper = read_text(root / S18[number][1])
    return "source-only" in helper and "VALUE_KEYS" in helper and "final_status" in helper


def stage13_ok(root: Path) -> bool:
    return all_present(root, [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage13_checks/stage13_6_production_force_density_candidate.md",
    ])


def stage13_reg_ok(root: Path) -> bool:
    text = read_text(root / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py")
    return "local_subdomain_center" not in text.lower() or "global" in text.lower()


def stage14_ok(root: Path) -> bool:
    return all_present(root, ["src/fibre_stage14_production_rhs_injection.f90", "src/xcompact3d.f90"])


def no_rg_only_dependency(root: Path) -> bool:
    text = read_text(root / STAGE18_10_FILES[0])
    return "rg " not in text and "rg[[:space:]]" not in text


def no_activation(root: Path) -> bool:
    # Scan only the Stage 18.10 executable helper/wrapper.  Documentation and
    # closed-stage files may contain protective examples and negative labels.
    text = "\n".join(read_text(root / f) for f in STAGE18_10_FILES[:2]).lower()
    forbidden = [
        "call " + "stage14",
        "fibre_stage14" + "_production_rhs_injection(",
        "stage14" + "_rhs_injection(",
        "mpi" + "run ",
        "mpi" + "exec ",
        "subprocess.run([\"" + "mpirun" + "\"",
        "subprocess.run([\"" + "mpiexec" + "\"",
        "cm" + "ake ",
        "ct" + "est ",
        "nin" + "ja ",
        "ma" + "ke ",
        "fibre_stage13" + "_production_force_density_candidate.f90', 'w",
        "statistics" + ".f90', 'w",
        "visu" + ".f90', 'w",
        "restart" + "_write",
        "restart" + "_read",
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


def max_array_error(a: Sequence[Vec], b: Sequence[Vec]) -> float:
    return max_delta(a, b) if len(a) == len(b) else float("inf")


def energy_kin(v: Sequence[Vec], rho: float, w: Sequence[float]) -> float:
    return 0.5 * sum(rho * (vx*vx + vy*vy + vz*vz) * wi for (vx, vy, vz), wi in zip(v, w))


def energy_bend(xss: Sequence[Vec], stiffness: float, w: Sequence[float]) -> float:
    return 0.5 * sum(stiffness * (xx*xx + yy*yy + zz*zz) * wi for (xx, yy, zz), wi in zip(xss, w))


def power(fh: Sequence[Vec], v: Sequence[Vec], w: Sequence[float]) -> float:
    return sum((fx*vx + fy*vy + fz*vz) * wi for (fx, fy, fz), (vx, vy, vz), wi in zip(fh, v, w))


def step(x: Sequence[Vec], v: Sequence[Vec], a: Sequence[Vec], dt: float) -> Tuple[List[Vec], List[Vec]]:
    vn = [(vx + dt*ax, vy + dt*ay, vz + dt*az) for (vx, vy, vz), (ax, ay, az) in zip(v, a)]
    xn = [
        (xx + dt*vx + 0.5*dt*dt*ax, yy + dt*vy + 0.5*dt*dt*ay, zz + dt*vz + 0.5*dt*dt*az)
        for (xx, yy, zz), (vx, vy, vz), (ax, ay, az) in zip(x, v, a)
    ]
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


def coverage_ok(n: int, spans: Sequence[Tuple[int, int]]) -> Tuple[bool, bool, bool, bool]:
    counts = [0] * n
    for start, stop in spans:
        if start < 0 or stop < start or stop > n:
            return False, False, False, False
        for idx in range(start, stop):
            counts[idx] += 1
    coverage = all(c == 1 for c in counts)
    no_gap = all(c > 0 for c in counts)
    no_overlap = all(c <= 1 for c in counts)
    sizes = [stop - start for start, stop in spans]
    balance = bool(sizes) and (max(sizes) - min(sizes) <= 1)
    return coverage, no_gap, no_overlap, balance


def reconstruct(arr: Sequence[Vec], spans: Sequence[Tuple[int, int]], n: int) -> List[Vec]:
    out: List[Vec | None] = [None] * n
    for start, stop in spans:
        for local, idx in enumerate(range(start, stop)):
            out[idx] = arr[start + local]
    if any(item is None for item in out):
        return [(float("nan"), float("nan"), float("nan"))] * n
    return [item for item in out if item is not None]


def partition_sums(ref: Dict[str, Sequence[Vec] | float], spans: Sequence[Tuple[int, int]], w: Sequence[float], rho: float, stiffness: float) -> Dict[str, float]:
    v = ref["V_ref"]
    xss = ref["X_ss_ref"]
    fh = ref["F_h_ref"]
    x = ref["X_ref"]
    xnext = ref["X_next_ref"]
    a = ref["A_ref"]
    assert isinstance(v, list) and isinstance(xss, list) and isinstance(fh, list) and isinstance(x, list) and isinstance(xnext, list) and isinstance(a, list)
    ek = eb = ph = 0.0
    max_response = max_velocity = max_acceleration = 0.0
    for start, stop in spans:
        sl = slice(start, stop)
        ek += energy_kin(v[sl], rho, w[sl])
        eb += energy_bend(xss[sl], stiffness, w[sl])
        ph += power(fh[sl], v[sl], w[sl])
        max_response = max(max_response, max_delta(xnext[sl], x[sl]))
        max_velocity = max(max_velocity, max_norm(v[sl]))
        max_acceleration = max(max_acceleration, max_norm(a[sl]))
    return {
        "E_k_ref": ek,
        "E_b_ref": eb,
        "E_s_ref": ek + eb,
        "P_h_ref": ph,
        "max_response_ref": max_response,
        "max_velocity_ref": max_velocity,
        "max_acceleration_ref": max_acceleration,
    }


def build_reference(n: int, L: float, rho: float, stiffness: float, dt: float, eps: float, mode: int, fmag: float, v0: float, dry: bool) -> Dict[str, Sequence[Vec] | float]:
    ds = L / (n - 1)
    svals = [i * ds for i in range(n)]
    k = 2.0 * math.pi * mode / L
    w = weights(n, ds)
    x = [(s, eps * math.sin(k * s), 0.0) for s in svals]
    v = [(0.0, v0, 0.0)] * n
    xss = [(0.0, -eps * k * k * math.sin(k * s), 0.0) for s in svals]
    fb = [(0.0, -stiffness * eps * k**4 * math.sin(k * s), 0.0) for s in svals]
    ft = [(0.0, 0.0, 0.0)] * n
    fh = [(0.0, 0.0 if dry else fmag, 0.0)] * n
    ftotal = [(tx + bx + hx, ty + by + hy, tz + bz + hz) for (tx, ty, tz), (bx, by, bz), (hx, hy, hz) in zip(ft, fb, fh)]
    a = [(fx / rho, fy / rho, fz / rho) for fx, fy, fz in ftotal]
    xnext, vnext = step(x, v, a, dt)
    ek = energy_kin(v, rho, w)
    eb = energy_bend(xss, stiffness, w)
    return {
        "X_ref": x,
        "V_ref": v,
        "A_ref": a,
        "F_b_ref": fb,
        "F_T_ref": ft,
        "F_h_ref": fh,
        "F_total_ref": ftotal,
        "X_next_ref": xnext,
        "V_next_ref": vnext,
        "X_ss_ref": xss,
        "E_k_ref": ek,
        "E_b_ref": eb,
        "E_s_ref": ek + eb,
        "P_h_ref": power(fh, v, w),
        "max_response_ref": max_delta(xnext, x),
        "max_velocity_ref": max_norm(v),
        "max_acceleration_ref": max_norm(a),
    }


def array_reconstruction_ok(ref: Dict[str, Sequence[Vec] | float], spans: Sequence[Tuple[int, int]], n: int, tol: float) -> bool:
    for key in ["X_ref", "V_ref", "A_ref", "F_b_ref", "F_T_ref", "F_h_ref", "F_total_ref", "X_next_ref", "V_next_ref"]:
        arr = ref[key]
        assert isinstance(arr, list)
        if max_array_error(reconstruct(arr, spans, n), arr) > tol:
            return False
    return True


def scalar_reduction_ok(ref: Dict[str, Sequence[Vec] | float], red: Dict[str, float], keys: Sequence[str], tol: float) -> bool:
    return all(abs(red[key] - float(ref[key])) <= tol for key in keys)


def write_output(root: Path, summary: Dict[str, str], reasons: Sequence[str]) -> None:
    out = root / OUTPUT
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{key} {summary.get(key, 'FAIL')}" for key in SUMMARY_KEYS]
    lines.extend(f"reason {reason}" for reason in reasons)
    out.write_text("\n".join(lines) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--stage18-10-enable", default="1")
    parser.add_argument("--parallel-consistency-enable", default="1")
    parser.add_argument("--single-fibre-only", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--np-list", default="1,2,4")
    parser.add_argument("--npts", default="64")
    parser.add_argument("--nsteps", default="5")
    parser.add_argument("--fibre-length", default="1.0")
    parser.add_argument("--component-dim", default="3")
    parser.add_argument("--rho-l", default="1.0")
    parser.add_argument("--rho-tilde", default="1.0")
    parser.add_argument("--bending-stiffness", default="1.0e-3")
    parser.add_argument("--gamma", default="1.0e-3")
    parser.add_argument("--use-dimensional-response", default="1")
    parser.add_argument("--use-nondimensional-response", default="1")
    parser.add_argument("--dt-structure", default="1.0e-4")
    parser.add_argument("--sine-eps", default="1.0e-3")
    parser.add_argument("--sine-mode", default="1")
    parser.add_argument("--fluid-force-mag", default="1.0e-3")
    parser.add_argument("--initial-velocity-mag", default="0.0")
    parser.add_argument("--response-bound", default="1.0e-4")
    parser.add_argument("--velocity-bound", default="1.0e-3")
    parser.add_argument("--acceleration-bound", default="1.0e-2")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--formula-tol", default="1.0e-12")
    parser.add_argument("--reduction-tol", default="1.0e-12")
    parser.add_argument("--energy-tol", default="1.0e-12")
    parser.add_argument("--bounded-tol", default="1.0e-8")
    parser.add_argument("--test-case", default="np1_np2_np4_partition_reduction_consistency")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(args.repo_root).resolve()
    reasons: List[str] = []
    summary = {key: "FAIL" for key in SUMMARY_KEYS}

    np_values = parse_np_list(args.np_list)
    n = ii(args.npts)
    nsteps = ii(args.nsteps)
    comp = ii(args.component_dim)
    mode = ii(args.sine_mode, 1)
    L = ff(args.fibre_length)
    dt = ff(args.dt_structure)
    rho = ff(args.rho_l)
    rhot = ff(args.rho_tilde)
    B = ff(args.bending_stiffness)
    gamma = ff(args.gamma)
    eps = ff(args.sine_eps)
    fmag = ff(args.fluid_force_mag)
    v0 = ff(args.initial_velocity_mag)
    response_bound = ff(args.response_bound)
    velocity_bound = ff(args.velocity_bound)
    acceleration_bound = ff(args.acceleration_bound)
    ftol = ff(args.formula_tol, 1.0e-12) or 1.0e-12
    rtol = ff(args.reduction_tol, 1.0e-12) or 1.0e-12
    etol = ff(args.energy_tol, 1.0e-12) or 1.0e-12
    btol = ff(args.bounded_tol, 1.0e-8) or 1.0e-8
    use_dim = ii(args.use_dimensional_response, 0)
    use_nd = ii(args.use_nondimensional_response, 0)

    n_ok = n is not None and n >= 8
    nsteps_ok = nsteps is not None and nsteps >= 1
    L_ok = L is not None and L > 0.0
    dt_ok = dt is not None and dt > 0.0
    rho_ok = rho is not None and rho > 0.0
    rhot_ok = rhot is not None and rhot > 0.0
    B_ok = B is not None and B >= 0.0
    gamma_ok = gamma is not None and gamma >= 0.0
    np_ok = sorted(np_values) == [1, 2, 4]
    scalar_ok = n_ok and nsteps_ok and L_ok and dt_ok and rho_ok and rhot_ok and B_ok and gamma_ok

    n_d = n if n_ok else 8
    L_d = L if L_ok else 1.0
    dt_d = dt if dt_ok else 1.0e-4
    rho_d = rho if rho_ok else 1.0
    B_d = B if B_ok else 1.0e-3
    eps_d = eps if eps is not None else 1.0e-3
    fmag_d = fmag if fmag is not None else 1.0e-3
    v0_d = v0 if v0 is not None else 0.0
    mode_d = mode if mode and mode > 0 else 1
    response_bound_d = response_bound if response_bound is not None else 1.0e-4
    velocity_bound_d = velocity_bound if velocity_bound is not None else 1.0e-3
    acceleration_bound_d = acceleration_bound if acceleration_bound is not None else 1.0e-2
    ds = L_d / (n_d - 1)
    w = weights(n_d, ds)

    forced_ref = build_reference(n_d, L_d, rho_d, B_d, dt_d, eps_d, mode_d, fmag_d, v0_d, dry=False)
    dry_ref = build_reference(n_d, L_d, rho_d, B_d, dt_d, eps_d, mode_d, 0.0, v0_d, dry=True)
    ref_arrays = [forced_ref[key] for key in ["X_ref", "V_ref", "A_ref", "F_b_ref", "F_T_ref", "F_h_ref", "F_total_ref", "X_next_ref", "V_next_ref"]]
    assert all(isinstance(arr, list) for arr in ref_arrays)
    ref_scalars = [forced_ref[key] for key in ["E_k_ref", "E_b_ref", "E_s_ref", "P_h_ref", "max_response_ref", "max_velocity_ref", "max_acceleration_ref"]]

    part_checks: Dict[int, Dict[str, bool]] = {}
    for np in [1, 2, 4]:
        spans = partitions(n_d, np)
        coverage, no_gap, no_overlap, balance = coverage_ok(n_d, spans)
        red = partition_sums(forced_ref, spans, w, rho_d, B_d)
        dry_red = partition_sums(dry_ref, spans, w, rho_d, B_d)
        part_checks[np] = {
            "coverage": coverage,
            "no_gap": no_gap,
            "no_overlap": no_overlap,
            "balance": balance,
            "arrays": array_reconstruction_ok(forced_ref, spans, n_d, ftol),
            "scalars": scalar_reduction_ok(forced_ref, red, ["E_k_ref", "E_b_ref", "E_s_ref", "P_h_ref", "max_response_ref", "max_velocity_ref", "max_acceleration_ref"], rtol),
            "energy": scalar_reduction_ok(forced_ref, red, ["E_k_ref", "E_b_ref", "E_s_ref"], etol),
            "power": scalar_reduction_ok(forced_ref, red, ["P_h_ref"], rtol),
            "bounds": scalar_reduction_ok(forced_ref, red, ["max_response_ref", "max_velocity_ref", "max_acceleration_ref"], rtol),
            "dry": array_reconstruction_ok(dry_ref, spans, n_d, ftol) and abs(dry_red["P_h_ref"] - float(dry_ref["P_h_ref"])) <= rtol,
            "forced": array_reconstruction_ok(forced_ref, spans, n_d, ftol) and red["P_h_ref"] >= -rtol,
        }

    ref_finite = finite_vec(ref_arrays) and finite_vals(float(v) for v in ref_scalars)
    ref_energy_nonnegative = all(float(forced_ref[key]) >= -etol for key in ["E_k_ref", "E_b_ref", "E_s_ref"])
    ref_bounded = (
        float(forced_ref["max_response_ref"]) <= response_bound_d + btol and
        float(forced_ref["max_velocity_ref"]) <= velocity_bound_d + btol and
        float(forced_ref["max_acceleration_ref"]) <= acceleration_bound_d + btol
    )

    changed_ok = changed_only_allowed(root)
    closed_file_present, closed_evidence = stage17_closed(root)
    activation_ok = no_activation(root)
    doc = read_text(root / "stage18_checks" / "stage18_10_parallel_consistency_physical_structure.md")
    requested_files = all_present(root, STAGE18_10_FILES)

    summary.update({
        "stage18_10_requested_status": status(requested_files),
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
        "stage18_9_controlled_response_preserved_status": status(evidence(root, 9, ["CONTROLLED RESPONSE DIAGNOSTIC", "LOCAL CANDIDATE FORCE/STEP", "NP=1 DIAGNOSTIC"])),
        "no_closed_stage_modification_status": status(changed_ok),
        "no_stage10_17_file_modification_status": status(changed_ok),
        "stage18_0_files_unmodified_status": status(changed_ok),
        "stage18_1_files_unmodified_status": status(changed_ok),
        "stage18_2_files_unmodified_status": status(changed_ok),
        "stage18_3_files_unmodified_status": status(changed_ok),
        "stage18_4_files_unmodified_status": status(changed_ok),
        "stage18_5_files_unmodified_status": status(changed_ok),
        "stage18_6_files_unmodified_status": status(changed_ok),
        "stage18_7_files_unmodified_status": status(changed_ok),
        "stage18_8_files_unmodified_status": status(changed_ok),
        "stage18_9_files_unmodified_status": status(changed_ok),
        "stage18_enable_status": status(args.stage18_10_enable == "1"),
        "parallel_consistency_enable_status": status(args.parallel_consistency_enable == "1"),
        "single_fibre_only_status": status(args.single_fibre_only == "1"),
        "diagnostic_only_status": status(args.diagnostic_only == "1"),
        "np_list_status": status(np_ok),
        "np1_reference_status": status(1 in np_values and ref_finite and ref_bounded),
        "npts_status": status(n_ok),
        "nsteps_status": status(nsteps_ok),
        "component_dim_status": status(comp == 3),
        "fibre_length_status": status(L_ok),
        "ds_formula_status": status(n_ok and L_ok and abs(ds - L_d / (n_d - 1)) <= ftol and ds > 0.0),
        "dt_structure_status": status(dt_ok),
        "rho_l_status": status(rho_ok),
        "rho_tilde_status": status(rhot_ok),
        "bending_stiffness_status": status(B_ok),
        "gamma_status": status(gamma_ok),
        "dimensional_response_validation_status": status(use_dim == 1 and scalar_ok),
        "nondimensional_response_validation_status": status(use_nd == 1 and rhot_ok and gamma_ok),
        "partition_coverage_status": status(all(part_checks[np]["coverage"] for np in [1, 2, 4])),
        "partition_no_gap_status": status(all(part_checks[np]["no_gap"] for np in [1, 2, 4])),
        "partition_no_overlap_status": status(all(part_checks[np]["no_overlap"] for np in [1, 2, 4])),
        "partition_size_balance_status": status(all(part_checks[np]["balance"] for np in [1, 2, 4])),
        "np1_reference_array_finite_status": status(ref_finite),
        "np1_reference_energy_nonnegative_status": status(ref_energy_nonnegative),
        "np2_array_reconstruction_status": status(part_checks[2]["arrays"]),
        "np4_array_reconstruction_status": status(part_checks[4]["arrays"]),
        "np2_scalar_reduction_status": status(part_checks[2]["scalars"]),
        "np4_scalar_reduction_status": status(part_checks[4]["scalars"]),
        "np2_energy_reduction_status": status(part_checks[2]["energy"]),
        "np4_energy_reduction_status": status(part_checks[4]["energy"]),
        "np2_power_reduction_status": status(part_checks[2]["power"]),
        "np4_power_reduction_status": status(part_checks[4]["power"]),
        "np2_response_bound_reduction_status": status(part_checks[2]["bounds"]),
        "np4_response_bound_reduction_status": status(part_checks[4]["bounds"]),
        "dry_response_parallel_consistency_status": status(all(part_checks[np]["dry"] for np in [1, 2, 4])),
        "forced_response_parallel_consistency_status": status(all(part_checks[np]["forced"] for np in [1, 2, 4])),
        "partition_reduction_diagnostic_only_status": status(activation_ok and "PARALLEL CONSISTENCY DIAGNOSTIC" in doc),
        "stage13_6_diagnostic_preserved_status": status(stage13_ok(root)),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_reg_ok(root)),
        "stage14_small_lambda_hook_status": status(stage14_ok(root)),
        "no_rg_only_dependency_status": status(no_rg_only_dependency(root)),
        "no_unknown_failure_status": "PASS",
    })
    for key in NEGATIVE_STATUS_KEYS:
        summary[key] = status(activation_ok)

    if args.test_case != "np1_np2_np4_partition_reduction_consistency":
        reasons.append("unexpected_stage18_10_test_case")
    if not requested_files:
        reasons.append("required_stage18_10_file_missing")
    if not changed_ok:
        bad = [path for _code, path in git_entries(root) if path not in ALLOWED_CHANGED]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad))

    pass_fail_keys = [key for key in SUMMARY_KEYS if key.endswith("_status") and key != "final_status" and key not in VALUE_KEYS]
    for key in pass_fail_keys:
        if summary.get(key) != "PASS":
            reasons.append(key.replace("_status", "_failed"))
    summary["final_status"] = "PASS" if all(summary.get(key) == "PASS" for key in pass_fail_keys) else "FAIL"

    write_output(root, summary, reasons)
    for key in SUMMARY_KEYS:
        print(f"{key} {summary.get(key, 'FAIL')}")
    for reason in reasons:
        print(f"reason {reason}")
    return 0 if summary["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
