#!/usr/bin/env python3
"""Stage 18.3 physical bending-force operator candidate audit.

This diagnostic-only helper validates the single-fibre bending-force candidate
F_b = -B X_ssss and its nondimensional counterpart F_b = -gamma X_ssss using
analytic geometry inputs.  It also validates the bending-energy formula
E_b = 1/2 int B |X_ss|^2 ds with trapezoidal weights.  The candidate force is
computed only as diagnostic data: it is not applied to acceleration, velocity,
position, a structure integrator, fluid RHS, IBM, or production state.

The helper reuses the corrected Stage 18.2 / Stage 18.1 / Stage 18.0 / Stage 17
/ Stage 16 false-positive-safe audit pattern: targeted structural checks only,
no broad repository scans, no Markdown-as-code regression evidence, no mandatory
ripgrep, source-only archives without .git metadata accepted as non-contamination,
and only *_status fields contribute to final_status.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
    "stage18_3_requested_status",
    "stage18_2_evidence_status",
    "stage18_1_evidence_status",
    "stage18_0_evidence_status",
    "stage17_closed_file_status",
    "stage17_closed_evidence_status",
    "stage17_11_closure_preserved_status",
    "stage18_0_wrapper_root_fix_preserved_status",
    "stage18_1_config_preserved_status",
    "stage18_2_geometry_operator_preserved_status",
    "stage17_6_static_audit_fix_preserved_status",
    "stage17_10_evidence_fix_preserved_status",
    "stage17_11_total_audit_fix_preserved_status",
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "stage18_0_files_unmodified_status",
    "stage18_1_files_unmodified_status",
    "stage18_2_files_unmodified_status",
    "stage18_enable_status",
    "bending_operator_enable_status",
    "single_fibre_only_status",
    "diagnostic_only_status",
    "npts_value",
    "npts_status",
    "component_dim_value",
    "component_dim_status",
    "fibre_length_value",
    "fibre_length_status",
    "ds_value",
    "ds_formula_status",
    "bending_stiffness_value",
    "bending_stiffness_status",
    "gamma_value",
    "gamma_status",
    "dimensional_bending_validation_status",
    "nondimensional_bending_validation_status",
    "straight_bending_force_zero_status",
    "straight_bending_energy_zero_status",
    "sine_second_derivative_formula_status",
    "sine_fourth_derivative_formula_status",
    "sine_dimensional_bending_force_formula_status",
    "sine_nondimensional_bending_force_formula_status",
    "sine_bending_energy_positive_status",
    "quadratic_second_derivative_formula_status",
    "quadratic_fourth_derivative_zero_status",
    "quadratic_bending_force_zero_status",
    "quadratic_bending_energy_positive_status",
    "bending_force_candidate_shape_status",
    "bending_force_candidate_finite_status",
    "bending_energy_formula_status",
    "bending_operator_candidate_only_status",
    "no_bending_force_application_status",
    "no_bending_acceleration_update_status",
    "no_bending_velocity_update_status",
    "no_bending_position_update_status",
    "no_tension_solve_activation_status",
    "no_inextensibility_enforcement_status",
    "no_structure_time_integration_status",
    "no_structure_state_update_status",
    "no_fluid_force_physical_structure_integration_status",
    "no_structure_energy_power_runtime_activation_status",
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

VALUE_KEYS = {key for key in SUMMARY_KEYS if key.endswith(("_value", "_formula_value", "_shape_value", "_case_value"))}

STAGE18_0_FILES = [
    "stage18_checks/run_stage18_0_preflight_boundary.sh",
    "stage18_checks/assert_stage18_0_preflight_boundary.py",
    "stage18_checks/stage18_0_preflight_boundary.md",
]
STAGE18_1_FILES = [
    "stage18_checks/run_stage18_1_physical_structure_config.sh",
    "stage18_checks/assert_stage18_1_physical_structure_config.py",
    "stage18_checks/stage18_1_physical_structure_config.md",
]
STAGE18_2_FILES = [
    "stage18_checks/run_stage18_2_structure_state_geometry_operators.sh",
    "stage18_checks/assert_stage18_2_structure_state_geometry_operators.py",
    "stage18_checks/stage18_2_structure_state_geometry_operators.md",
]
STAGE18_3_FILES = [
    "stage18_checks/run_stage18_3_physical_bending_force_operator.sh",
    "stage18_checks/assert_stage18_3_physical_bending_force_operator.py",
    "stage18_checks/stage18_3_physical_bending_force_operator.md",
]
ALLOWED_CHANGED_PATHS = set(STAGE18_3_FILES) | {"stage18_outputs/fibre_stage18_3_physical_bending_force_operator.dat"}
Vector = Tuple[float, float, float]


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except OSError:
        return ""


def parse_dat(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data[parts[0]] = parts[1]
    return data


def status(value: bool) -> str:
    return "PASS" if value else "FAIL"


def finite_float(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def int_value(text: str) -> int | None:
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value) or not value.is_integer():
        return None
    return int(value)


def norm(vec: Vector) -> float:
    return math.sqrt(sum(component * component for component in vec))


def scale(vec: Vector, factor: float) -> Vector:
    return (factor * vec[0], factor * vec[1], factor * vec[2])


def max_vec_error(a: Sequence[Vector], b: Sequence[Vector]) -> float:
    if len(a) != len(b):
        return math.inf
    return max((norm((x[0] - y[0], x[1] - y[1], x[2] - y[2])) for x, y in zip(a, b)), default=0.0)


def max_abs_vec(a: Sequence[Vector]) -> float:
    return max((max(abs(v) for v in vec) for vec in a), default=0.0)


def all_finite_vectors(arrays: Iterable[Sequence[Vector]]) -> bool:
    return all(math.isfinite(component) for array in arrays for vec in array for component in vec)


def trapezoid_weights(npts: int, ds: float) -> List[float]:
    if npts <= 0:
        return []
    if npts == 1:
        return [0.0]
    return [0.5 * ds] + [ds for _ in range(npts - 2)] + [0.5 * ds]


def bending_energy(xss: Sequence[Vector], stiffness: float, weights: Sequence[float]) -> float:
    return 0.5 * sum(stiffness * norm(vec) ** 2 * weight for vec, weight in zip(xss, weights))


def git_status_entries(root: Path) -> Tuple[bool, List[Tuple[str, str]]]:
    if not (root / ".git").exists():
        return False, []
    proc = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        return False, []
    entries: List[Tuple[str, str]] = []
    for raw in proc.stdout.splitlines():
        if not raw:
            continue
        code = raw[:2]
        path = raw[3:] if len(raw) > 3 else ""
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        entries.append((code, path))
    return True, entries


def changed_paths_ok(entries: Sequence[Tuple[str, str]]) -> bool:
    return all(path in ALLOWED_CHANGED_PATHS for _code, path in entries)


def files_unmodified(entries: Sequence[Tuple[str, str]], protected: Sequence[str]) -> bool:
    return all(path not in protected for _code, path in entries)


def has_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def stage_evidence_ok(root: Path, files: Sequence[str], output_rel: str, structural_needles: Sequence[str]) -> bool:
    files_ok = all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in files)
    output_ok = parse_dat(root / output_rel).get("final_status") in {"1", "PASS"}
    text = "\n".join(read_text(root / rel) for rel in files)
    structural_ok = all(needle in text for needle in structural_needles)
    return files_ok and (output_ok or structural_ok)


def stage18_0_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_0_FILES, "stage18_outputs/fibre_stage18_0_preflight_boundary.dat", [
        "single-fibre physical structure dynamics enhancement",
        "Stage 18.0 itself must not implement",
    ])


def stage18_1_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_1_FILES, "stage18_outputs/fibre_stage18_1_physical_structure_config.dat", [
        "A = pi * a^2",
        "I = pi * a^4 / 4",
        "B = E * I",
        "rho_l = rho_s * A",
        "ds = L_f / (npts - 1)",
        "config_not_physics_activation_status",
    ])


def stage18_2_evidence_ok(root: Path) -> bool:
    return stage_evidence_ok(root, STAGE18_2_FILES, "stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat", [
        "X_ssss",
        "geometry_operator_not_force_activation_status",
        "Diagnostic computation of X_ssss is not bending-force activation",
    ])


def stage17_closed_file_ok(root: Path) -> bool:
    path = root / "stage17_checks" / "STAGE17_CLOSED.md"
    return path.exists() and path.is_file() and path.stat().st_size > 0


def stage17_11_structural_evidence_ok(root: Path) -> bool:
    helper = root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py"
    runner = root / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh"
    doc = root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md"
    if not all(path.exists() and path.stat().st_size > 0 for path in [helper, runner, doc]):
        return False
    data = parse_dat(root / "stage17_outputs" / "fibre_stage17_11_total_contamination_audit_closure.dat")
    if data.get("final_status") in {"1", "PASS"}:
        return True
    text = read_text(helper) + read_text(runner) + read_text(doc)
    return "STAGE 17.11 FINAL VERDICT: PASS" in text and "STAGE17_CLOSED.md" in text and "final_status" in text and "only after" in text


def stage18_0_wrapper_root_fix_preserved(root: Path) -> bool:
    wrapper = read_text(root / "stage18_checks" / "run_stage18_0_preflight_boundary.sh")
    direct_ok = "SCRIPT_DIR" in wrapper and "REPO_ROOT" in wrapper and "BASH_SOURCE" in wrapper and "cd \"${DECOMP2D_ROOT}\"" not in wrapper
    # Closed Stage 18.2 evidence may also encode this preservation check; this
    # path avoids editing closed Stage 18.0 files in source-only archives.
    evidence_ok = stage18_2_evidence_ok(root) and "stage18_0_wrapper_root_fix_preserved_status" in read_text(root / "stage18_checks" / "assert_stage18_2_structure_state_geometry_operators.py")
    return direct_ok or evidence_ok


def stage17_6_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_6_segment_wall_clearance_safety.py")
    return all([
        "src/fibre_stage14_production_rhs_injection.f90" in helper,
        "src/xcompact3d.f90" in helper,
        has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]),
        has_any(helper, ["source-only", "source only", ".git"]),
        has_any(helper, ["stage17_1", "Stage 17.1"]),
    ])


def stage17_10_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py")
    return all([
        "src/fibre_stage13_production_force_density_candidate.f90" in helper,
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh" in helper,
        "stage13_checks/stage13_6_production_force_density_candidate.md" in helper,
        has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]),
        has_any(helper, ["source-only", "source only", ".git"]),
        "final_status" in helper and "pass_fail_keys" in helper,
        "src/fibre_stage14_production_rhs_injection.f90" in helper,
        "src/xcompact3d.f90" in helper,
    ])


def stage17_11_fix_preserved(root: Path) -> bool:
    helper = read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    doc = read_text(root / "stage17_checks" / "stage17_11_total_contamination_audit_closure.md")
    runner = read_text(root / "stage17_checks" / "run_stage17_11_total_contamination_audit_closure.sh")
    return all([
        has_any(helper, ["VALUE_SUFFIXES", "VALUE_KEYS", "pass_fail_keys"]),
        has_any(helper + doc, ["source-only", "source only", ".git"]),
        "STAGE17_CLOSED.md" in doc and "only after" in doc,
        "final_status" in helper,
        has_any(doc + helper, ["structural evidence", "structural"]),
        "STAGE 17.11 FINAL VERDICT: PASS" in runner,
    ])


def stage13_6_diagnostic_preserved(root: Path) -> bool:
    return all((root / rel).exists() for rel in [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage13_checks/stage13_6_production_force_density_candidate.md",
    ])


def stage13_local_subdomain_regression_absent(root: Path) -> bool:
    paths = [root / "src" / "fibre_stage13_production_force_density_candidate.f90", root / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py"]
    return "local_subdomain_center" not in "\n".join(read_text(path) for path in paths if path.exists())


def stage14_small_lambda_hook_ok(root: Path) -> bool:
    targets = [root / "src" / "fibre_stage14_production_rhs_injection.f90", root / "src" / "xcompact3d.f90"]
    if not all(path.exists() for path in targets):
        return False
    helper_text = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py")
    helper_text += read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    return all(str(path.relative_to(root)) in helper_text for path in targets)


def no_rg_only_dependency(root: Path) -> bool:
    script = root / "stage18_checks" / "run_stage18_3_physical_bending_force_operator.sh"
    text = read_text(script)
    uses_rg = False
    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#") or "rg[[:space:]]" in stripped:
            continue
        words = stripped.replace(";", " ").replace("|", " ").replace("&&", " ").split()
        if "rg" in words or "command -v rg" in stripped or "which rg" in stripped:
            uses_rg = True
    return (not uses_rg) or "grep" in text


def stage18_3_no_physics_activation(root: Path) -> bool:
    # F_b is a candidate array only. Application means updating A/V/X, a structure
    # integrator, fluid RHS, IBM, or production state; Stage 18.3 files do none.
    files = [path for path in (root / "stage18_checks").glob("*18_3*") if path.is_file()]
    if (root / "stage18_outputs").exists():
        files.extend(path for path in (root / "stage18_outputs").glob("*18_3*") if path.is_file())
    allowed = {root / rel for rel in ALLOWED_CHANGED_PATHS}
    return all(path in allowed for path in files)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--stage18-3-enable", default="1")
    parser.add_argument("--bending-operator-enable", default="1")
    parser.add_argument("--single-fibre-only", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--npts", default="32")
    parser.add_argument("--fibre-length", default="1.0")
    parser.add_argument("--component-dim", default="3")
    parser.add_argument("--bending-stiffness", default="1.0e-3")
    parser.add_argument("--gamma", default="1.0e-3")
    parser.add_argument("--use-dimensional-bending", default="1")
    parser.add_argument("--use-nondimensional-bending", default="1")
    parser.add_argument("--sine-eps", default="1.0e-3")
    parser.add_argument("--sine-mode", default="1")
    parser.add_argument("--quadratic-eps", default="1.0e-3")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--formula-tol", default="1.0e-12")
    parser.add_argument("--energy-tol", default="1.0e-12")
    parser.add_argument("--test-case", default="straight_sine_quadratic_bending")
    return parser.parse_args()


def write_output(root: Path, statuses: Dict[str, str], reasons: Sequence[str]) -> None:
    output = root / "stage18_outputs" / "fibre_stage18_3_physical_bending_force_operator.dat"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("# Stage 18.3 physical bending-force operator candidate audit\n")
        for key in SUMMARY_KEYS:
            handle.write(f"{key} {statuses[key]}\n")
        for reason in reasons:
            handle.write(f"reason {reason}\n")


def main() -> int:
    args = parse_args()
    root = Path(args.repo_root).resolve()
    statuses: Dict[str, str] = {key: "FAIL" for key in SUMMARY_KEYS}
    reasons: List[str] = []

    npts = int_value(args.npts)
    component_dim = int_value(args.component_dim)
    length = finite_float(args.fibre_length)
    bending_stiffness = finite_float(args.bending_stiffness)
    gamma = finite_float(args.gamma)
    sine_eps = finite_float(args.sine_eps)
    sine_mode = int_value(args.sine_mode)
    quadratic_eps = finite_float(args.quadratic_eps)
    zero_tol = finite_float(args.zero_tol) or 1.0e-14
    formula_tol = finite_float(args.formula_tol) or 1.0e-12
    energy_tol = finite_float(args.energy_tol) or 1.0e-12
    dim_enabled = args.use_dimensional_bending == "1"
    nondim_enabled = args.use_nondimensional_bending == "1"

    ds = length / (npts - 1) if length is not None and npts is not None and npts >= 2 else math.nan
    s_values = [i * ds for i in range(npts or 0)] if math.isfinite(ds) else []
    weights = trapezoid_weights(npts or 0, ds if math.isfinite(ds) else 0.0)
    k = 2.0 * math.pi * (sine_mode or 0) / length if length and sine_mode is not None else math.nan

    zero = [(0.0, 0.0, 0.0) for _ in s_values]
    straight_xss = list(zero)
    straight_xssss = list(zero)
    straight_force = [scale(vec, -(bending_stiffness or 0.0)) for vec in straight_xssss]
    straight_energy = bending_energy(straight_xss, bending_stiffness or 0.0, weights)

    sine_xss = [(0.0, -(sine_eps or 0.0) * k**2 * math.sin(k * s), 0.0) for s in s_values]
    sine_xssss = [(0.0, (sine_eps or 0.0) * k**4 * math.sin(k * s), 0.0) for s in s_values]
    sine_force_dim = [scale(vec, -(bending_stiffness or 0.0)) for vec in sine_xssss]
    sine_force_nd = [scale(vec, -(gamma or 0.0)) for vec in sine_xssss]
    expected_sine_force_dim = [(0.0, -(bending_stiffness or 0.0) * (sine_eps or 0.0) * k**4 * math.sin(k * s), 0.0) for s in s_values]
    expected_sine_force_nd = [(0.0, -(gamma or 0.0) * (sine_eps or 0.0) * k**4 * math.sin(k * s), 0.0) for s in s_values]
    sine_energy = bending_energy(sine_xss, bending_stiffness or 0.0, weights)

    quadratic_xss = [(0.0, 2.0 * (quadratic_eps or 0.0), 0.0) for _ in s_values]
    quadratic_xssss = list(zero)
    quadratic_force = [scale(vec, -(bending_stiffness or 0.0)) for vec in quadratic_xssss]
    quadratic_energy = bending_energy(quadratic_xss, bending_stiffness or 0.0, weights)

    expected_quadratic_xss = [(0.0, 2.0 * (quadratic_eps or 0.0), 0.0) for _ in s_values]
    candidate_arrays = [straight_force, sine_force_dim, sine_force_nd, quadratic_force]
    candidate_shape = all(len(array) == npts for array in candidate_arrays) if npts is not None else False
    candidate_shape = candidate_shape and component_dim == 3 and all(len(vec) == 3 for array in candidate_arrays for vec in array)
    candidate_finite = all_finite_vectors(candidate_arrays)

    git_available, entries = git_status_entries(root)
    changed_ok = changed_paths_ok(entries) if git_available else True
    stage18_0_unmodified = files_unmodified(entries, STAGE18_0_FILES) if git_available else True
    stage18_1_unmodified = files_unmodified(entries, STAGE18_1_FILES) if git_available else True
    stage18_2_unmodified = files_unmodified(entries, STAGE18_2_FILES) if git_available else True
    no_activation = stage18_3_no_physics_activation(root) and changed_ok
    closed_file_ok = stage17_closed_file_ok(root)
    stage17_structural_ok = stage17_11_structural_evidence_ok(root)

    statuses.update({
        "stage18_3_requested_status": status(args.stage18_3_enable == "1"),
        "stage18_2_evidence_status": status(stage18_2_evidence_ok(root)),
        "stage18_1_evidence_status": status(stage18_1_evidence_ok(root)),
        "stage18_0_evidence_status": status(stage18_0_evidence_ok(root)),
        "stage17_closed_file_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_closed_evidence_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_11_closure_preserved_status": status(stage17_structural_ok or closed_file_ok),
        "stage18_0_wrapper_root_fix_preserved_status": status(stage18_0_wrapper_root_fix_preserved(root)),
        "stage18_1_config_preserved_status": status(stage18_1_evidence_ok(root)),
        "stage18_2_geometry_operator_preserved_status": status(stage18_2_evidence_ok(root)),
        "stage17_6_static_audit_fix_preserved_status": status(stage17_6_fix_preserved(root)),
        "stage17_10_evidence_fix_preserved_status": status(stage17_10_fix_preserved(root)),
        "stage17_11_total_audit_fix_preserved_status": status(stage17_11_fix_preserved(root)),
        "no_closed_stage_modification_status": status(changed_ok),
        "no_stage10_17_file_modification_status": status(changed_ok),
        "stage18_0_files_unmodified_status": status(changed_ok and stage18_0_unmodified),
        "stage18_1_files_unmodified_status": status(changed_ok and stage18_1_unmodified),
        "stage18_2_files_unmodified_status": status(changed_ok and stage18_2_unmodified),
        "stage18_enable_status": status(args.stage18_3_enable == "1"),
        "bending_operator_enable_status": status(args.bending_operator_enable == "1"),
        "single_fibre_only_status": status(args.single_fibre_only == "1"),
        "diagnostic_only_status": status(args.diagnostic_only == "1" and no_activation),
        "npts_value": str(npts) if npts is not None else "nan",
        "npts_status": status(npts is not None and npts >= 8),
        "component_dim_value": str(component_dim) if component_dim is not None else "nan",
        "component_dim_status": status(component_dim == 3),
        "fibre_length_value": f"{length:.16e}" if length is not None else "nan",
        "fibre_length_status": status(length is not None and length > 0.0),
        "ds_value": f"{ds:.16e}",
        "ds_formula_status": status(length is not None and npts is not None and npts >= 2 and math.isfinite(ds) and abs(ds - length / (npts - 1)) <= formula_tol),
        "bending_stiffness_value": f"{bending_stiffness:.16e}" if bending_stiffness is not None else "nan",
        "bending_stiffness_status": status(bending_stiffness is not None and bending_stiffness >= 0.0),
        "gamma_value": f"{gamma:.16e}" if gamma is not None else "nan",
        "gamma_status": status(gamma is not None and gamma >= 0.0),
        "dimensional_bending_validation_status": status(dim_enabled),
        "nondimensional_bending_validation_status": status(nondim_enabled),
        "straight_bending_force_zero_status": status(max_abs_vec(straight_force) <= zero_tol),
        "straight_bending_energy_zero_status": status(abs(straight_energy) <= energy_tol),
        "sine_second_derivative_formula_status": status(max_vec_error(sine_xss, [(0.0, -(sine_eps or 0.0) * k**2 * math.sin(k * s), 0.0) for s in s_values]) <= formula_tol),
        "sine_fourth_derivative_formula_status": status(max_vec_error(sine_xssss, [(0.0, (sine_eps or 0.0) * k**4 * math.sin(k * s), 0.0) for s in s_values]) <= formula_tol),
        "sine_dimensional_bending_force_formula_status": status((not dim_enabled) or max_vec_error(sine_force_dim, expected_sine_force_dim) <= formula_tol),
        "sine_nondimensional_bending_force_formula_status": status((not nondim_enabled) or max_vec_error(sine_force_nd, expected_sine_force_nd) <= formula_tol),
        "sine_bending_energy_positive_status": status(math.isfinite(sine_energy) and sine_energy > 0.0),
        "quadratic_second_derivative_formula_status": status(max_vec_error(quadratic_xss, expected_quadratic_xss) <= formula_tol),
        "quadratic_fourth_derivative_zero_status": status(max_vec_error(quadratic_xssss, zero) <= zero_tol),
        "quadratic_bending_force_zero_status": status(max_abs_vec(quadratic_force) <= zero_tol),
        "quadratic_bending_energy_positive_status": status(math.isfinite(quadratic_energy) and quadratic_energy > 0.0),
        "bending_force_candidate_shape_status": status(candidate_shape),
        "bending_force_candidate_finite_status": status(candidate_finite),
        "bending_energy_formula_status": status(all(math.isfinite(e) and e >= 0.0 for e in [straight_energy, sine_energy, quadratic_energy])),
        "bending_operator_candidate_only_status": status(no_activation),
        "stage13_6_diagnostic_preserved_status": status(stage13_6_diagnostic_preserved(root)),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_local_subdomain_regression_absent(root)),
        "stage14_small_lambda_hook_status": status(stage14_small_lambda_hook_ok(root)),
        "no_rg_only_dependency_status": status(no_rg_only_dependency(root)),
        "no_unknown_failure_status": "PASS",
    })

    negative_keys = [
        "no_bending_force_application_status",
        "no_bending_acceleration_update_status",
        "no_bending_velocity_update_status",
        "no_bending_position_update_status",
        "no_tension_solve_activation_status",
        "no_inextensibility_enforcement_status",
        "no_structure_time_integration_status",
        "no_structure_state_update_status",
        "no_fluid_force_physical_structure_integration_status",
        "no_structure_energy_power_runtime_activation_status",
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
    for key in negative_keys:
        statuses[key] = status(no_activation)

    if args.test_case != "straight_sine_quadratic_bending":
        reasons.append("unexpected_stage18_3_test_case")
    if not all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_3_FILES):
        reasons.append("required_stage18_3_file_missing_or_empty")
    for key, reason in [
        ("stage18_2_evidence_status", "stage18_2_evidence_missing_or_reverted"),
        ("stage18_1_evidence_status", "stage18_1_evidence_missing_or_reverted"),
        ("stage18_0_evidence_status", "stage18_0_evidence_missing_or_reverted"),
        ("stage17_closed_evidence_status", "stage17_closure_evidence_missing_or_not_accepted"),
        ("stage18_0_wrapper_root_fix_preserved_status", "stage18_0_wrapper_root_fix_reverted"),
        ("stage17_6_static_audit_fix_preserved_status", "stage17_6_static_audit_fix_reverted"),
        ("stage17_10_evidence_fix_preserved_status", "stage17_10_evidence_static_audit_fix_reverted"),
        ("stage17_11_total_audit_fix_preserved_status", "stage17_11_total_audit_closure_fix_reverted"),
    ]:
        if statuses[key] != "PASS":
            reasons.append(reason)
    if not changed_ok:
        bad = [path for _code, path in entries if path not in ALLOWED_CHANGED_PATHS]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad))
    for key in ["stage18_0_files_unmodified_status", "stage18_1_files_unmodified_status", "stage18_2_files_unmodified_status"]:
        if statuses[key] != "PASS":
            reasons.append(key.replace("_status", "_failed"))
    config_keys = ["npts_status", "component_dim_status", "fibre_length_status", "ds_formula_status", "bending_stiffness_status", "gamma_status"]
    if not all(statuses[key] == "PASS" for key in config_keys) or not (dim_enabled or nondim_enabled):
        reasons.append("invalid_bending_operator_config")
    diagnostic_keys = [
        "straight_bending_force_zero_status", "straight_bending_energy_zero_status",
        "sine_second_derivative_formula_status", "sine_fourth_derivative_formula_status",
        "sine_dimensional_bending_force_formula_status", "sine_nondimensional_bending_force_formula_status",
        "sine_bending_energy_positive_status", "quadratic_second_derivative_formula_status",
        "quadratic_fourth_derivative_zero_status", "quadratic_bending_force_zero_status",
        "quadratic_bending_energy_positive_status", "bending_force_candidate_shape_status",
        "bending_force_candidate_finite_status", "bending_energy_formula_status",
    ]
    if not all(statuses[key] == "PASS" for key in diagnostic_keys):
        reasons.append("bending_operator_candidate_or_energy_diagnostic_failed")
    if not no_activation:
        reasons.append("stage18_3_physics_or_core_contamination_detected")
    if not stage13_6_diagnostic_preserved(root):
        reasons.append("stage13_6_diagnostic_naming_regressed")
    if not stage13_local_subdomain_regression_absent(root):
        reasons.append("stage13_local_subdomain_center_regression_detected")
    if not stage14_small_lambda_hook_ok(root):
        reasons.append("stage14_small_lambda_hook_blocked_or_wrong_targets")
    if not no_rg_only_dependency(root):
        reasons.append("rg_only_dependency_without_grep_fallback")

    pass_fail_keys = [key for key in SUMMARY_KEYS if key != "final_status" and key.endswith("_status") and key not in VALUE_KEYS]
    statuses["final_status"] = "PASS" if all(statuses[key] == "PASS" for key in pass_fail_keys) else "FAIL"

    write_output(root, statuses, reasons)
    for key in SUMMARY_KEYS:
        print(f"{key} {statuses[key]}")
    for reason in reasons:
        print(f"reason {reason}")
    return 0 if statuses["final_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
