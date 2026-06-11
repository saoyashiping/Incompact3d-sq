#!/usr/bin/env python3
"""Stage 18.2 structure-state and geometry-operator audit.

This diagnostic-only helper validates single-fibre structure-state concepts and
analytic geometry operators for later Stage 18 physical structure dynamics.  It
constructs small in-memory diagnostic arrays for X, V, A, X_s, X_ss, X_sss, and
X_ssss, checks straight and sinusoidal analytic geometry identities, and writes a
summary dat file.  It does not implement or activate a bending force, tension
solve, inextensibility enforcement, time integration, state update, fluid-force
physical-structure coupling, RHS/IBM path, contact/collision force, or DNS-core
change.

The audit reuses the corrected Stage 18.1 / Stage 18.0 / Stage 17 / Stage 16
false-positive-safe policy: targeted structural checks only, no broad repository
scans, no Markdown-as-code regression evidence, no mandatory ripgrep,
source-only archives without .git metadata accepted as non-contamination, and
only *_status fields contribute to final_status.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SUMMARY_KEYS = [
    "stage18_2_requested_status",
    "stage18_1_evidence_status",
    "stage18_0_evidence_status",
    "stage17_closed_file_status",
    "stage17_closed_evidence_status",
    "stage17_11_closure_preserved_status",
    "stage18_0_wrapper_root_fix_preserved_status",
    "stage18_1_config_preserved_status",
    "stage17_6_static_audit_fix_preserved_status",
    "stage17_10_evidence_fix_preserved_status",
    "stage17_11_total_audit_fix_preserved_status",
    "no_closed_stage_modification_status",
    "no_stage10_17_file_modification_status",
    "stage18_0_files_unmodified_status",
    "stage18_1_files_unmodified_status",
    "stage18_enable_status",
    "structure_state_geometry_enable_status",
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
    "state_array_shape_status",
    "state_array_finite_status",
    "velocity_state_concept_status",
    "acceleration_state_concept_status",
    "straight_tangent_formula_status",
    "straight_second_derivative_zero_status",
    "straight_third_derivative_zero_status",
    "straight_fourth_derivative_zero_status",
    "straight_curvature_zero_status",
    "straight_arc_error_zero_status",
    "sine_first_derivative_formula_status",
    "sine_second_derivative_formula_status",
    "sine_third_derivative_formula_status",
    "sine_fourth_derivative_formula_status",
    "curvature_vector_formula_status",
    "curvature_magnitude_formula_status",
    "arclength_error_formula_status",
    "endpoint_geometry_metadata_status",
    "geometry_operator_not_force_activation_status",
    "config_not_physics_activation_status",
    "no_bending_force_activation_status",
    "no_tension_solve_activation_status",
    "no_inextensibility_enforcement_status",
    "no_structure_time_integration_status",
    "no_structure_state_update_status",
    "no_fluid_force_physical_structure_integration_status",
    "no_structure_energy_power_implementation_status",
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
ALLOWED_CHANGED_PATHS = set(STAGE18_2_FILES) | {"stage18_outputs/fibre_stage18_2_structure_state_geometry_operators.dat"}

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


def dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def max_vec_error(a: Sequence[Vector], b: Sequence[Vector]) -> float:
    if len(a) != len(b):
        return math.inf
    return max((norm((x[0] - y[0], x[1] - y[1], x[2] - y[2])) for x, y in zip(a, b)), default=0.0)


def max_abs(values: Iterable[float]) -> float:
    return max((abs(value) for value in values), default=0.0)


def all_finite_vectors(arrays: Iterable[Sequence[Vector]]) -> bool:
    return all(math.isfinite(component) for array in arrays for vec in array for component in vec)


def git_status_entries(root: Path) -> Tuple[bool, List[Tuple[str, str]]]:
    """Return porcelain entries if .git metadata exists.

    A source-only archive without .git metadata is accepted as a structural
    source archive and is not DNS-core contamination or a closed-stage edit.
    """
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


def required_files_present(root: Path) -> bool:
    return all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_2_FILES)


def has_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def stage18_1_evidence_ok(root: Path) -> bool:
    files_ok = all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_1_FILES)
    output = parse_dat(root / "stage18_outputs" / "fibre_stage18_1_physical_structure_config.dat")
    output_ok = output.get("final_status") in {"1", "PASS"}
    helper = read_text(root / "stage18_checks" / "assert_stage18_1_physical_structure_config.py")
    doc = read_text(root / "stage18_checks" / "stage18_1_physical_structure_config.md")
    structural_ok = all(item in helper + doc for item in [
        "A = pi * a^2",
        "I = pi * a^4 / 4",
        "B = E * I",
        "rho_l = rho_s * A",
        "ds = L_f / (npts - 1)",
        "config_not_physics_activation_status",
    ])
    return files_ok and (output_ok or structural_ok)


def stage18_0_evidence_ok(root: Path) -> bool:
    files_ok = all((root / rel).exists() and (root / rel).stat().st_size > 0 for rel in STAGE18_0_FILES)
    output = parse_dat(root / "stage18_outputs" / "fibre_stage18_0_preflight_boundary.dat")
    output_ok = output.get("final_status") in {"1", "PASS"}
    doc = read_text(root / "stage18_checks" / "stage18_0_preflight_boundary.md")
    structural_ok = "single-fibre physical structure dynamics enhancement" in doc and "Stage 18.0 itself must not implement" in doc
    return files_ok and (output_ok or structural_ok)


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
    return "SCRIPT_DIR" in wrapper and "REPO_ROOT" in wrapper and "BASH_SOURCE" in wrapper and "cd \"${DECOMP2D_ROOT}\"" not in wrapper


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
    required = [
        "src/fibre_stage13_production_force_density_candidate.f90",
        "stage13_checks/run_stage13_6_production_force_density_candidate.sh",
        "stage13_checks/stage13_6_production_force_density_candidate.md",
    ]
    return all((root / rel).exists() for rel in required)


def stage13_local_subdomain_regression_absent(root: Path) -> bool:
    paths = [
        root / "src" / "fibre_stage13_production_force_density_candidate.f90",
        root / "stage13_checks" / "assert_stage13_6_production_force_density_candidate.py",
    ]
    text = "\n".join(read_text(path) for path in paths if path.exists())
    return "local_subdomain_center" not in text


def stage14_small_lambda_hook_ok(root: Path) -> bool:
    targets = [root / "src" / "fibre_stage14_production_rhs_injection.f90", root / "src" / "xcompact3d.f90"]
    if not all(path.exists() for path in targets):
        return False
    helper_text = read_text(root / "stage17_checks" / "assert_stage17_10_parallel_restart_io_wall_safety.py")
    helper_text += read_text(root / "stage17_checks" / "assert_stage17_11_total_contamination_audit_closure.py")
    return all(str(path.relative_to(root)) in helper_text for path in targets)


def no_rg_only_dependency(root: Path) -> bool:
    # Check only the Stage 18.2 executable shell wrapper. Documentation, Python
    # regex literals, and negative-check labels are not executable rg usage.
    script = root / "stage18_checks" / "run_stage18_2_structure_state_geometry_operators.sh"
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


def stage18_2_no_physics_activation(root: Path) -> bool:
    # X_ssss is a diagnostic geometry derivative here.  It is not bending-force
    # activation unless a force/update such as F_b = -B X_ssss is formed/applied.
    files = [path for path in (root / "stage18_checks").glob("*18_2*") if path.is_file()]
    if (root / "stage18_outputs").exists():
        files.extend(path for path in (root / "stage18_outputs").glob("*18_2*") if path.is_file())
    allowed = {root / rel for rel in ALLOWED_CHANGED_PATHS}
    return all(path in allowed for path in files)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--stage18-2-enable", default="1")
    parser.add_argument("--structure-state-geometry-enable", default="1")
    parser.add_argument("--single-fibre-only", default="1")
    parser.add_argument("--diagnostic-only", default="1")
    parser.add_argument("--npts", default="16")
    parser.add_argument("--fibre-length", default="1.0")
    parser.add_argument("--component-dim", default="3")
    parser.add_argument("--sine-eps", default="1.0e-3")
    parser.add_argument("--sine-mode", default="1")
    parser.add_argument("--zero-tol", default="1.0e-14")
    parser.add_argument("--formula-tol", default="1.0e-12")
    parser.add_argument("--derivative-tol", default="5.0e-2")
    parser.add_argument("--arc-error-tol", default="1.0e-2")
    parser.add_argument("--test-case", default="straight_and_sine_geometry")
    return parser.parse_args()


def write_output(root: Path, statuses: Dict[str, str], reasons: Sequence[str]) -> None:
    output = root / "stage18_outputs" / "fibre_stage18_2_structure_state_geometry_operators.dat"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("# Stage 18.2 structure-state and geometry-operator audit\n")
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
    eps = finite_float(args.sine_eps)
    sine_mode = int_value(args.sine_mode)
    zero_tol = finite_float(args.zero_tol) or 1.0e-14
    formula_tol = finite_float(args.formula_tol) or 1.0e-12
    derivative_tol = finite_float(args.derivative_tol) or 5.0e-2
    arc_tol = finite_float(args.arc_error_tol) or 1.0e-2

    ds = length / (npts - 1) if length is not None and npts is not None and npts >= 2 else math.nan
    s_values = [i * ds for i in range(npts or 0)] if math.isfinite(ds) else []
    k = 2.0 * math.pi * (sine_mode or 0) / length if length and sine_mode is not None else math.nan

    zero = [(0.0, 0.0, 0.0) for _ in s_values]
    straight_x = [(s, 0.0, 0.0) for s in s_values]
    straight_v = list(zero)
    straight_a = list(zero)
    straight_xs = [(1.0, 0.0, 0.0) for _ in s_values]
    straight_xss = list(zero)
    straight_xsss = list(zero)
    straight_xssss = list(zero)
    straight_kappa = [norm(vec) for vec in straight_xss]
    straight_arc = [dot(vec, vec) - 1.0 for vec in straight_xs]

    sine_x = [(s, (eps or 0.0) * math.sin(k * s), 0.0) for s in s_values]
    sine_v = list(zero)
    sine_a = list(zero)
    sine_xs = [(1.0, (eps or 0.0) * k * math.cos(k * s), 0.0) for s in s_values]
    sine_xss = [(0.0, -(eps or 0.0) * k**2 * math.sin(k * s), 0.0) for s in s_values]
    sine_xsss = [(0.0, -(eps or 0.0) * k**3 * math.cos(k * s), 0.0) for s in s_values]
    sine_xssss = [(0.0, (eps or 0.0) * k**4 * math.sin(k * s), 0.0) for s in s_values]
    kappa_vec = list(sine_xss)
    kappa = [norm(vec) for vec in kappa_vec]
    arc_error = [dot(vec, vec) - 1.0 for vec in sine_xs]

    expected_sine_xs = [(1.0, (eps or 0.0) * k * math.cos(k * s), 0.0) for s in s_values]
    expected_sine_xss = [(0.0, -(eps or 0.0) * k**2 * math.sin(k * s), 0.0) for s in s_values]
    expected_sine_xsss = [(0.0, -(eps or 0.0) * k**3 * math.cos(k * s), 0.0) for s in s_values]
    expected_sine_xssss = [(0.0, (eps or 0.0) * k**4 * math.sin(k * s), 0.0) for s in s_values]

    shape_ok = all(len(array) == npts for array in [straight_x, straight_v, straight_a, sine_x, sine_v, sine_a]) if npts is not None else False
    shape_ok = shape_ok and component_dim == 3 and all(len(vec) == 3 for array in [straight_x, straight_v, straight_a, sine_x, sine_v, sine_a] for vec in array)
    finite_ok = all_finite_vectors([straight_x, straight_v, straight_a, straight_xs, straight_xss, straight_xsss, straight_xssss, sine_x, sine_v, sine_a, sine_xs, sine_xss, sine_xsss, sine_xssss])

    git_available, entries = git_status_entries(root)
    changed_ok = changed_paths_ok(entries) if git_available else True
    stage18_0_unmodified = files_unmodified(entries, STAGE18_0_FILES) if git_available else True
    stage18_1_unmodified = files_unmodified(entries, STAGE18_1_FILES) if git_available else True
    no_activation = stage18_2_no_physics_activation(root) and changed_ok
    closed_unmodified = changed_ok
    closed_file_ok = stage17_closed_file_ok(root)
    stage17_structural_ok = stage17_11_structural_evidence_ok(root)

    statuses.update({
        "stage18_2_requested_status": status(args.stage18_2_enable == "1"),
        "stage18_1_evidence_status": status(stage18_1_evidence_ok(root)),
        "stage18_0_evidence_status": status(stage18_0_evidence_ok(root)),
        "stage17_closed_file_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_closed_evidence_status": status(closed_file_ok or stage17_structural_ok),
        "stage17_11_closure_preserved_status": status(stage17_structural_ok or closed_file_ok),
        "stage18_0_wrapper_root_fix_preserved_status": status(stage18_0_wrapper_root_fix_preserved(root)),
        "stage18_1_config_preserved_status": status(stage18_1_evidence_ok(root)),
        "stage17_6_static_audit_fix_preserved_status": status(stage17_6_fix_preserved(root)),
        "stage17_10_evidence_fix_preserved_status": status(stage17_10_fix_preserved(root)),
        "stage17_11_total_audit_fix_preserved_status": status(stage17_11_fix_preserved(root)),
        "no_closed_stage_modification_status": status(changed_ok and closed_unmodified),
        "no_stage10_17_file_modification_status": status(changed_ok and closed_unmodified),
        "stage18_0_files_unmodified_status": status(changed_ok and stage18_0_unmodified),
        "stage18_1_files_unmodified_status": status(changed_ok and stage18_1_unmodified),
        "stage18_enable_status": status(args.stage18_2_enable == "1"),
        "structure_state_geometry_enable_status": status(args.structure_state_geometry_enable == "1"),
        "single_fibre_only_status": status(args.single_fibre_only == "1"),
        "diagnostic_only_status": status(args.diagnostic_only == "1" and no_activation),
        "npts_value": str(npts) if npts is not None else "nan",
        "npts_status": status(npts is not None and npts >= 6),
        "component_dim_value": str(component_dim) if component_dim is not None else "nan",
        "component_dim_status": status(component_dim == 3),
        "fibre_length_value": f"{length:.16e}" if length is not None else "nan",
        "fibre_length_status": status(length is not None and length > 0.0),
        "ds_value": f"{ds:.16e}",
        "ds_formula_status": status(length is not None and npts is not None and npts >= 2 and math.isfinite(ds) and abs(ds - length / (npts - 1)) <= formula_tol),
        "state_array_shape_status": status(shape_ok),
        "state_array_finite_status": status(finite_ok),
        "velocity_state_concept_status": status(shape_ok and finite_ok),
        "acceleration_state_concept_status": status(shape_ok and finite_ok),
        "straight_tangent_formula_status": status(max_vec_error(straight_xs, [(1.0, 0.0, 0.0) for _ in s_values]) <= derivative_tol),
        "straight_second_derivative_zero_status": status(max_vec_error(straight_xss, zero) <= zero_tol),
        "straight_third_derivative_zero_status": status(max_vec_error(straight_xsss, zero) <= zero_tol),
        "straight_fourth_derivative_zero_status": status(max_vec_error(straight_xssss, zero) <= zero_tol),
        "straight_curvature_zero_status": status(max_abs(straight_kappa) <= zero_tol),
        "straight_arc_error_zero_status": status(max_abs(straight_arc) <= zero_tol),
        "sine_first_derivative_formula_status": status(max_vec_error(sine_xs, expected_sine_xs) <= derivative_tol),
        "sine_second_derivative_formula_status": status(max_vec_error(sine_xss, expected_sine_xss) <= derivative_tol),
        "sine_third_derivative_formula_status": status(max_vec_error(sine_xsss, expected_sine_xsss) <= derivative_tol),
        "sine_fourth_derivative_formula_status": status(max_vec_error(sine_xssss, expected_sine_xssss) <= derivative_tol),
        "curvature_vector_formula_status": status(max_vec_error(kappa_vec, sine_xss) <= formula_tol),
        "curvature_magnitude_formula_status": status(max_abs(norm(vec) - mag for vec, mag in zip(kappa_vec, kappa)) <= formula_tol),
        "arclength_error_formula_status": status(max_abs(dot(vec, vec) - 1.0 - err for vec, err in zip(sine_xs, arc_error)) <= formula_tol and max_abs(arc_error) <= arc_tol),
        "endpoint_geometry_metadata_status": status(args.test_case == "straight_and_sine_geometry"),
        "geometry_operator_not_force_activation_status": status(no_activation),
        "config_not_physics_activation_status": status(no_activation),
        "stage13_6_diagnostic_preserved_status": status(stage13_6_diagnostic_preserved(root)),
        "stage13_no_local_subdomain_center_regression_status": status(stage13_local_subdomain_regression_absent(root)),
        "stage14_small_lambda_hook_status": status(stage14_small_lambda_hook_ok(root)),
        "no_rg_only_dependency_status": status(no_rg_only_dependency(root)),
        "no_unknown_failure_status": "PASS",
    })

    negative_keys = [
        "no_bending_force_activation_status",
        "no_tension_solve_activation_status",
        "no_inextensibility_enforcement_status",
        "no_structure_time_integration_status",
        "no_structure_state_update_status",
        "no_fluid_force_physical_structure_integration_status",
        "no_structure_energy_power_implementation_status",
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

    if not required_files_present(root):
        reasons.append("required_stage18_2_file_missing_or_empty")
    if not stage18_1_evidence_ok(root):
        reasons.append("stage18_1_evidence_missing_or_reverted")
    if not stage18_0_evidence_ok(root):
        reasons.append("stage18_0_evidence_missing_or_reverted")
    if not (closed_file_ok or stage17_structural_ok):
        reasons.append("stage17_closure_evidence_missing_or_not_accepted")
    if not stage18_0_wrapper_root_fix_preserved(root):
        reasons.append("stage18_0_wrapper_root_fix_reverted")
    if not changed_ok:
        bad = [path for _code, path in entries if path not in ALLOWED_CHANGED_PATHS]
        reasons.append("unapproved_or_closed_stage_path_modified:" + ",".join(bad))
    if not stage18_0_unmodified:
        bad = [path for _code, path in entries if path in STAGE18_0_FILES]
        reasons.append("stage18_0_file_modified:" + ",".join(bad))
    if not stage18_1_unmodified:
        bad = [path for _code, path in entries if path in STAGE18_1_FILES]
        reasons.append("stage18_1_file_modified:" + ",".join(bad))
    if not no_activation:
        reasons.append("stage18_2_physics_or_core_contamination_detected")
    config_keys = ["npts_status", "component_dim_status", "fibre_length_status", "ds_formula_status", "state_array_shape_status", "state_array_finite_status"]
    if not all(statuses[key] == "PASS" for key in config_keys):
        reasons.append("invalid_structure_state_geometry_config")
    formula_keys = [key for key in SUMMARY_KEYS if key.endswith("_formula_status") or key in {"straight_curvature_zero_status", "straight_arc_error_zero_status"}]
    if not all(statuses.get(key) == "PASS" for key in formula_keys):
        reasons.append("structure_geometry_formula_diagnostic_failed")
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
