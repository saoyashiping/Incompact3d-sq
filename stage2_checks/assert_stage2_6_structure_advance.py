#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys


def parse_kv(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip()] = value.strip()
    return data


def get_float(data: dict[str, str], key: str) -> float:
    if key not in data:
        raise KeyError(key)
    return float(data[key].split()[0])


def get_int(data: dict[str, str], key: str) -> int:
    if key not in data:
        raise KeyError(key)
    return int(data[key].split()[0])


def require(cond: bool, msg: str) -> None:
    if not cond:
        print(f"STAGE 2.6 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_structure_advance_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    rest_preservation_error_max = get_float(data, "rest_preservation_error_max")
    rest_final_bending_energy = get_float(data, "rest_final_bending_energy")
    rest_final_kinetic_energy = get_float(data, "rest_final_kinetic_energy")
    rest_final_max_curvature = get_float(data, "rest_final_max_curvature")
    rest_final_length_error = get_float(data, "rest_final_length_error")
    rest_solver_failure_count = get_int(data, "rest_solver_failure_count")

    translation_center_error_norm = get_float(data, "translation_center_error_norm")
    translation_velocity_error_norm = get_float(data, "translation_velocity_error_norm")
    translation_shape_error_max = get_float(data, "translation_shape_error_max")
    translation_final_bending_energy = get_float(data, "translation_final_bending_energy")
    translation_final_length_error = get_float(data, "translation_final_length_error")
    translation_momentum_error_norm = get_float(data, "translation_momentum_error_norm")
    translation_solver_failure_count = get_int(data, "translation_solver_failure_count")

    stretch_initial_length_error = get_float(data, "stretch_initial_length_error")
    stretch_final_length_error = get_float(data, "stretch_final_length_error")
    stretch_tension_norm = get_float(data, "stretch_tension_norm")
    stretch_solver_failure_count = get_int(data, "stretch_solver_failure_count")

    curved_initial_bending_energy = get_float(data, "curved_initial_bending_energy")
    curved_final_bending_energy = get_float(data, "curved_final_bending_energy")
    curved_final_max_curvature = get_float(data, "curved_final_max_curvature")
    curved_nan_detected = get_int(data, "curved_nan_detected")
    curved_solver_failure_count = get_int(data, "curved_solver_failure_count")
    curved_endpoint_fixed_constraint_detected = get_int(data, "curved_endpoint_fixed_constraint_detected")

    advance_matrix_apply_maxdiff = get_float(data, "advance_matrix_apply_maxdiff")

    require(rest_preservation_error_max <= 1e-10, f"rest_preservation_error_max={rest_preservation_error_max}")
    require(rest_final_bending_energy <= 1e-12, f"rest_final_bending_energy={rest_final_bending_energy}")
    require(rest_final_kinetic_energy <= 1e-12, f"rest_final_kinetic_energy={rest_final_kinetic_energy}")
    require(rest_final_max_curvature <= 1e-10, f"rest_final_max_curvature={rest_final_max_curvature}")
    require(rest_final_length_error <= 1e-12, f"rest_final_length_error={rest_final_length_error}")
    require(rest_solver_failure_count == 0, f"rest_solver_failure_count={rest_solver_failure_count}")

    require(translation_center_error_norm <= 1e-10,
            f"translation_center_error_norm={translation_center_error_norm}")
    require(translation_velocity_error_norm <= 1e-10,
            f"translation_velocity_error_norm={translation_velocity_error_norm}")
    require(translation_shape_error_max <= 1e-10,
            f"translation_shape_error_max={translation_shape_error_max}")
    require(translation_final_bending_energy <= 1e-12,
            f"translation_final_bending_energy={translation_final_bending_energy}")
    require(translation_final_length_error <= 1e-12,
            f"translation_final_length_error={translation_final_length_error}")
    require(translation_momentum_error_norm <= 1e-10,
            f"translation_momentum_error_norm={translation_momentum_error_norm}")
    require(translation_solver_failure_count == 0,
            f"translation_solver_failure_count={translation_solver_failure_count}")

    require(stretch_initial_length_error > 0.0,
            f"stretch_initial_length_error={stretch_initial_length_error}")
    require(stretch_tension_norm > 0.0, f"stretch_tension_norm={stretch_tension_norm}")
    require(stretch_final_length_error < stretch_initial_length_error,
            f"stretch_final_length_error={stretch_final_length_error}, initial={stretch_initial_length_error}")
    require(stretch_solver_failure_count == 0,
            f"stretch_solver_failure_count={stretch_solver_failure_count}")

    require(curved_initial_bending_energy > 0.0,
            f"curved_initial_bending_energy={curved_initial_bending_energy}")
    require(curved_final_bending_energy >= 0.0,
            f"curved_final_bending_energy={curved_final_bending_energy}")
    require(curved_final_bending_energy <= 10.0 * curved_initial_bending_energy,
            f"curved_final_bending_energy={curved_final_bending_energy}, initial={curved_initial_bending_energy}")
    require(curved_final_max_curvature >= 0.0,
            f"curved_final_max_curvature={curved_final_max_curvature}")
    require(curved_nan_detected == 0, f"curved_nan_detected={curved_nan_detected}")
    require(curved_solver_failure_count == 0,
            f"curved_solver_failure_count={curved_solver_failure_count}")
    require(curved_endpoint_fixed_constraint_detected == 0,
            f"curved_endpoint_fixed_constraint_detected={curved_endpoint_fixed_constraint_detected}")

    require(advance_matrix_apply_maxdiff <= 1e-10,
            f"advance_matrix_apply_maxdiff={advance_matrix_apply_maxdiff}")

    print("STAGE 2.6 STRUCTURE ADVANCE CHECK PASSED")


if __name__ == "__main__":
    main()
