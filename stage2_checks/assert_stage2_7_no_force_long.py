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
        print(f"STAGE 2.7 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_no_force_long_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    rest_long_preservation_error_max = get_float(data, "rest_long_preservation_error_max")
    rest_long_final_bending_energy = get_float(data, "rest_long_final_bending_energy")
    rest_long_final_kinetic_energy = get_float(data, "rest_long_final_kinetic_energy")
    rest_long_final_max_curvature = get_float(data, "rest_long_final_max_curvature")
    rest_long_final_length_error = get_float(data, "rest_long_final_length_error")
    rest_long_momentum_norm = get_float(data, "rest_long_momentum_norm")
    rest_long_center_drift_norm = get_float(data, "rest_long_center_drift_norm")
    rest_long_solver_failure_count = get_int(data, "rest_long_solver_failure_count")

    translation_long_center_error_norm = get_float(data, "translation_long_center_error_norm")
    translation_long_velocity_error_norm = get_float(data, "translation_long_velocity_error_norm")
    translation_long_shape_error_max = get_float(data, "translation_long_shape_error_max")
    translation_long_bending_energy = get_float(data, "translation_long_bending_energy")
    translation_long_length_error = get_float(data, "translation_long_length_error")
    translation_long_momentum_error_norm = get_float(data, "translation_long_momentum_error_norm")
    translation_long_solver_failure_count = get_int(data, "translation_long_solver_failure_count")

    curved_long_initial_bending_energy = get_float(data, "curved_long_initial_bending_energy")
    curved_long_final_bending_energy = get_float(data, "curved_long_final_bending_energy")
    curved_long_max_total_energy = get_float(data, "curved_long_max_total_energy")
    curved_long_initial_length_error = get_float(data, "curved_long_initial_length_error")
    curved_long_final_length_error = get_float(data, "curved_long_final_length_error")
    curved_long_max_length_error = get_float(data, "curved_long_max_length_error")
    curved_long_final_momentum_norm = get_float(data, "curved_long_final_momentum_norm")
    curved_long_final_momentum_relative = get_float(data, "curved_long_final_momentum_relative")
    curved_long_center_drift_norm = get_float(data, "curved_long_center_drift_norm")
    curved_long_final_max_curvature = get_float(data, "curved_long_final_max_curvature")
    curved_long_nan_detected = get_int(data, "curved_long_nan_detected")
    curved_long_solver_failure_count = get_int(data, "curved_long_solver_failure_count")
    curved_long_endpoint_fixed_constraint_detected = get_int(data, "curved_long_endpoint_fixed_constraint_detected")

    require(rest_long_preservation_error_max <= 1e-9,
            f"rest_long_preservation_error_max={rest_long_preservation_error_max}")
    require(rest_long_final_bending_energy <= 1e-10,
            f"rest_long_final_bending_energy={rest_long_final_bending_energy}")
    require(rest_long_final_kinetic_energy <= 1e-10,
            f"rest_long_final_kinetic_energy={rest_long_final_kinetic_energy}")
    require(rest_long_final_max_curvature <= 1e-8,
            f"rest_long_final_max_curvature={rest_long_final_max_curvature}")
    require(rest_long_final_length_error <= 1e-10,
            f"rest_long_final_length_error={rest_long_final_length_error}")
    require(rest_long_momentum_norm <= 1e-8,
            f"rest_long_momentum_norm={rest_long_momentum_norm}")
    require(rest_long_center_drift_norm <= 1e-9,
            f"rest_long_center_drift_norm={rest_long_center_drift_norm}")
    require(rest_long_solver_failure_count == 0,
            f"rest_long_solver_failure_count={rest_long_solver_failure_count}")

    require(translation_long_center_error_norm <= 1e-8,
            f"translation_long_center_error_norm={translation_long_center_error_norm}")
    require(translation_long_velocity_error_norm <= 1e-8,
            f"translation_long_velocity_error_norm={translation_long_velocity_error_norm}")
    require(translation_long_shape_error_max <= 1e-8,
            f"translation_long_shape_error_max={translation_long_shape_error_max}")
    require(translation_long_bending_energy <= 1e-10,
            f"translation_long_bending_energy={translation_long_bending_energy}")
    require(translation_long_length_error <= 1e-10,
            f"translation_long_length_error={translation_long_length_error}")
    require(translation_long_momentum_error_norm <= 1e-8,
            f"translation_long_momentum_error_norm={translation_long_momentum_error_norm}")
    require(translation_long_solver_failure_count == 0,
            f"translation_long_solver_failure_count={translation_long_solver_failure_count}")

    require(curved_long_initial_bending_energy > 0.0,
            f"curved_long_initial_bending_energy={curved_long_initial_bending_energy}")
    require(curved_long_final_bending_energy >= 0.0,
            f"curved_long_final_bending_energy={curved_long_final_bending_energy}")
    require(curved_long_max_total_energy <= 10.0 * curved_long_initial_bending_energy,
            f"curved_long_max_total_energy={curved_long_max_total_energy}, initial={curved_long_initial_bending_energy}")
    require(curved_long_initial_length_error <= 1e-12,
            f"curved_long_initial_length_error={curved_long_initial_length_error}")
    require(curved_long_final_length_error <= 1e-6,
            f"curved_long_final_length_error={curved_long_final_length_error}")
    require(curved_long_max_length_error <= 1e-6,
            f"curved_long_max_length_error={curved_long_max_length_error}")
    require(curved_long_final_momentum_norm <= 1e-4,
            f"curved_long_final_momentum_norm={curved_long_final_momentum_norm}")
    require(curved_long_final_momentum_relative <= 1e-6,
            f"curved_long_final_momentum_relative={curved_long_final_momentum_relative}")
    require(curved_long_center_drift_norm <= 1e-6,
            f"curved_long_center_drift_norm={curved_long_center_drift_norm}")
    require(curved_long_final_max_curvature >= 0.0,
            f"curved_long_final_max_curvature={curved_long_final_max_curvature}")
    require(curved_long_nan_detected == 0,
            f"curved_long_nan_detected={curved_long_nan_detected}")
    require(curved_long_solver_failure_count == 0,
            f"curved_long_solver_failure_count={curved_long_solver_failure_count}")
    require(curved_long_endpoint_fixed_constraint_detected == 0,
            f"curved_long_endpoint_fixed_constraint_detected={curved_long_endpoint_fixed_constraint_detected}")

    print("STAGE 2.7 NO-FORCE LONG-RUN CHECK PASSED")


if __name__ == "__main__":
    main()
