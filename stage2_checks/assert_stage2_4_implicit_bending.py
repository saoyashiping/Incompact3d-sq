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
        print(f"STAGE 2.4 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_implicit_bending_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    d4_diff = get_float(data, "d4_matrix_operator_maxdiff")

    straight_preservation_error_max = get_float(data, "straight_preservation_error_max")
    straight_final_bending_energy = get_float(data, "straight_final_bending_energy")
    straight_final_max_curvature = get_float(data, "straight_final_max_curvature")

    sine_initial_bending_energy = get_float(data, "sine_initial_bending_energy")
    sine_final_bending_energy = get_float(data, "sine_final_bending_energy")
    sine_energy_increase_count = get_int(data, "sine_energy_increase_count")
    sine_max_energy_increase = get_float(data, "sine_max_energy_increase")
    sine_initial_max_curvature = get_float(data, "sine_initial_max_curvature")
    sine_final_max_curvature = get_float(data, "sine_final_max_curvature")
    sine_final_freefree_boundary_residual = get_float(data, "sine_final_freefree_boundary_residual")

    large_dt_initial_bending_energy = get_float(data, "large_dt_initial_bending_energy")
    large_dt_final_bending_energy = get_float(data, "large_dt_final_bending_energy")
    large_dt_energy_increase_count = get_int(data, "large_dt_energy_increase_count")
    large_dt_solver_failure_count = get_int(data, "large_dt_solver_failure_count")

    endpoint_fixed_constraint_detected = get_int(data, "endpoint_fixed_constraint_detected")

    require(d4_diff <= 1e-10, f"d4_matrix_operator_maxdiff={d4_diff}")

    require(straight_preservation_error_max <= 1e-12,
            f"straight_preservation_error_max={straight_preservation_error_max}")
    require(straight_final_bending_energy <= 1e-14,
            f"straight_final_bending_energy={straight_final_bending_energy}")
    require(straight_final_max_curvature <= 1e-14,
            f"straight_final_max_curvature={straight_final_max_curvature}")

    require(sine_initial_bending_energy > 0.0,
            f"sine_initial_bending_energy={sine_initial_bending_energy}")
    require(sine_final_bending_energy < sine_initial_bending_energy,
            f"sine_final_bending_energy={sine_final_bending_energy}, initial={sine_initial_bending_energy}")
    require(sine_energy_increase_count == 0,
            f"sine_energy_increase_count={sine_energy_increase_count}")
    require(sine_max_energy_increase <= 1e-12,
            f"sine_max_energy_increase={sine_max_energy_increase}")
    require(sine_final_max_curvature < sine_initial_max_curvature,
            f"sine_final_max_curvature={sine_final_max_curvature}, initial={sine_initial_max_curvature}")
    require(sine_final_freefree_boundary_residual <= 1e-10,
            f"sine_final_freefree_boundary_residual={sine_final_freefree_boundary_residual}")

    require(large_dt_final_bending_energy < large_dt_initial_bending_energy,
            f"large_dt_final_bending_energy={large_dt_final_bending_energy}, initial={large_dt_initial_bending_energy}")
    require(large_dt_energy_increase_count == 0,
            f"large_dt_energy_increase_count={large_dt_energy_increase_count}")
    require(large_dt_solver_failure_count == 0,
            f"large_dt_solver_failure_count={large_dt_solver_failure_count}")

    require(endpoint_fixed_constraint_detected == 0,
            f"endpoint_fixed_constraint_detected={endpoint_fixed_constraint_detected}")

    print("STAGE 2.4 IMPLICIT BENDING CHECK PASSED")


if __name__ == "__main__":
    main()
