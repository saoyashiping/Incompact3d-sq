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


def require(cond: bool, msg: str) -> None:
    if not cond:
        print(f"STAGE 2.3 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_energy_diagnostics_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    straight_rest_bending_energy = get_float(data, "straight_rest_bending_energy")
    straight_rest_kinetic_energy = get_float(data, "straight_rest_kinetic_energy")
    straight_rest_total_energy = get_float(data, "straight_rest_total_energy")
    straight_rest_linear_momentum_norm = get_float(data, "straight_rest_linear_momentum_norm")
    straight_rest_max_curvature = get_float(data, "straight_rest_max_curvature")
    straight_rest_stretch_ratio = get_float(data, "straight_rest_stretch_ratio")
    straight_rest_total_length_relative_error = get_float(data, "straight_rest_total_length_relative_error")

    translation_bending_energy = get_float(data, "translation_bending_energy")
    translation_kinetic_energy = get_float(data, "translation_kinetic_energy")
    translation_total_energy = get_float(data, "translation_total_energy")
    translation_momentum_consistency_error = get_float(data, "translation_momentum_consistency_error")
    translation_relative_kinetic_energy = get_float(data, "translation_relative_kinetic_energy")
    translation_kinetic_split_error = get_float(data, "translation_kinetic_split_error")

    sine_bending_energy = get_float(data, "sine_bending_energy")
    sine_kinetic_energy = get_float(data, "sine_kinetic_energy")
    sine_total_energy = get_float(data, "sine_total_energy")
    sine_max_curvature = get_float(data, "sine_max_curvature")
    sine_rms_curvature = get_float(data, "sine_rms_curvature")
    sine_total_length_relative_error = get_float(data, "sine_total_length_relative_error")

    require(straight_rest_bending_energy <= 1e-14, f"straight_rest_bending_energy={straight_rest_bending_energy}")
    require(straight_rest_kinetic_energy <= 1e-14, f"straight_rest_kinetic_energy={straight_rest_kinetic_energy}")
    require(straight_rest_total_energy <= 1e-14, f"straight_rest_total_energy={straight_rest_total_energy}")
    require(straight_rest_linear_momentum_norm <= 1e-14,
            f"straight_rest_linear_momentum_norm={straight_rest_linear_momentum_norm}")
    require(straight_rest_max_curvature <= 1e-14, f"straight_rest_max_curvature={straight_rest_max_curvature}")
    require(abs(straight_rest_stretch_ratio - 1.0) <= 1e-14,
            f"straight_rest_stretch_ratio={straight_rest_stretch_ratio}")
    require(straight_rest_total_length_relative_error <= 1e-14,
            f"straight_rest_total_length_relative_error={straight_rest_total_length_relative_error}")

    require(translation_bending_energy <= 1e-14, f"translation_bending_energy={translation_bending_energy}")
    require(translation_kinetic_energy > 0.0, f"translation_kinetic_energy={translation_kinetic_energy}")
    require(translation_total_energy > 0.0, f"translation_total_energy={translation_total_energy}")
    require(translation_momentum_consistency_error <= 1e-14,
            f"translation_momentum_consistency_error={translation_momentum_consistency_error}")
    require(translation_relative_kinetic_energy <= 1e-14,
            f"translation_relative_kinetic_energy={translation_relative_kinetic_energy}")
    require(translation_kinetic_split_error <= 1e-14,
            f"translation_kinetic_split_error={translation_kinetic_split_error}")

    require(sine_bending_energy > 0.0, f"sine_bending_energy={sine_bending_energy}")
    require(sine_kinetic_energy <= 1e-14, f"sine_kinetic_energy={sine_kinetic_energy}")
    require(abs(sine_total_energy - sine_bending_energy) <= max(1e-14, 1e-12 * abs(sine_bending_energy)),
            f"sine_total_energy={sine_total_energy}, sine_bending_energy={sine_bending_energy}")
    require(sine_max_curvature > 0.0, f"sine_max_curvature={sine_max_curvature}")
    require(sine_rms_curvature > 0.0, f"sine_rms_curvature={sine_rms_curvature}")
    require(sine_total_length_relative_error > 0.0,
            f"sine_total_length_relative_error={sine_total_length_relative_error}")

    print("STAGE 2.3 ENERGY DIAGNOSTICS CHECK PASSED")


if __name__ == "__main__":
    main()
