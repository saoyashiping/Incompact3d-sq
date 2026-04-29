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
        print(f"STAGE 2.5 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_tension_solver_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    straight_rhs_norm = get_float(data, "straight_rhs_norm")
    straight_tension_norm = get_float(data, "straight_tension_norm")
    straight_residual_norm = get_float(data, "straight_residual_norm")

    translation_rhs_norm = get_float(data, "translation_rhs_norm")
    translation_tension_norm = get_float(data, "translation_tension_norm")
    translation_residual_norm = get_float(data, "translation_residual_norm")

    stretch_rhs_norm = get_float(data, "stretch_rhs_norm")
    stretch_tension_norm = get_float(data, "stretch_tension_norm")
    stretch_residual_norm = get_float(data, "stretch_residual_norm")
    stretch_relative_residual = get_float(data, "stretch_relative_residual")

    curved_rhs_norm = get_float(data, "curved_rhs_norm")
    curved_residual_norm = get_float(data, "curved_residual_norm")
    curved_relative_residual = get_float(data, "curved_relative_residual")
    curved_endpoint_change_norm = get_float(data, "curved_endpoint_change_norm")

    tension_matrix_apply_maxdiff = get_float(data, "tension_matrix_apply_maxdiff")

    require(straight_rhs_norm <= 1e-12, f"straight_rhs_norm={straight_rhs_norm}")
    require(straight_tension_norm <= 1e-10, f"straight_tension_norm={straight_tension_norm}")
    require(straight_residual_norm <= 1e-12, f"straight_residual_norm={straight_residual_norm}")

    # The translation RHS is theoretically zero, but the Huang constraint-correction term contains dt^{-2};
    # roundoff-level differences in |t_n|^2 and |t_old|^2 can be amplified to O(1e-10).
    # This is acceptable if solved tension and residual remain small.
    require(translation_rhs_norm <= 1e-8, f"translation_rhs_norm={translation_rhs_norm}")
    require(translation_tension_norm <= 1e-9, f"translation_tension_norm={translation_tension_norm}")
    require(translation_residual_norm <= 1e-12, f"translation_residual_norm={translation_residual_norm}")

    require(stretch_rhs_norm > 0.0, f"stretch_rhs_norm={stretch_rhs_norm}")
    require(stretch_tension_norm > 0.0, f"stretch_tension_norm={stretch_tension_norm}")
    require(stretch_residual_norm <= 1e-8 or stretch_relative_residual <= 1e-10,
            f"stretch_residual_norm={stretch_residual_norm}, stretch_relative_residual={stretch_relative_residual}")

    require(curved_rhs_norm > 0.0, f"curved_rhs_norm={curved_rhs_norm}")
    require(curved_residual_norm <= 1e-8 or curved_relative_residual <= 1e-10,
            f"curved_residual_norm={curved_residual_norm}, curved_relative_residual={curved_relative_residual}")
    require(curved_endpoint_change_norm <= 1e-14,
            f"curved_endpoint_change_norm={curved_endpoint_change_norm}")

    require(tension_matrix_apply_maxdiff <= 1e-10,
            f"tension_matrix_apply_maxdiff={tension_matrix_apply_maxdiff}")

    print("STAGE 2.5 TENSION SOLVER CHECK PASSED")


if __name__ == "__main__":
    main()
