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
        print(f"STAGE 2.2 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_bending_operator_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    straight_max_abs_d4 = get_float(data, "straight_max_abs_d4")
    straight_max_abs_fb = get_float(data, "straight_max_abs_fb")
    cubic_inner_max_abs_d4 = get_float(data, "cubic_inner_max_abs_d4")
    sine_inner_relative_l2_error_d4 = get_float(data, "sine_inner_relative_l2_error_d4")
    d4_left_endpoint_change_norm = get_float(data, "d4_left_endpoint_change_norm")
    d4_right_endpoint_change_norm = get_float(data, "d4_right_endpoint_change_norm")
    sine_l2_error_nl33 = get_float(data, "sine_l2_error_nl33")
    sine_l2_error_nl65 = get_float(data, "sine_l2_error_nl65")
    sine_l2_error_nl129 = get_float(data, "sine_l2_error_nl129")

    require(straight_max_abs_d4 <= 1e-10, f"straight_max_abs_d4={straight_max_abs_d4}")
    require(straight_max_abs_fb <= 1e-10, f"straight_max_abs_fb={straight_max_abs_fb}")
    require(cubic_inner_max_abs_d4 <= 1e-8, f"cubic_inner_max_abs_d4={cubic_inner_max_abs_d4}")
    require(sine_inner_relative_l2_error_d4 <= 0.15,
            f"sine_inner_relative_l2_error_d4={sine_inner_relative_l2_error_d4}")
    require(d4_left_endpoint_change_norm <= 1e-14,
            f"d4_left_endpoint_change_norm={d4_left_endpoint_change_norm}")
    require(d4_right_endpoint_change_norm <= 1e-14,
            f"d4_right_endpoint_change_norm={d4_right_endpoint_change_norm}")
    require(sine_l2_error_nl65 < sine_l2_error_nl33,
            f"sine_l2_error_nl65={sine_l2_error_nl65}, sine_l2_error_nl33={sine_l2_error_nl33}")
    require(sine_l2_error_nl129 < sine_l2_error_nl65,
            f"sine_l2_error_nl129={sine_l2_error_nl129}, sine_l2_error_nl65={sine_l2_error_nl65}")

    print("STAGE 2.2 BENDING OPERATOR CHECK PASSED")


if __name__ == "__main__":
    main()
