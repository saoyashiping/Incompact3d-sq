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
        print(f"STAGE 2.1 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_freefree_boundary_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    nl = get_int(data, "nl")
    length = get_float(data, "length")
    ds = get_float(data, "ds")
    rho_tilde = get_float(data, "rho_tilde")
    gamma = get_float(data, "gamma")
    center_of_mass_x = get_float(data, "center_of_mass_x")
    max_boundary = get_float(data, "max_freefree_boundary_residual")
    left_change = get_float(data, "left_endpoint_change_norm")
    right_change = get_float(data, "right_endpoint_change_norm")
    tension_size = get_int(data, "tension_size")

    tol = 1e-12
    require(nl == 33, f"nl={nl} != 33")
    require(abs(length - 1.0) <= tol, f"length={length}")
    require(abs(ds - 1.0 / 32.0) <= tol, f"ds={ds}")
    require(abs(rho_tilde - 1.0) <= tol, f"rho_tilde={rho_tilde}")
    require(abs(gamma - 1.0) <= tol, f"gamma={gamma}")
    require(abs(center_of_mass_x - 0.5) <= 5e-4, f"center_of_mass_x={center_of_mass_x}")
    require(max_boundary <= 1e-10, f"max_freefree_boundary_residual={max_boundary}")
    require(left_change <= 1e-14, f"left_endpoint_change_norm={left_change}")
    require(right_change <= 1e-14, f"right_endpoint_change_norm={right_change}")
    require(tension_size == nl - 1, f"tension_size={tension_size}, nl-1={nl-1}")

    print("STAGE 2.1 FREE-FREE BOUNDARY CHECK PASSED")


if __name__ == "__main__":
    main()
