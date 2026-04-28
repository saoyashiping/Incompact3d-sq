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
        print(f"STAGE 3.0 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage3_outputs/fibre_ibm_skeleton_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    require(get_int(data, "grid_nx") == 16, f"grid_nx={get_int(data, 'grid_nx')}")
    require(get_int(data, "grid_ny") == 12, f"grid_ny={get_int(data, 'grid_ny')}")
    require(get_int(data, "grid_nz") == 10, f"grid_nz={get_int(data, 'grid_nz')}")

    dx = get_float(data, "grid_dx")
    dy = get_float(data, "grid_dy")
    dz = get_float(data, "grid_dz")
    cell_volume = get_float(data, "grid_cell_volume")

    require(abs(dx - (2.0 / 16.0)) <= 1e-12, f"grid_dx={dx}")
    require(abs(dy - (2.0 / 12.0)) <= 1e-12, f"grid_dy={dy}")
    require(abs(dz - (1.0 / 10.0)) <= 1e-12, f"grid_dz={dz}")
    require(abs(cell_volume - dx * dy * dz) <= 1e-12,
            f"grid_cell_volume={cell_volume}, dx*dy*dz={dx*dy*dz}")

    require(get_int(data, "grid_periodic_x") == 1, f"grid_periodic_x={get_int(data, 'grid_periodic_x')}")
    require(get_int(data, "grid_periodic_y") == 0, f"grid_periodic_y={get_int(data, 'grid_periodic_y')}")
    require(get_int(data, "grid_periodic_z") == 1, f"grid_periodic_z={get_int(data, 'grid_periodic_z')}")

    lag_nl = get_int(data, "lag_nl")
    lag_weight_sum = get_float(data, "lag_weight_sum")
    expected_fibre_length = get_float(data, "expected_fibre_length")

    require(lag_nl == 33, f"lag_nl={lag_nl}")
    require(abs(lag_weight_sum - expected_fibre_length) <= 1e-12,
            f"lag_weight_sum={lag_weight_sum}, expected_fibre_length={expected_fibre_length}")
    require(get_int(data, "lag_inside_count") == lag_nl,
            f"lag_inside_count={get_int(data, 'lag_inside_count')}, lag_nl={lag_nl}")
    require(get_int(data, "lag_outside_count") == 0,
            f"lag_outside_count={get_int(data, 'lag_outside_count')}")

    require(get_int(data, "stencil_nl") == lag_nl,
            f"stencil_nl={get_int(data, 'stencil_nl')}, lag_nl={lag_nl}")
    require(get_int(data, "stencil_max_points_per_lag") == 0,
            f"stencil_max_points_per_lag={get_int(data, 'stencil_max_points_per_lag')}")
    require(get_int(data, "stencil_total_count") == 0,
            f"stencil_total_count={get_int(data, 'stencil_total_count')}")

    require(get_int(data, "ibm_skeleton_status") == 1,
            f"ibm_skeleton_status={get_int(data, 'ibm_skeleton_status')}")

    print("STAGE 3.0 IBM SKELETON CHECK PASSED")


if __name__ == "__main__":
    main()
