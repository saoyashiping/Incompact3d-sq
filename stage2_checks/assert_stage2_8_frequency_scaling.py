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
        print(f"STAGE 2.8 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_frequency_scaling_check.dat")
    require(path.exists(), f"missing output file: {path}")
    data = parse_kv(path)

    freq_gamma1_dt = get_float(data, "freq_gamma1_dt")
    freq_gamma1_dt_half = get_float(data, "freq_gamma1_dt_half")
    freq_gamma4_dt = get_float(data, "freq_gamma4_dt")
    freq_dt_relative_diff = get_float(data, "freq_dt_relative_diff")
    freq_gamma_scaling_ratio = get_float(data, "freq_gamma_scaling_ratio")
    expected_gamma_scaling_ratio = get_float(data, "expected_gamma_scaling_ratio")
    solver_failure_count_total = get_int(data, "solver_failure_count_total")
    nan_detected = get_int(data, "nan_detected")

    require(freq_gamma1_dt > 0.0, f"freq_gamma1_dt={freq_gamma1_dt}")
    require(freq_gamma1_dt_half > 0.0, f"freq_gamma1_dt_half={freq_gamma1_dt_half}")
    require(freq_gamma4_dt > 0.0, f"freq_gamma4_dt={freq_gamma4_dt}")
    require(freq_dt_relative_diff <= 5e-2, f"freq_dt_relative_diff={freq_dt_relative_diff}")
    require(abs(freq_gamma_scaling_ratio - expected_gamma_scaling_ratio) <= 1.5e-1,
            f"freq_gamma_scaling_ratio={freq_gamma_scaling_ratio}, expected={expected_gamma_scaling_ratio}")
    require(solver_failure_count_total == 0, f"solver_failure_count_total={solver_failure_count_total}")
    require(nan_detected == 0, f"nan_detected={nan_detected}")

    print("STAGE 2.8 FREQUENCY SCALING CHECK PASSED")


if __name__ == "__main__":
    main()
