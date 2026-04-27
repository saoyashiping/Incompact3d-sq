#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys


def require(cond: bool, msg: str) -> None:
    if not cond:
        print(f"STAGE 2.8 CHECK FAILED: {msg}")
        sys.exit(1)


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


def read_time_signal(path: Path) -> tuple[list[float], list[float]]:
    times: list[float] = []
    sig: list[float] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        times.append(float(parts[0]))
        sig.append(float(parts[1]))
    return times, sig


def estimate_frequency(times: list[float], signal: list[float], case_name: str) -> tuple[float, int, int]:
    n = len(times)
    require(n == len(signal) and n >= 10, f"{case_name}: insufficient time-series samples ({n})")

    start = int(0.2 * n)
    if start >= n - 2:
        start = 0

    t = times[start:]
    y = signal[start:]
    mean_y = sum(y) / len(y)
    y = [v - mean_y for v in y]

    crossings: list[float] = []
    for i in range(len(y) - 1):
      y0 = y[i]
      y1 = y[i + 1]
      if y0 <= 0.0 and y1 > 0.0:
            t0 = t[i]
            t1 = t[i + 1]
            denom = (y1 - y0)
            if denom == 0.0:
                tc = t0
            else:
                tc = t0 + (0.0 - y0) * (t1 - t0) / denom
            crossings.append(tc)

    periods: list[float] = []
    for i in range(len(crossings) - 1):
        dt = crossings[i + 1] - crossings[i]
        if dt > 0.0:
            periods.append(dt)

    n_cross = len(crossings)
    n_periods = len(periods)
    require(n_periods >= 3, f"{case_name}: need >=3 periods from upward crossings, got {n_periods}")

    mean_period = sum(periods) / n_periods
    require(mean_period > 0.0, f"{case_name}: non-positive mean period {mean_period}")

    freq = 1.0 / mean_period
    return freq, n_cross, n_periods


def main() -> None:
    out_dir = Path("stage2_outputs")
    summary_path = out_dir / "fibre_frequency_scaling_driver_summary.dat"

    file_map = {
        "base": out_dir / "frequency_base.dat",
        "gamma4": out_dir / "frequency_gamma4.dat",
        "rho4": out_dir / "frequency_rho4.dat",
        "length2": out_dir / "frequency_length2.dat",
    }

    for k, p in file_map.items():
        require(p.exists(), f"missing {k} time series file: {p}")
    require(summary_path.exists(), f"missing summary file: {summary_path}")

    summary = parse_kv(summary_path)
    for case in ("base", "gamma4", "rho4", "length2"):
        require(get_int(summary, f"{case}_nan_detected") == 0,
                f"{case}_nan_detected={get_int(summary, f'{case}_nan_detected')}")
        require(get_int(summary, f"{case}_solver_failure_count") == 0,
                f"{case}_solver_failure_count={get_int(summary, f'{case}_solver_failure_count')}")
        require(get_float(summary, f"{case}_initial_length_error") <= 1e-12,
                f"{case}_initial_length_error={get_float(summary, f'{case}_initial_length_error')}")
        require(get_float(summary, f"{case}_final_length_error") <= 1e-5,
                f"{case}_final_length_error={get_float(summary, f'{case}_final_length_error')}")

    f_base, n_cross_base, _ = estimate_frequency(*read_time_signal(file_map["base"]), case_name="base")
    f_gamma4, n_cross_gamma4, _ = estimate_frequency(*read_time_signal(file_map["gamma4"]), case_name="gamma4")
    f_rho4, n_cross_rho4, _ = estimate_frequency(*read_time_signal(file_map["rho4"]), case_name="rho4")
    f_length2, n_cross_length2, _ = estimate_frequency(*read_time_signal(file_map["length2"]), case_name="length2")

    require(n_cross_base >= 3, f"n_cross_base={n_cross_base}")
    require(n_cross_gamma4 >= 3, f"n_cross_gamma4={n_cross_gamma4}")
    require(n_cross_rho4 >= 3, f"n_cross_rho4={n_cross_rho4}")
    require(n_cross_length2 >= 3, f"n_cross_length2={n_cross_length2}")

    require(f_base > 0.0, f"f_base={f_base}")
    require(f_gamma4 > 0.0, f"f_gamma4={f_gamma4}")
    require(f_rho4 > 0.0, f"f_rho4={f_rho4}")
    require(f_length2 > 0.0, f"f_length2={f_length2}")

    ratio_gamma4_to_base = f_gamma4 / f_base
    ratio_rho4_to_base = f_rho4 / f_base
    ratio_length2_to_base = f_length2 / f_base

    target_ratio_gamma4_to_base = 2.0
    target_ratio_rho4_to_base = 0.5
    target_ratio_length2_to_base = 0.25

    relerr_ratio_gamma4 = abs(ratio_gamma4_to_base - target_ratio_gamma4_to_base) / target_ratio_gamma4_to_base
    relerr_ratio_rho4 = abs(ratio_rho4_to_base - target_ratio_rho4_to_base) / target_ratio_rho4_to_base
    relerr_ratio_length2 = abs(ratio_length2_to_base - target_ratio_length2_to_base) / target_ratio_length2_to_base

    require(relerr_ratio_gamma4 <= 0.20,
            f"relerr_ratio_gamma4={relerr_ratio_gamma4}, ratio_gamma4_to_base={ratio_gamma4_to_base}")
    require(relerr_ratio_rho4 <= 0.20,
            f"relerr_ratio_rho4={relerr_ratio_rho4}, ratio_rho4_to_base={ratio_rho4_to_base}")
    require(relerr_ratio_length2 <= 0.25,
            f"relerr_ratio_length2={relerr_ratio_length2}, ratio_length2_to_base={ratio_length2_to_base}")

    output_path = out_dir / "fibre_frequency_scaling_check.dat"
    with output_path.open("w", encoding="utf-8") as f:
        f.write(f"f_base = {f_base:.16e}\n")
        f.write(f"f_gamma4 = {f_gamma4:.16e}\n")
        f.write(f"f_rho4 = {f_rho4:.16e}\n")
        f.write(f"f_length2 = {f_length2:.16e}\n")
        f.write(f"n_cross_base = {n_cross_base}\n")
        f.write(f"n_cross_gamma4 = {n_cross_gamma4}\n")
        f.write(f"n_cross_rho4 = {n_cross_rho4}\n")
        f.write(f"n_cross_length2 = {n_cross_length2}\n")
        f.write(f"ratio_gamma4_to_base = {ratio_gamma4_to_base:.16e}\n")
        f.write(f"ratio_rho4_to_base = {ratio_rho4_to_base:.16e}\n")
        f.write(f"ratio_length2_to_base = {ratio_length2_to_base:.16e}\n")
        f.write(f"target_ratio_gamma4_to_base = {target_ratio_gamma4_to_base:.16e}\n")
        f.write(f"target_ratio_rho4_to_base = {target_ratio_rho4_to_base:.16e}\n")
        f.write(f"target_ratio_length2_to_base = {target_ratio_length2_to_base:.16e}\n")
        f.write(f"relerr_ratio_gamma4 = {relerr_ratio_gamma4:.16e}\n")
        f.write(f"relerr_ratio_rho4 = {relerr_ratio_rho4:.16e}\n")
        f.write(f"relerr_ratio_length2 = {relerr_ratio_length2:.16e}\n")

    print("STAGE 2.8 FREQUENCY SCALING CHECK PASSED")


if __name__ == "__main__":
    main()
