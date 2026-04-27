#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import math
import sys


def parse_scalar(lines: list[str], key: str) -> float:
    for line in lines:
      if line.startswith(f"{key} ="):
        return float(line.split("=", 1)[1].strip().split()[0])
    raise KeyError(key)


def parse_center(lines: list[str]) -> tuple[float, float, float]:
    for line in lines:
      if line.startswith("center_of_mass ="):
        rhs = line.split("=", 1)[1].strip().split()
        return float(rhs[0]), float(rhs[1]), float(rhs[2])
    raise KeyError("center_of_mass")


def require(cond: bool, msg: str) -> None:
    if not cond:
        print(f"STAGE 2.0 CHECK FAILED: {msg}")
        sys.exit(1)


def main() -> None:
    path = Path("stage2_outputs/fibre_geometry_check.dat")
    require(path.exists(), f"missing output file: {path}")
    lines = path.read_text().splitlines()

    nl = int(parse_scalar(lines, "nl"))
    length = parse_scalar(lines, "length")
    ds = parse_scalar(lines, "ds")
    max_seg_err = parse_scalar(lines, "max_segment_length_error")
    min_seg = parse_scalar(lines, "min_segment_length")
    max_seg = parse_scalar(lines, "max_segment_length")
    end_to_end = parse_scalar(lines, "end_to_end_distance")
    max_curv = parse_scalar(lines, "max_curvature")
    cmx, cmy, cmz = parse_center(lines)

    tol = 1e-12
    require(nl == 33, f"nl={nl} != 33")
    require(abs(length - 1.0) <= tol, f"length={length}")
    require(abs(ds - 1.0 / 32.0) <= tol, f"ds={ds}")
    require(max_seg_err <= 1e-12, f"max_segment_length_error={max_seg_err}")
    require(abs(min_seg - ds) <= tol, f"min_segment_length={min_seg}, ds={ds}")
    require(abs(max_seg - ds) <= tol, f"max_segment_length={max_seg}, ds={ds}")
    require(abs(end_to_end - length) <= tol, f"end_to_end_distance={end_to_end}, length={length}")
    require(abs(cmx - 0.5) <= tol and abs(cmy) <= tol and abs(cmz) <= tol,
            f"center_of_mass=({cmx}, {cmy}, {cmz})")
    require(max_curv <= 1e-12, f"max_curvature={max_curv}")

    print("STAGE 2.0 GEOMETRY CHECK PASSED")


if __name__ == "__main__":
    main()
