#!/usr/bin/env python3
"""Extract channel mean profile from statistics text output.

Usage:
  python scripts/postprocess_2p7_mean_profile.py --input statistics/umean.dat0000400 --output umean_profile.dat
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to statistics file (e.g. umean.dat*)")
    ap.add_argument("--output", required=True, help="Output profile file")
    args = ap.parse_args()

    arr = np.loadtxt(args.input)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    # Expecting (x,y,z,value) text export. Fall back to first/last columns when needed.
    if arr.shape[1] >= 4:
        y = arr[:, 1]
        u = arr[:, -1]
    elif arr.shape[1] >= 2:
        y = arr[:, 0]
        u = arr[:, -1]
    else:
        raise ValueError("Input statistics file must have at least 2 columns.")

    y_unique = np.unique(y)
    prof = np.zeros((y_unique.size, 2), dtype=float)
    prof[:, 0] = y_unique
    for i, yi in enumerate(y_unique):
        prof[i, 1] = u[y == yi].mean()

    header = "y mean_value"
    np.savetxt(args.output, prof, header=header, comments="")


if __name__ == "__main__":
    main()
