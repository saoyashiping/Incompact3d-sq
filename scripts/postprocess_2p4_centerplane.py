#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_xy_scalar(path: Path):
    data = np.loadtxt(path)
    if data.shape[1] < 3:
        raise ValueError(f"{path} must contain at least columns: x y u")
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    u = data[:, 2].reshape(len(y), len(x))
    return x, y, u


def render_case(report_dir: Path, prefix: str):
    # Expected pre-extracted center-plane files at z = Lz/2:
    # ux_centerplane_tXXXX.dat with columns x y u
    files = sorted(report_dir.glob("ux_centerplane_t*.dat"))
    if len(files) < 2:
        raise FileNotFoundError(f"need at least two center-plane files in {report_dir}")
    pick = [files[0], files[-1]]
    tags = ["t0", "tend"]
    for f, tag in zip(pick, tags):
        x, y, u = load_xy_scalar(f)
        fig, ax = plt.subplots(figsize=(7, 3))
        im = ax.imshow(
            u,
            origin="lower",
            extent=[x.min(), x.max(), y.min(), y.max()],
            aspect="auto",
            cmap="viridis",
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"{prefix} centerplane u ({tag})")
        fig.colorbar(im, ax=ax, label="u")
        out = report_dir.parent / f"velocity_centerplane_{prefix}_{tag}.png"
        fig.tight_layout()
        fig.savefig(out, dpi=180)
        plt.close(fig)


def main():
    run_root = Path("/home/sq/runs24_2p4_round4_codefix")
    render_case(run_root / "parallel_report", "parallel")
    render_case(run_root / "wallnormal_report", "wallnormal")
    render_case(run_root / "deg45_report", "45deg")


if __name__ == "__main__":
    main()
