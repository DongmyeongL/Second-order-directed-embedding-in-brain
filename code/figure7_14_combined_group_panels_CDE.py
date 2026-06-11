"""Panel set C-E from the combined cross-species group figure."""

from __future__ import annotations

import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

import figure_style as fs
import figure7_14_combined_group_panels as combo


OUT_FIG = PROJECT_ROOT / "figures" / "figure7_14_combined_group_panels_CDE.png"
OUT_OUTPUT = PROJECT_ROOT / "outputs" / "figure7_14_combined_group_panels_CDE.png"


def main() -> None:
    fs.apply_main_figure_style()
    fig = plt.figure(figsize=(8.0, 6.3))
    gs = fig.add_gridspec(
        3,
        3,
        height_ratios=[1.0, 1.0, 1.0],
        left=0.100,
        right=0.985,
        bottom=0.135,
        top=0.940,
        wspace=0.50,
        hspace=0.58,
    )
    axes = combo.draw_functional_group_rows(fig, gs, start_row=0)
    combo.shift_axes(axes[:3], dy=-0.00)
    combo.shift_axes(axes[3:6], dy=0.08)
    combo.shift_axes(axes[6:9], dy=0.15)
    for ax in axes:
        ax.yaxis.set_label_coords(-0.20, 0.5)

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    OUT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight", transparent=False)
    fig.savefig(OUT_OUTPUT, dpi=600, bbox_inches="tight", transparent=False)
    plt.close(fig)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_OUTPUT}")


if __name__ == "__main__":
    main()
