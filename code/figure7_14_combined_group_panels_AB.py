"""Panel set A-B from the combined cross-species group figure."""

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


OUT_FIG = PROJECT_ROOT / "figures" / "figure7_14_combined_group_panels_AB.png"
OUT_OUTPUT = PROJECT_ROOT / "outputs" / "figure7_14_combined_group_panels_AB.png"


def main() -> None:
    fs.apply_main_figure_style()
    fig = plt.figure(figsize=(8.0, 5))
    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.0, 1.0],
        left=0.100,
        right=0.985,
        bottom=0.145,
        top=0.955,
        wspace=0.50,
        hspace=0.62,
    )
    axes = combo.draw_anatomy_rows(fig, gs)
    combo.shift_axes(axes[:3], dy=0.010)
    combo.shift_axes(axes[3:6], dy=-0.012)
    #for ax in axes:
    #    ax.yaxis.set_label_coords(-0.3, 0.5)
    # Per-panel y-label offsets. Axes order: A=0, B=1, C=2, D=3, E=4, F=5.
    for idx in [0, 3]:
        axes[idx].yaxis.set_label_coords(-0.25, 0.5)
    
    for idx in [1, 4]:
        axes[idx].yaxis.set_label_coords(-0.28, 0.5)
        
    for idx in [2, 5]:
        axes[idx].yaxis.set_label_coords(-0.25, 0.5)

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    OUT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight", transparent=False)
    fig.savefig(OUT_OUTPUT, dpi=600, bbox_inches="tight", transparent=False)
    plt.close(fig)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_OUTPUT}")


if __name__ == "__main__":
    main()
