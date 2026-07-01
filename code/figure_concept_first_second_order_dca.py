import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PROJECT_ROOT.parent
OUT_FINAL = PROJECT_ROOT / "figures"
OUT_MANUSCRIPT = WORKSPACE_ROOT / "figures"

TEAL = "#2B7585"
GREEN = "#20DB19"
MAGENTA = "#D927C9"
GRAY = "#D6D6D6"
HIGH_BG = "#F4F62A"
LOW_BG = "#FF2329"
RED_TOP = "#FF604C"
RED_BOTTOM = "#B91612"
BLUE_TOP = "#50BDF2"
BLUE_BOTTOM = "#0E80B7"


def setup_style():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )


def rounded_panel(ax, x, y, w, h, color):
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.025,rounding_size=0.14",
            facecolor=color,
            edgecolor="none",
            zorder=0,
        )
    )


def arrow(ax, start, end, color, width=7.0, head=24, alpha=0.96, z=3):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            color=color,
            linewidth=width,
            mutation_scale=head,
            shrinkA=3,
            shrinkB=9,
            capstyle="butt",
            joinstyle="miter",
            alpha=alpha,
            zorder=z,
        )
    )


def node(ax, xy, r=0.055):
    ax.add_patch(Circle(xy, r, facecolor=GRAY, edgecolor=GRAY, linewidth=1.0, zorder=5))


def gradient_circle(ax, xy, r, top, bottom, edge, label, fs=30):
    n = 256
    grad = np.linspace(0, 1, n)[:, None]
    top_rgb = np.array(mpl.colors.to_rgb(top))
    bot_rgb = np.array(mpl.colors.to_rgb(bottom))
    img = (top_rgb * (1 - grad) + bot_rgb * grad)[:, None, :]
    im = ax.imshow(img, extent=[xy[0] - r, xy[0] + r, xy[1] - r, xy[1] + r], origin="upper", zorder=7)
    clip = Circle(xy, r, facecolor="none", edgecolor=edge, linewidth=1.0, zorder=8)
    ax.add_patch(clip)
    im.set_clip_path(clip)
    ax.text(xy[0], xy[1] - 0.003, label, ha="center", va="center", color="black", fontsize=fs, zorder=9)


def draw_inputs(ax, a):
    starts = [(a[0] - 0.28, a[1] + 0.13), (a[0] - 0.31, a[1]), (a[0] - 0.28, a[1] - 0.13)]
    ends = [(a[0] - 0.060, a[1] + 0.055), (a[0] - 0.060, a[1]), (a[0] - 0.060, a[1] - 0.055)]
    for start, end in zip(starts, ends):
        arrow(ax, start, end, TEAL, width=6.0, head=24, alpha=0.95, z=2)


def draw_panel(ax, x0, y0, w, h, title, bg, high=True):
    rounded_panel(ax, x0, y0, w, h, bg)
    ax.text(x0 + 0.060, y0 + h - 0.050, title, ha="left", va="center", fontsize=22, color="black", zorder=10)

    a = (x0 + 0.31, y0 + 0.36)
    b = (x0 + 0.63, y0 + 0.36)
    ar = 0.070
    br = 0.070

    draw_inputs(ax, a)
    arrow(ax, (a[0] + ar * 1.00, a[1]), (b[0] - br * 1.00, b[1]), TEAL, width=6.5, head=26, alpha=0.96, z=2)

    if high:
        ax.text(x0 + 0.200, y0 + 0.510, "First-order Topology", ha="left", va="center", fontsize=14, color="black", zorder=10)
        ax.text(x0 + 0.505, y0 + 0.070, "Second-order embedding", ha="left", va="center", fontsize=14, color="black", zorder=10)
        out_nodes = [
            ((b[0] + 0.060, b[1] + 0.205), (b[0] + 0.012, b[1] + br * 0.95)),
            ((b[0] + 0.200, b[1] + 0.115), (b[0] + br * 0.78, b[1] + 0.040)),
            ((b[0] + 0.215, b[1]), (b[0] + br * 0.96, b[1])),
            ((b[0] + 0.160, b[1] - 0.145), (b[0] + br * 0.62, b[1] - 0.048)),
            ((b[0], b[1] - 0.200), (b[0], b[1] - br * 0.95)),
        ]
        in_nodes = [((b[0] + 0.200, b[1] + 0.090), (b[0] + br * 0.78, b[1] + 0.038))]
    else:
        out_nodes = [
            ((b[0] + 0.060, b[1] + 0.195), (b[0] + 0.012, b[1] + br * 0.95)),
            ((b[0] + 0.215, b[1]), (b[0] + br * 0.96, b[1])),
            ((b[0], b[1] - 0.195), (b[0], b[1] - br * 0.95)),
        ]
        in_nodes = [
            ((b[0] + 0.200, b[1] + 0.120), (b[0] + br * 0.78, b[1] + 0.045)),
            ((b[0] + 0.160, b[1] - 0.145), (b[0] + br * 0.62, b[1] - 0.050)),
        ]

    for pt, start in out_nodes:
        node(ax, pt, r=0.052)
        arrow(ax, start, pt, GREEN, width=5.5, head=23, alpha=0.95, z=3)

    for pt, end in in_nodes:
        node(ax, pt, r=0.052)
        arrow(ax, pt, end, MAGENTA, width=5.5, head=23, alpha=0.92, z=4)

    gradient_circle(ax, a, ar, RED_TOP, RED_BOTTOM, "#D94A3E", "A", fs=28)
    gradient_circle(ax, b, br, BLUE_TOP, BLUE_BOTTOM, "#218FC4", "B", fs=28)


def main():
    setup_style()
    OUT_FINAL.mkdir(parents=True, exist_ok=True)
    OUT_MANUSCRIPT.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(10.8, 4.0), facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 2.10)
    ax.set_ylim(0, 0.78)
    ax.set_aspect("equal")
    ax.axis("off")

    draw_panel(ax, 0.030, 0.035, 0.955, 0.705, "High second-order embedding", HIGH_BG, high=True)
    draw_panel(ax, 1.115, 0.035, 0.955, 0.705, "Low second-order embedding", LOW_BG, high=False)

    for out_dir in (OUT_FINAL, OUT_MANUSCRIPT):
        fig.savefig(out_dir / "figure_concept_first_second_order_dca.png", dpi=600, facecolor="black")
        fig.savefig(out_dir / "figure_concept_first_second_order_dca.pdf", facecolor="black")
    plt.close(fig)


if __name__ == "__main__":
    main()
