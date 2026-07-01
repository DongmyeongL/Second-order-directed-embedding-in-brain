#!/usr/bin/env python3
"""Conceptual schematic for first-order vs second-order output embedding."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Ellipse
import matplotlib.patheffects as pe


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
MIRROR_DIR = ROOT.parents[0] / "figures"


COLORS = {
    "blue": "#1554B7",
    "blue_dark": "#173E8F",
    "light_blue": "#DDEEFF",
    "orange": "#F59A23",
    "orange_dark": "#E55317",
    "orange_light": "#FFE6C6",
    "green": "#177A3A",
    "green_dark": "#17652A",
    "gray": "#6E747A",
    "light_gray": "#F4F6F8",
    "node_edge": "#202326",
}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
)


def setup_ax(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def add_title(ax, text, color, x0=0.045, y0=0.895, w=0.91, h=0.070):
    box = FancyBboxPatch(
        (x0, y0),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor=color,
        edgecolor="none",
        transform=ax.transAxes,
    )
    box.set_path_effects([pe.SimplePatchShadow(offset=(0, -0.8), alpha=0.10), pe.Normal()])
    ax.add_patch(box)
    ax.text(
        x0 + w / 2,
        y0 + h / 2,
        text,
        color="white",
        fontsize=17,
        fontweight="bold",
        ha="center",
        va="center",
        transform=ax.transAxes,
    )


def add_node(ax, xy, r=0.027, fc="#FFFFFF", ec=COLORS["node_edge"], lw=1.4, label=None, fs=14):
    # Use a marker instead of a data-space Circle so nodes remain circular even
    # when an axis is not square.
    size = (r * 830) ** 2
    circ = ax.scatter(
        [xy[0]],
        [xy[1]],
        s=size,
        marker="o",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=5,
    )
    circ.set_path_effects([pe.SimplePatchShadow(offset=(0.5, -0.5), alpha=0.12), pe.Normal()])
    if label is not None:
        ax.text(
            xy[0],
            xy[1],
            label,
            color="white",
            fontsize=fs,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=6,
        )
    return circ


def add_arrow(
    ax,
    start,
    end,
    color=COLORS["gray"],
    lw=1.2,
    ms=17,
    alpha=1.0,
    zorder=4,
    shrinkA=10,
    shrinkB=14,
):
    arr = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        connectionstyle="arc3,rad=0.0",
        mutation_scale=ms,
        lw=lw,
        color=color,
        alpha=alpha,
        shrinkA=shrinkA,
        shrinkB=shrinkB,
        zorder=zorder,
    )
    ax.add_patch(arr)
    return arr


def add_dotted_box(ax, xy, w, h, color, lw=1.15):
    box = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.050",
        facecolor="none",
        edgecolor=color,
        lw=lw,
        linestyle=(0, (4.5, 4.0)),
        zorder=1,
    )
    ax.add_patch(box)
    return box


def panel_first_order(ax):
    add_title(ax, "First-order Topology", COLORS["blue_dark"])

    # Local context clouds.
    ax.add_patch(Ellipse((0.28, 0.49), 0.18, 0.45, facecolor="#EAF2FF", edgecolor="none", alpha=0.72))
    ax.add_patch(Ellipse((0.72, 0.49), 0.18, 0.45, facecolor="#FFF0E1", edgecolor="none", alpha=0.70))

    s = (0.50, 0.49)
    add_node(ax, s, r=0.050, fc=COLORS["blue"], label="S", fs=18)

    upstream = [(0.08, 0.65), (0.05, 0.53), (0.08, 0.40), (0.14, 0.29)]
    inputs = [(0.21, 0.58), (0.20, 0.47), (0.21, 0.34)]
    for u in upstream:
        add_node(ax, u, r=0.025, fc="#EAF4FF")
    for i, xy in enumerate(inputs):
        add_node(ax, xy, r=0.025, fc="#EAF4FF")
        add_arrow(ax, upstream[i if i < 3 else 2], xy, color="#8B8E90", lw=1.0, ms=14, shrinkB=12)
        add_arrow(ax, xy, s, color="#55585A", lw=1.15, ms=15, shrinkB=24)
    add_arrow(ax, upstream[-1], inputs[-1], color="#8B8E90", lw=1.0, ms=14, shrinkB=12)

    targets = [(0.76, 0.59), (0.76, 0.48), (0.76, 0.34)]
    further = [(0.91, 0.65), (0.94, 0.53), (0.91, 0.40), (0.91, 0.29)]
    for t in targets:
        add_node(ax, t, r=0.026, fc=COLORS["orange"])
        add_arrow(ax, s, t, color="#55585A", lw=1.15, ms=15, shrinkA=24, shrinkB=13)
    for f in further:
        add_node(ax, f, r=0.024, fc="#F4F4F4", ec="#555555")
    add_arrow(ax, targets[0], further[0], color="#8B8E90", lw=1.0, ms=14, shrinkB=12)
    add_arrow(ax, targets[1], further[1], color="#8B8E90", lw=1.0, ms=14, shrinkB=12)
    add_arrow(ax, targets[1], further[2], color="#8B8E90", lw=1.0, ms=14, shrinkB=12)
    add_arrow(ax, targets[2], further[3], color="#8B8E90", lw=1.0, ms=14, shrinkB=12)


def add_upstream_to_source(ax, s):
    upstream = [(0.10, 0.65), (0.07, 0.53), (0.10, 0.40), (0.12, 0.28)]
    for u in upstream:
        add_node(ax, u, r=0.025, fc="#EAF4FF")
        add_arrow(ax, u, s, color="#5A5D60", lw=1.0, ms=14, shrinkB=24)


def panel_high_embedding(ax):
    add_title(ax, "High 2nd-order embedding", COLORS["orange_dark"])
    s = (0.22, 0.49)
    t = (0.42, 0.51)
    add_upstream_to_source(ax, s)
    add_node(ax, s, r=0.050, fc=COLORS["blue"], label="S", fs=18)
    add_node(ax, t, r=0.052, fc=COLORS["orange_dark"], label="", fs=12)
    add_arrow(ax, s, t, color="#4F555B", lw=1.35, ms=17, shrinkA=24, shrinkB=24)

    add_dotted_box(ax, (0.55, 0.19), 0.40, 0.57, COLORS["orange_dark"])
    mid_nodes = [(0.68, 0.66), (0.70, 0.55), (0.68, 0.42), (0.68, 0.29)]
    end_nodes = [(0.89, 0.71), (0.89, 0.61), (0.90, 0.57), (0.89, 0.48), (0.89, 0.36), (0.89, 0.24)]
    for m in mid_nodes:
        add_node(ax, m, r=0.027, fc=COLORS["orange"])
        add_arrow(ax, t, m, color=COLORS["orange_dark"], lw=1.25, ms=16, shrinkA=24, shrinkB=13)
    for e in end_nodes:
        add_node(ax, e, r=0.024, fc="#F4F4F4", ec="#555555")

    pairs = [(0, 0), (0, 1), (1, 2), (1, 3), (2, 4), (3, 5)]
    for mi, ei in pairs:
        add_arrow(ax, mid_nodes[mi], end_nodes[ei], color="#85898C", lw=1.0, ms=14, shrinkB=12)
    ax.text(0.965, 0.64, "$\\cdots$", rotation=90, fontsize=16, ha="center", va="center")
    ax.text(0.965, 0.23, "$\\cdots$", rotation=90, fontsize=16, ha="center", va="center")


def panel_low_embedding(ax):
    add_title(ax, "Low 2nd-order embedding", COLORS["green_dark"])
    s = (0.22, 0.49)
    t = (0.42, 0.51)
    add_upstream_to_source(ax, s)
    add_node(ax, s, r=0.050, fc=COLORS["blue"], label="S", fs=18)
    add_node(ax, t, r=0.050, fc="#F2A13B", label="", fs=12)
    add_arrow(ax, s, t, color="#4F555B", lw=1.35, ms=17, shrinkA=24, shrinkB=23)

    add_dotted_box(ax, (0.55, 0.19), 0.40, 0.57, COLORS["green_dark"])
    input_nodes = [(0.66, 0.66), (0.70, 0.26)]
    mid_nodes = [(0.74, 0.58), (0.76, 0.47), (0.74, 0.35)]
    end_nodes = [(0.88, 0.66), (0.91, 0.55), (0.90, 0.44), (0.88, 0.32)]
    input_arrow_ends = [(0.465, 0.545), (0.465, 0.475)]
    for inp, arrow_end in zip(input_nodes, input_arrow_ends):
        add_node(ax, inp, r=0.024, fc="#F2A13B")
        add_arrow(
            ax,
            inp,
            arrow_end,
            color="#6D7073",
            lw=1.15,
            ms=15,
            zorder=4,
            shrinkA=8,
            shrinkB=0,
        )
    for m in mid_nodes:
        add_node(ax, m, r=0.027, fc=COLORS["orange"])
        add_arrow(ax, t, m, color="#6D7073", lw=1.15, ms=15, shrinkA=23, shrinkB=13)
            
    for e in end_nodes:
        add_node(ax, e, r=0.024, fc="#F4F4F4", ec="#555555")
    for mi, ei in [(0, 0), (0, 1), (1, 2), (2, 3)]:
        add_arrow(ax, mid_nodes[mi], end_nodes[ei], color="#85898C", lw=1.0, ms=14, shrinkB=12)
    ax.text(0.955, 0.50, "$\\cdots$", rotation=90, fontsize=16, ha="center", va="center")
    ax.text(0.955, 0.35, "$\\cdots$", rotation=90, fontsize=16, ha="center", va="center")


def legend_item(ax, x, y, kind, text):
    if kind == "S":
        add_node(ax, (x, y), r=0.030, fc=COLORS["blue"], label="S", fs=15)
    #elif kind == "Tout":
    #    add_node(ax, (x, y), r=0.031, fc=COLORS["orange_dark"], label="T$_{out}$", fs=9)
    #elif kind == "Tbal":
    #    add_node(ax, (x, y), r=0.031, fc="#F39A3B", label="T$_{bal}$", fs=9)
    elif kind == "up":
        add_node(ax, (x, y), r=0.020, fc="#EAF4FF")
    elif kind == "down":
        add_node(ax, (x, y), r=0.020, fc=COLORS["orange"])
    elif kind == "far":
        add_node(ax, (x, y), r=0.020, fc="#F4F4F4", ec="#555555")
    ax.text(x + 0.040, y, text, fontsize=9.5, ha="left", va="center", color="#111111")


def add_legend(fig):
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.13])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = FancyBboxPatch(
        (0.01, 0.06),
        0.98,
        0.86,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        facecolor="white",
        edgecolor="#9A9A9A",
        lw=0.9,
    )
    ax.add_patch(box)
    legend_item(ax, 0.045, 0.53, "S", "Source node\n(upstream)")
    legend_item(ax, 0.195, 0.53, "Tout", "Output-biased\nintermediate node")
    legend_item(ax, 0.370, 0.53, "Tbal", "Balanced\nintermediate node")
    legend_item(ax, 0.535, 0.53, "up", "Upstream nodes\n(inputs)")
    legend_item(ax, 0.670, 0.53, "down", "Downstream nodes\n(outputs)")
    legend_item(ax, 0.820, 0.53, "far", "Further downstream\nnodes")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MIRROR_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(14, 7.8), dpi=300)
    fig.patch.set_facecolor("white")
    for ax in axes:
        setup_ax(ax)

    panel_first_order(axes[0])
    panel_high_embedding(axes[1])
    panel_low_embedding(axes[2])

    for x in [1 / 3, 2 / 3]:
        fig.add_artist(
            plt.Line2D([x, x], [0.055, 0.965], transform=fig.transFigure, color="#D7D7D7", lw=0.9)
        )

    #add_legend(fig)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.03, wspace=0.06)

    for out_dir in [OUT_DIR, MIRROR_DIR]:
        fig.savefig(
            out_dir / "figure_concept_second_order_embedding_reference_style.png",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.03,
        )
        fig.savefig(
            out_dir / "figure_concept_second_order_embedding_reference_style.pdf",
            bbox_inches="tight",
            pad_inches=0.03,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
