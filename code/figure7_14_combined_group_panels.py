"""Combined cross-species group panels from Figures 7 and 14/15.

This figure is redrawn from the underlying tables and plotting utilities rather
than assembled by cropping existing PNG files.
"""

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
import figure14_15_celegans_drosophila_AFHI_combined as fig14
import figure7_by_species_with_figure5_abc as fig7

OUT_FIG = PROJECT_ROOT / "figures" / "figure7_14_combined_group_panels.png"
OUT_OUTPUT = PROJECT_ROOT / "outputs" / "figure7_14_combined_group_panels.png"

# Direct row-position adjustments in figure coordinates. Positive values move
# rows upward; negative values move rows downward.
CE_ROW_DY = 0.030
FLY_ROW_DY = -0.005
FUNCTIONAL_ROW_DY = -0.020

SHORT_ANATOMY_LABELS = {
    "olf./chemo": "olf.",
    "sensory neurons": "sensory",
    "interneurons": "inter",
    "motor / command": "motor",
    "state-modulatory": "state",
    "olfactory neuropil": "olf.",
    "mushroom body": "MB",
    "optic neuropil": "optic",
    "central complex": "CX",
    "protocerebrum": "proto.",
    "premotor interface": "premotor",
}

SHORT_FUNCTIONAL_LABELS = {
    "olf./chemo": "olf.",
    "assoc./state": "assoc.",
    "other sensory": "sensory",
    "sensorimotor": "sensori.",
    "integrative": "integ.",
}


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.18, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fs.PANEL_LABEL_FS_2COL,
        fontweight="bold",
    )


def shift_axes(axes: list[plt.Axes], dx: float = 0.0, dy: float = 0.0) -> None:
    for ax in axes:
        pos = ax.get_position()
        ax.set_position([pos.x0 + dx, pos.y0 + dy, pos.width, pos.height])


def shorten_xticklabels(ax: plt.Axes, mapping: dict[str, str], rotation: float | None = None) -> None:
    labels = [tick.get_text() for tick in ax.get_xticklabels()]
    ax.set_xticklabels(
        [mapping.get(label, label) for label in labels],
        rotation=rotation if rotation is not None else ax.get_xticklabels()[0].get_rotation() if labels else 0,
        ha="right",
    )


def draw_anatomy_rows(fig: plt.Figure, gs) -> list[plt.Axes]:
    """Draw Fig. 14/15 group panels C-E and H-J."""
    fig14.apply_local_group_colors()
    row_specs = [
        ("C. elegans", "A", 0),
        ("Drosophila", "D", 1),
    ]
    panel_labels = [["A", "B", "C"], ["D", "E", "F"]]
    measures = [
        ("fcv", "FCV"),
        ("PostDCA", fig14.DCA_POST_LABEL),
        ("PreDCA", fig14.DCA_PRE_LABEL),
    ]
    axes: list[plt.Axes] = []
    for species, row_label, row in row_specs:
        for col, (measure, label) in enumerate(measures):
            ax = fig.add_subplot(gs[row, col])
            axes.append(ax)
            if measure == "fcv":
                fig14.draw_fcv_anatomy_panel(ax, species, show_ylabel=True)
            else:
                fig14.draw_sc_anatomy_panel(ax, species, measure, label, show_ylabel=True)
            shorten_xticklabels(ax, SHORT_ANATOMY_LABELS, rotation=38)
            if col == -1:
                ax.text(
                    -0.42,
                    0.50,
                    species,
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=fs.AXIS_LABEL_FS_2COL,
                    fontweight="bold",
                    style="italic" if species.startswith("C.") else "normal",
                )
            label_row = 0 if species == "C. elegans" else 1
            add_panel_label(ax, panel_labels[label_row][col],x=-0.45)
    return axes


def draw_functional_group_rows(fig: plt.Figure, gs, start_row: int = 2) -> list[plt.Axes]:
    """Draw Fig. 7 functional-group synthesis panels D-F."""
    tables = fig7.use_bilateral_zebrafish_ob_fcv(fig7.collapse_metric_tables(fig7.metric_tables()))
    tables, effects, group_signature, group_perm = fig7.load_or_compute_statistics(
        tables,
        recompute=False,
        effects_path=fig7.OUT_EFFECTS,
        group_path=fig7.OUT_GROUP,
        group_perm_path=fig7.OUT_GROUP_PERM,
    )
    effects, group_signature, group_perm = fig7.relabel_collapsed_group_tables(
        effects,
        group_signature,
        group_perm,
    )

    axes_grid = [[fig.add_subplot(gs[row, col]) for col in range(3)] for row in range(start_row, start_row + 3)]
    fig7.draw_metric_species_grid(axes_grid, tables, group_perm)
    fig7.style_signature_axes([ax for row in axes_grid for ax in row])
    for ax in [ax for row in axes_grid for ax in row]:
        shorten_xticklabels(ax, SHORT_FUNCTIONAL_LABELS, rotation=32)
    for ax, label in zip([row[0] for row in axes_grid], ["A", "B", "C"], strict=False):
        add_panel_label(ax, label,x=-0.35)
    return [ax for row in axes_grid for ax in row]


def main() -> None:
    fs.apply_main_figure_style()
    fig = plt.figure(figsize=(8, 12))
    gs = fig.add_gridspec(
        5,
        3,
        height_ratios=[0.78, 0.78, 0.78, 0.78, 0.78],
        left=0.100,
        right=0.985,
        bottom=0.085,
        top=0.975,
        wspace=0.50,
        hspace=0.46,
    )

    anatomy_axes = draw_anatomy_rows(fig, gs)
    functional_axes = draw_functional_group_rows(fig, gs)

    # Fine-tune row spacing by moving axes directly. This is easier to control
    # than changing GridSpec hspace when only one gap needs adjustment.
    ce_axes = anatomy_axes[:3]
    fly_axes = anatomy_axes[3:6]
    shift_axes(ce_axes, dy=CE_ROW_DY)
    shift_axes(fly_axes, dy=FLY_ROW_DY+0.003)
    
    shift_axes(functional_axes[0:3], dy=-0.015)
    shift_axes(functional_axes[3:6], dy=0.01)
    shift_axes(functional_axes[6:9], dy=0.04)
    
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    OUT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight", transparent=False)
    fig.savefig(OUT_OUTPUT, dpi=600, bbox_inches="tight", transparent=False)
    plt.close(fig)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_OUTPUT}")


if __name__ == "__main__":
    main()
