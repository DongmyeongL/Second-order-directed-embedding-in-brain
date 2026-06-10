"""Compact comparative figure combining the core messages of old Figs. 5--6.

The figure is regenerated from source tables and plotting functions, not by
cropping or pasting existing PNG panels.
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
import figure14_15_celegans_drosophila_AFHI_combined as fig5
import figure7_by_species_with_figure5_abc as fig6


OUT_FIG = PROJECT_ROOT / "figures" / "figure5_compact_comparative.png"
OUT_EFFECTS = PROJECT_ROOT / "data" / "final_summary_tables" / "figure5_compact_comparative_effects.csv"
OUT_GROUP = PROJECT_ROOT / "data" / "final_summary_tables" / "figure5_compact_comparative_group_signature.csv"
OUT_GROUP_PERM = PROJECT_ROOT / "data" / "final_summary_tables" / "figure5_compact_comparative_group_permutation_tests.csv"


def add_label(ax: plt.Axes, label: str, x: float = -0.18, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=fs.PANEL_LABEL_FS_2COL,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def main() -> None:
    fs.apply_main_figure_style()
    fig5.apply_local_group_colors()

    tables = fig6.use_bilateral_zebrafish_ob_fcv(fig6.collapse_metric_tables(fig6.metric_tables()))
    tables, effects, group_signature, group_perm = fig6.load_or_compute_statistics(
        tables,
        recompute=False,
        effects_path=OUT_EFFECTS,
        group_path=OUT_GROUP,
        group_perm_path=OUT_GROUP_PERM,
    )
    effects, group_signature, group_perm = fig6.relabel_collapsed_group_tables(
        effects,
        group_signature,
        group_perm,
    )
    effects.to_csv(OUT_EFFECTS, index=False)
    group_signature.to_csv(OUT_GROUP, index=False)
    group_perm.to_csv(OUT_GROUP_PERM, index=False)

    fig = plt.figure(figsize=(8.6, 8.25))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.05, 2.45],
        left=0.065,
        right=0.985,
        bottom=0.115,
        top=0.930,
        hspace=0.46,
    )
    top_gs = outer[0].subgridspec(1, 5, wspace=0.50)
    sig_gs = outer[1].subgridspec(3, 3, wspace=0.36, hspace=0.46)

    ax_ce_fcv = fig.add_subplot(top_gs[0, 0])
    ax_ce_post = fig.add_subplot(top_gs[0, 1])
    ax_fly_fcv = fig.add_subplot(top_gs[0, 2])
    ax_fly_post = fig.add_subplot(top_gs[0, 3])
    ax_fly_pre = fig.add_subplot(top_gs[0, 4])

    fig5.draw_fcv_anatomy_panel(ax_ce_fcv, "C. elegans", show_ylabel=True)
    fig5.draw_sc_anatomy_panel(
        ax_ce_post,
        "C. elegans",
        "PostDCA",
        fig5.DCA_POST_LABEL,
        show_ylabel=True,
    )
    fig5.draw_fcv_anatomy_panel(ax_fly_fcv, "Drosophila", show_ylabel=True)
    fig5.draw_sc_anatomy_panel(
        ax_fly_post,
        "Drosophila",
        "PostDCA",
        fig5.DCA_POST_LABEL,
        show_ylabel=True,
    )
    fig5.draw_sc_anatomy_panel(
        ax_fly_pre,
        "Drosophila",
        "PreDCA",
        fig5.DCA_PRE_LABEL,
        show_ylabel=True,
    )

    for ax, title, color in [
        (ax_ce_fcv, "FCV", fig6.SPECIES_COLORS["C. elegans"]),
        (ax_ce_post, r"$\mathrm{DCA}_{\mathrm{post}}$", fig6.SPECIES_COLORS["C. elegans"]),
        (ax_fly_fcv, "FCV", fig6.SPECIES_COLORS["Drosophila"]),
        (ax_fly_post, r"$\mathrm{DCA}_{\mathrm{post}}$", fig6.SPECIES_COLORS["Drosophila"]),
        (ax_fly_pre, r"$\mathrm{DCA}_{\mathrm{pre}}$", fig6.SPECIES_COLORS["Drosophila"]),
    ]:
        ax.set_title(title, fontsize=fs.AXIS_LABEL_FS_2COL - 0.5, fontweight="bold", color=color, pad=3)
        ax.tick_params(axis="both", labelsize=fs.TICK_FS_2COL - 1)

    top_axes = [ax_ce_fcv, ax_ce_post, ax_fly_fcv, ax_fly_post, ax_fly_pre]
    for ax, label in zip(top_axes, list("ABCDE"), strict=True):
        add_label(ax, label, x=-0.18, y=1.12)

    signature_axes_grid = [
        [fig.add_subplot(sig_gs[row, col]) for col in range(3)]
        for row in range(3)
    ]
    signature_axes = [ax for row in signature_axes_grid for ax in row]
    fig6.draw_metric_species_grid(signature_axes_grid, tables, group_perm)
    fig6.style_signature_axes(signature_axes)
    for ax, label in zip([row[0] for row in signature_axes_grid], ["F", "G", "H"], strict=True):
        add_label(ax, label)

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight", transparent=False)
    plt.close(fig)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_EFFECTS}")
    print(f"wrote {OUT_GROUP}")
    print(f"wrote {OUT_GROUP_PERM}")


if __name__ == "__main__":
    main()
