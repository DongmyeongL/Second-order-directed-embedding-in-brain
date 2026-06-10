"""
Combined C. elegans and Drosophila figure.

Rows:
  C. elegans: network | DCA-FCV scatter | FCV box | Post-DCA box | Pre-DCA box
  Drosophila: network | DCA-FCV scatter | FCV box | Post-DCA box | Pre-DCA box
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import itertools

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

import figure14_celegans_full_combined_FINAL as celegans
import figure15_drosophila_full_combined_FINAL as drosophila
import figure1_cross_species_fc_dynamics_measures as fc_dynamics
import figure2_structural_connectivity_predictors_of_fcv as sc_predictors
import figure4_post_dca_distribution_summary as post_dca_scatter
import figure_style as fs


OUT_PNG = PROJECT_ROOT / "figures" / "figure14_15_celegans_drosophila_AFHI_combined.png"
SC_VALUES_CSV = PROJECT_ROOT / "data" / "final_summary_tables" / "sc_four_measures_vs_fcv_all_species_values.csv"
INHERITED_PANEL_LABELS = {"A", "B", "C", "D", "E", "F", "H", "I"}
PANEL_LABELS_BY_ROW = [["A", "B", "C", "D", "E"], ["F", "G", "H", "I", "J"]]
BOTTOM_ROW_SHIFT_Y = -0.020
SPECIES_TITLE_OFFSET_Y = 0.090
PANEL_LABEL_Y_PAD = 0.016
DCA_SCATTER_SHIFT_X = -0.018
FCV_PANEL_SHIFT_X = -0.010
DCA_POST_LABEL = celegans.DCA_POST_LABEL
DCA_PRE_LABEL = celegans.DCA_PRE_LABEL
AXIS_FS = fs.AXIS_LABEL_FS_2COL
TICK_FS = fs.TICK_FS_2COL
PANEL_FS = fs.PANEL_LABEL_FS_2COL
STAT_FS = fs.STAT_FS_2COL
STAR_FS = fs.STAR_FS_2COL
GROUP_LABEL_FS = fs.TICK_FS_2COL - 2
NETWORK_LABEL_FS = fs.TICK_FS_2COL - 2
SIG_BAR_LW = 0.65
SIG_BAR_COLOR = "#333333"
SIG_BAR_Y_START = 0.035
SIG_BAR_STEP = 0.062
SIG_BAR_HEIGHT_FRAC = 0.22

DIVISION_COLORS = fs.ZEBRAFISH_DIVISION_COLORS
ADDED_GROUP_COLORS = {
    "associative": "#B65F73",
    "modulatory": "#8E6BBE",
    "other": fs.MAIN_COLORS["neutral_light"],
}
ANATOMY_COLORS = {
    "olf./chemo": DIVISION_COLORS["Tel"],
    "sensory neurons": DIVISION_COLORS["Di"],
    "interneurons": DIVISION_COLORS["Mes"],
    "motor / command": DIVISION_COLORS["Hind"],
    "state-modulatory": ADDED_GROUP_COLORS["modulatory"],
    "olfactory neuropil": DIVISION_COLORS["Tel"],
    "optic neuropil": DIVISION_COLORS["Hind"],
    "mushroom body": ADDED_GROUP_COLORS["associative"],
    "central complex": DIVISION_COLORS["Di"],
    "protocerebrum": DIVISION_COLORS["Mes"],
    "premotor interface": fs.MAIN_COLORS["neutral"],
    "other": ADDED_GROUP_COLORS["other"],
}
CLASS_COLORS = {
    0: DIVISION_COLORS["Tel"],               # olfactory / chemosensory
    1: DIVISION_COLORS["Hind"],              # visual / optic
    2: DIVISION_COLORS["Di"],                # other sensory
    3: DIVISION_COLORS["Hind"],              # sensorimotor output
    4: DIVISION_COLORS["Mes"],               # integrative relay
    5: ADDED_GROUP_COLORS["associative"],    # associative / learning
    6: ADDED_GROUP_COLORS["modulatory"],     # state-dependent / modulatory
    -1: ADDED_GROUP_COLORS["other"],
}
DROSOPHILA_GROUP_PANEL_ORDER = [
    "olfactory neuropil",
    "mushroom body",
    "optic neuropil",
    "central complex",
    "protocerebrum",
    "premotor interface",
    "other",
]


def anatomy_order_for_group_panel(species: str) -> list[str]:
    if species == "Drosophila":
        return DROSOPHILA_GROUP_PANEL_ORDER
    return fc_dynamics.ANATOMY_ORDER[species]


def apply_local_group_colors() -> None:
    fc_dynamics.ANATOMY_COLORS.update(ANATOMY_COLORS)
    post_dca_scatter.CLASS_COLORS.update(CLASS_COLORS)
    sc_predictors.ANATOMY_COLORS.update(ANATOMY_COLORS)
    sc_predictors.add_pairwise_sig_bars = add_permutation_holm_bars
    if hasattr(celegans, "rep") and hasattr(celegans.rep, "ANATOMY_COLORS"):
        celegans.rep.ANATOMY_COLORS.update(
            {key: ANATOMY_COLORS[key] for key in celegans.rep.ANATOMY_ORDER if key in ANATOMY_COLORS}
        )


def p_to_stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*"


def draw_significance_bars(
    ax: plt.Axes,
    significant_pairs: list[tuple[tuple[int, int], float]],
    y_start: float = SIG_BAR_Y_START,
    y_step: float = SIG_BAR_STEP,
) -> None:
    if not significant_pairs:
        return
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1.0
    step = y_range * y_step
    bar_h = step * SIG_BAR_HEIGHT_FRAC
    for level, ((i, j), p_value) in enumerate(significant_pairs):
        y = y_max + y_range * y_start + level * step
        ax.plot(
            [i, i, j, j],
            [y, y + bar_h, y + bar_h, y],
            lw=SIG_BAR_LW,
            c=SIG_BAR_COLOR,
            clip_on=False,
        )
        ax.text(
            (i + j) / 2,
            y + bar_h,
            p_to_stars(p_value),
            ha="center",
            va="bottom",
            fontsize=STAR_FS,
            color=SIG_BAR_COLOR,
            clip_on=False,
        )
    ax.set_ylim(y_min, y_max)


def significant_pair_list(
    pair_rows: list[dict[str, object]],
    alpha: float,
    max_bars: int | None = None,
) -> list[tuple[tuple[int, int], float]]:
    pairs = [
        ((int(row["group_i"]), int(row["group_j"])), float(row["holm_p"]))
        for row in pair_rows
        if np.isfinite(row["holm_p"]) and float(row["holm_p"]) < alpha
    ]
    pairs = sorted(pairs, key=lambda item: (item[1], item[0][1] - item[0][0]))
    if max_bars is not None:
        pairs = pairs[:max_bars]
    return sorted(pairs, key=lambda item: (item[0][1] - item[0][0], item[1]))


def add_permutation_holm_bars(
    ax: plt.Axes,
    values_by_group: list[np.ndarray],
    seed: int,
    alpha: float = 0.05,
) -> tuple[float, list[dict[str, object]]]:
    global_p, pair_rows = sc_predictors.permutation_group_tests(values_by_group, seed=seed)
    draw_significance_bars(ax, significant_pair_list(pair_rows, alpha))
    return global_p, pair_rows
PANEL_LABEL_X_OFFSETS = {
    "B": -0.030,
    "G": -0.030,
    "C": -0.032,
    "D": -0.040,
    "E": -0.040,
    "H": -0.032,
    "I": -0.040,
    "J": -0.040,
}


def remove_inherited_panel_labels(axes: list[plt.Axes]) -> None:
    for ax in axes:
        for text in list(ax.texts):
            if text.get_text() in INHERITED_PANEL_LABELS:
                text.remove()


def add_aligned_panel_labels(fig: plt.Figure, label_axes_by_row: list[list[plt.Axes]]) -> None:
    fig.canvas.draw()
    for axes, panel_labels in zip(label_axes_by_row, PANEL_LABELS_BY_ROW):
        row_y = max(
            ax.get_position().y1
            for ax, label in zip(axes, panel_labels)
            if label
        ) + PANEL_LABEL_Y_PAD
        for ax, label in zip(axes, panel_labels):
            if not label:
                continue
            pos = ax.get_position()
            dx = PANEL_LABEL_X_OFFSETS.get(label, -0.012)
            fig.text(
                pos.x0 + dx,
                row_y,
                label,
                ha="right",
                va="bottom",
                fontsize=PANEL_FS,
                fontweight="bold",
            )


def add_species_title(fig: plt.Figure, axes: list[plt.Axes], text: str) -> None:
    left = min(ax.get_position().x0 for ax in axes)
    right = max(ax.get_position().x1 for ax in axes)
    top = max(ax.get_position().y1 for ax in axes)
    fig.text(
        (left + right) / 2,
        top + SPECIES_TITLE_OFFSET_Y,
        text,
        ha="center",
        va="bottom",
        fontsize=AXIS_FS,
        fontweight="bold",
        style="italic" if text.startswith("C.") else "normal",
    )


def shift_axes_y(axes: list[plt.Axes], dy: float) -> None:
    for ax in axes:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0 + dy, pos.width, pos.height])


def shift_axes_x(axes: list[plt.Axes], dx: float) -> None:
    for ax in axes:
        pos = ax.get_position()
        ax.set_position([pos.x0 + dx, pos.y0, pos.width, pos.height])


def match_network_dca_label(ax: plt.Axes) -> None:
    for text in ax.texts:
        if text.get_text() == "Post-DCA":
            text.set_text(DCA_POST_LABEL)
        elif text.get_text() == "Pre-DCA":
            text.set_text(DCA_PRE_LABEL)


def style_network_axis_labels(ax: plt.Axes) -> None:
    for text in ax.texts:
        label = text.get_text()
        if label in {"high", "low"}:
            text.set_fontsize(TICK_FS)
            text.set_fontweight("normal")
            text.set_color("#222222")
        elif label in {"Post-DCA", "Pre-DCA", DCA_POST_LABEL, DCA_PRE_LABEL}:
            text.set_fontsize(AXIS_FS)
            text.set_fontweight("bold")
            text.set_color("#222222")


def style_measure_ticks(
    top_axes: list[plt.Axes],
    bottom_axes: list[plt.Axes],
) -> None:
    scatter_axes = [top_axes[1], bottom_axes[1]]
    group_axes = top_axes[2:] + bottom_axes[2:]
    for ax in scatter_axes:
        ax.tick_params(axis="both", labelsize=TICK_FS)
        ax.xaxis.label.set_size(AXIS_FS)
        ax.yaxis.label.set_size(AXIS_FS)
    for ax in group_axes:
        ax.tick_params(axis="y", labelsize=TICK_FS)
        ax.tick_params(axis="x", labelsize=GROUP_LABEL_FS)
        ax.xaxis.label.set_size(AXIS_FS)
        ax.yaxis.label.set_size(AXIS_FS)

    # Reduce tick density for the species scatter panels.
    top_axes[1].xaxis.set_major_locator(MaxNLocator(nbins=4))
    top_axes[1].yaxis.set_major_locator(MaxNLocator(nbins=4))
    bottom_axes[1].xaxis.set_major_locator(MaxNLocator(nbins=4))
    bottom_axes[1].yaxis.set_major_locator(MaxNLocator(nbins=4))


def draw_dca_scatter_panel(ax: plt.Axes, species: str, predictor: str, title: str) -> None:
    summary = post_dca_scatter.replace_fcv_with_recording_zscore(
        pd.read_csv(post_dca_scatter.TARGET_SUMMARY)
    )
    values = summary[summary["EdgeStdFCV"].notna()][["species", "node", "PostDCA", "EdgeStdFCV"]].copy()
    if predictor == "PreDCA":
        sc_values = pd.read_csv(SC_VALUES_CSV)
        pre = sc_values[["species", "node", "PreDCA"]].drop_duplicates(["species", "node"])
        values = values.drop(columns=["PreDCA"], errors="ignore").merge(pre, on=["species", "node"], how="left")
    stats = pd.DataFrame(
        [post_dca_scatter.trend_stats(values, species, predictor)]
    )
    if predictor == "PostDCA":
        post_dca_scatter.plot_scatter(ax, values, stats, species, "B")
    else:
        plot_dca_scatter(ax, values, stats, species, predictor)
    ax.set_title("")
    ax.set_xlabel(DCA_POST_LABEL if predictor == "PostDCA" else DCA_PRE_LABEL, fontsize=AXIS_FS)
    ax.set_ylabel("FCV", fontsize=AXIS_FS)
    for text in list(ax.texts):
        text.remove()
    stat = stats.iloc[0]
    p = float(stat["pearson_p"])
    p_text = "p<0.001" if p < 0.001 else f"p={p:.3f}"
    ax.text(
        0.05,
        0.95,
        f"r={float(stat['pearson_r']):.2f}\n{p_text}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=STAT_FS,
    )


def plot_dca_scatter(ax: plt.Axes, values: pd.DataFrame, stats: pd.DataFrame, species: str, predictor: str) -> None:
    sub = post_dca_scatter.add_functional_classes(
        values.loc[values["species"].eq(species), ["species", "node", predictor, "EdgeStdFCV"]]
    )
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=[predictor, "EdgeStdFCV"])
    x = sub[predictor].to_numpy(float)
    for class_order, group in sub.groupby("class_order", sort=True):
        key = int(class_order)
        ax.scatter(
            group[predictor].to_numpy(float),
            group["EdgeStdFCV"].to_numpy(float),
            s=28,
            marker=post_dca_scatter.CLASS_MARKERS.get(key, "o"),
            color=CLASS_COLORS.get(key, CLASS_COLORS[-1]),
            edgecolor="#202020",
            linewidth=0.35,
            alpha=0.88,
        )
    stat = stats[stats["species"].eq(species)].iloc[0]
    if len(x) >= 3 and np.nanstd(x) > 0:
        slope = float(stat["theilsen_slope"])
        intercept = float(stat["theilsen_intercept"])
        grid = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ax.plot(grid, slope * grid + intercept, color="black", lw=1.0)
    ax.set_xlabel(DCA_PRE_LABEL, fontsize=AXIS_FS)
    ax.set_ylabel("FCV", fontsize=AXIS_FS)
    ax.tick_params(labelsize=TICK_FS, width=0.7, length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_mannwhitney_holm_bars(
    ax: plt.Axes,
    values_by_group: list[np.ndarray],
    alpha: float = 0.05,
    max_bars: int | None = 6,
    reference_idx: int | None = None,
    reference_indices: list[int] | None = None,
    display_pairs: list[tuple[int, int]] | None = None,
    y_start: float = SIG_BAR_Y_START,
    y_step: float = SIG_BAR_STEP,
) -> list[dict[str, object]]:
    pairs = list(itertools.combinations(range(len(values_by_group)), 2))
    raw_p: list[float] = []
    for i, j in pairs:
        a = np.asarray(values_by_group[i], dtype=float)
        b = np.asarray(values_by_group[j], dtype=float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if len(a) < 2 or len(b) < 2:
            raw_p.append(np.nan)
            continue
        try:
            raw_p.append(fc_dynamics.stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except ValueError:
            raw_p.append(np.nan)
    holm_p = fc_dynamics.holm_adjust(raw_p)
    rows = [
        {"group_i": i, "group_j": j, "mannwhitney_p": raw, "holm_p": adj}
        for (i, j), raw, adj in zip(pairs, raw_p, holm_p, strict=False)
    ]
    display_rows = rows
    if display_pairs is not None:
        pair_set = {tuple(sorted(pair)) for pair in display_pairs}
        display_rows = [
            row
            for row in rows
            if tuple(sorted((int(row["group_i"]), int(row["group_j"])))) in pair_set
        ]
    elif reference_idx is not None:
        reference_indices = [reference_idx]
    if display_pairs is None and reference_indices is not None:
        ref_set = set(reference_indices)
        display_rows = [
            row
            for row in rows
            if int(row["group_i"]) in ref_set or int(row["group_j"]) in ref_set
        ]
    draw_significance_bars(
        ax,
        significant_pair_list(display_rows, alpha, max_bars=max_bars),
        y_start=y_start,
        y_step=y_step,
    )
    return rows


def draw_fcv_anatomy_panel(ax: plt.Axes, species: str, show_ylabel: bool = True) -> None:
    spec = next(spec for spec in fc_dynamics.SPECIES if spec["species"] == species)
    rec_df = fc_dynamics.load_species_recording_table(spec)
    order = [
        group
        for group in anatomy_order_for_group_panel(species)
        if group in set(rec_df["anatomy_group"])
    ]
    values_by_group: list[np.ndarray] = []
    rng = np.random.default_rng(20260530)
    for idx, group in enumerate(order):
        vals = (
            rec_df.loc[rec_df["anatomy_group"].eq(group), "EdgeStdFCV"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
        )
        values_by_group.append(vals)
        if len(vals):
            ax.scatter(
                np.full(len(vals), idx) + rng.normal(0, 0.070, size=len(vals)),
                vals,
                s=9.0,
                marker=fc_dynamics.ANATOMY_MARKERS.get(group, "o"),
                facecolor=ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]),
                edgecolor="#222222",
                alpha=0.1,
                linewidth=0.25,
                rasterized=True,
                zorder=3,
            )
    box = ax.boxplot(
        values_by_group,
        positions=np.arange(len(order)),
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 0.8},
        boxprops={"linewidth": 0.75, "edgecolor": "#222222"},
        whiskerprops={"linewidth": 0.65, "color": "#222222"},
        capprops={"linewidth": 0.65, "color": "#222222"},
    )
    for patch, group in zip(box["boxes"], order, strict=False):
        color = ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"])
        patch.set_facecolor(color)
        patch.set_alpha(0.34)
        patch.set_edgecolor(color)
        patch.set_linewidth(1.05)
        patch.set_zorder(1)

    display_pairs = None
    if species == "Drosophila":
        candidate_pairs = [
            ("mushroom body", "olfactory neuropil"),
            ("mushroom body", "optic neuropil"),
            ("mushroom body", "central complex"),
            ("mushroom body", "protocerebrum"),
            ("mushroom body", "premotor interface"),
            ("olfactory neuropil", "optic neuropil"),
            ("olfactory neuropil", "protocerebrum"),
        ]
        display_pairs = [
            (order.index(group_a), order.index(group_b))
            for group_a, group_b in candidate_pairs
            if group_a in order and group_b in order
        ]
    add_mannwhitney_holm_bars(
        ax,
        values_by_group,
        alpha=0.05,
        max_bars=None if display_pairs else 6,
        display_pairs=display_pairs,
        y_start=-0.15 if display_pairs else SIG_BAR_Y_START,
        y_step=0.060 if display_pairs else SIG_BAR_STEP,
    )
    ax.axhline(0, color="#777777", lw=0.55, ls=":", zorder=0)
    ax.set_title("")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=38, ha="right", fontsize=GROUP_LABEL_FS)
    for tick, group in zip(ax.get_xticklabels(), order, strict=False):
        tick.set_color(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
        tick.set_fontweight("bold")
    if show_ylabel:
        ax.set_ylabel("FCV", fontsize=AXIS_FS)
    else:
        ax.set_ylabel("FCV", fontsize=AXIS_FS)
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="both", direction="out", pad=1.3, length=2.2, width=0.55)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_sc_anatomy_panel(ax: plt.Axes, species: str, measure: str, title: str, show_ylabel: bool) -> None:
    if SC_VALUES_CSV.exists():
        values = pd.read_csv(SC_VALUES_CSV)
    else:
        values = sc_predictors.load_values()
    values = values.loc[values["species"].eq(species)].copy()
    order = [
        group
        for group in anatomy_order_for_group_panel(species)
        if group in set(values["anatomy_group"].dropna())
    ]

    values_by_group: list[np.ndarray] = []
    rng_seed = sum(ord(ch) for ch in f"sc-{species}-{measure}") % (2**32)
    rng = np.random.default_rng(rng_seed)
    for idx, group in enumerate(order):
        vals = (
            values.loc[values["anatomy_group"].eq(group), measure]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
        )
        values_by_group.append(vals)
        if len(vals):
            ax.scatter(
                np.full(len(vals), idx) + rng.normal(0, 0.070, size=len(vals)),
                vals,
                s=9.0,
                marker=fc_dynamics.ANATOMY_MARKERS.get(group, "o"),
                facecolor=ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]),
                edgecolor="#222222",
                alpha=0.58,
                linewidth=0.25,
                rasterized=True,
                zorder=3,
            )

    box = ax.boxplot(
        values_by_group,
        positions=np.arange(len(order)),
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 0.8},
        boxprops={"linewidth": 0.75, "edgecolor": "#222222"},
        whiskerprops={"linewidth": 0.65, "color": "#222222"},
        capprops={"linewidth": 0.65, "color": "#222222"},
    )
    for patch, group in zip(box["boxes"], order, strict=False):
        color = ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"])
        patch.set_facecolor(color)
        patch.set_alpha(0.34)
        patch.set_edgecolor(color)
        patch.set_linewidth(1.05)
        patch.set_zorder(1)

    if (species, measure) in sc_predictors.GROUP_Y_LIMITS:
        ax.set_ylim(*sc_predictors.GROUP_Y_LIMITS[(species, measure)])
    perm_seed = sum(ord(ch) for ch in f"perm-{species}-{measure}") % (2**32)
    add_permutation_holm_bars(ax, values_by_group, seed=perm_seed, alpha=0.05)
    ax.axhline(0, color="#777777", lw=0.55, ls=":", zorder=0)
    ax.set_title("")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=38, ha="right", fontsize=GROUP_LABEL_FS)
    for tick, group in zip(ax.get_xticklabels(), order, strict=False):
        tick.set_color(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
        tick.set_fontweight("bold")
    ax.set_ylabel(title, fontsize=AXIS_FS)
    ax.tick_params(axis="both", direction="out", pad=1.3, length=2.2, width=0.55)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_drosophila_ncomms_network(ax: plt.Axes) -> None:
    values = pd.read_csv(SC_VALUES_CSV)
    values = values[
        values["species"].eq("Drosophila")
        & values["PreDCA"].notna()
        & values["anatomy_group"].notna()
    ][["node", "PreDCA", "anatomy_group"]].drop_duplicates("node")

    matched = drosophila.matched_node_table()
    sc = drosophila.aggregate_sc_to_matched_nodes(matched)
    nodes = matched[["atlas_region"]].drop_duplicates().rename(columns={"atlas_region": "node"})
    nodes = nodes.merge(values, on="node", how="inner")
    present_groups = set(nodes["anatomy_group"].dropna())
    order = [
        group
        for group in anatomy_order_for_group_panel("Drosophila")
        if group != "central complex" and group in present_groups
    ]
    nodes["anatomy_order"] = nodes["anatomy_group"].map({name: idx for idx, name in enumerate(order)})
    nodes = nodes.dropna(subset=["anatomy_order"]).copy()
    nodes["anatomy_order"] = nodes["anatomy_order"].astype(int)
    nodes = nodes.sort_values(["anatomy_order", "PreDCA", "node"], ascending=[True, False, True]).reset_index(drop=True)

    regions = nodes["node"].astype(str).tolist()
    sc = sc.loc[regions, regions].copy()

    predca = nodes["PreDCA"].to_numpy(float)
    x_min, x_max = float(np.nanmin(predca)), float(np.nanmax(predca))
    denom = x_max - x_min if x_max > x_min else 1.0
    x_scaled = 0.06 + 0.96 * (x_max - predca) / denom
    y_base = np.array([len(order) - 1 - int(idx) for idx in nodes["anatomy_order"]], dtype=float)
    y_jitter = np.zeros(len(nodes), dtype=float)
    for group, idxs in nodes.groupby("anatomy_group", sort=False).groups.items():
        idxs = np.asarray(list(idxs), dtype=int)
        y_jitter[idxs] = np.linspace(-0.30, 0.30, len(idxs)) if len(idxs) > 1 else 0.0
    y_scaled = y_base + y_jitter
    x_pos = pd.Series(x_scaled, index=regions)
    y_pos = pd.Series(y_scaled, index=regions)

    edge_rows = []
    arr = sc.to_numpy(dtype=float)
    for i, src in enumerate(regions):
        row = arr[i].copy()
        row[i] = 0
        nz = np.where(row > 0)[0]
        if len(nz) > 12:
            keep = nz[np.argsort(row[nz])[-12:]]
        else:
            keep = nz
        for j in keep:
            edge_rows.append((src, regions[j], row[j]))
    max_edge = max([edge[2] for edge in edge_rows], default=1.0)
    for src, dst, weight in edge_rows:
        x1, y1 = float(x_pos[src]), float(y_pos[src])
        x2, y2 = float(x_pos[dst]), float(y_pos[dst])
        if abs(x1 - x2) < 0.02 and abs(y1 - y2) < 0.02:
            continue
        ax.add_patch(
            FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                arrowstyle="-|>",
                mutation_scale=4.8,
                linewidth=0.14 + 0.32 * np.log1p(weight) / np.log1p(max_edge),
                color="#3f3f3f",
                alpha=0.16,
                shrinkA=3.0,
                shrinkB=3.0,
                zorder=1,
            )
        )

    for group in order:
        sub = nodes[nodes["anatomy_group"].eq(group)]
        if sub.empty:
            continue
        xy = np.asarray([[x_pos[node], y_pos[node]] for node in sub["node"]])
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=24,
            c=ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]),
            edgecolor="white",
            linewidth=0.45,
            alpha=0.94,
            zorder=3,
            label=group,
        )

    y_ticks = [len(order) - 1 - i for i in range(len(order))]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(order, fontsize=NETWORK_LABEL_FS)
    for tick, group in zip(ax.get_yticklabels(), order, strict=False):
        tick.set_color(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
        tick.set_fontweight("bold")

    ax.set_xlim(-0.06, 1.08)
    ax.set_ylim(-0.72, len(order) - 0.28)
    arrow_y = -0.50
    arrow_vertices = [
        (0.08, arrow_y - 0.042),
        (0.86, arrow_y - 0.010),
        (0.86, arrow_y - 0.034),
        (0.94, arrow_y),
        (0.86, arrow_y + 0.034),
        (0.86, arrow_y + 0.010),
        (0.08, arrow_y + 0.042),
    ]
    ax.add_patch(Polygon(arrow_vertices, closed=True, facecolor="#5b5b5b", edgecolor="none", alpha=0.82, clip_on=False, zorder=6))
    ax.text(0.02, arrow_y - 0.15, "high", ha="left", va="top", fontsize=TICK_FS, color="#222222")
    ax.text(0.92, arrow_y - 0.15, "low", ha="left", va="top", fontsize=TICK_FS, color="#222222")
    ax.text(0.42, arrow_y - 0.24, r"$\mathrm{DCA}_{\mathrm{pre}}$", ha="left", va="top", fontsize=AXIS_FS, color="#222222", fontweight="bold")
    style_network_axis_labels(ax)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)


def draw_celegans_row(fig: plt.Figure, row: dict) -> tuple[list[plt.Axes], list[plt.Axes]]:
    celegans.configure_shared_style()
    _, _, _, high_example, low_example, network_nodes, network_sc = celegans.load_representation_examples()

    ax_network = fig.add_subplot(row["network"])
    celegans.rep.add_network_panel(ax_network, network_nodes, network_sc, high_example, low_example)
    match_network_dca_label(ax_network)
    style_network_axis_labels(ax_network)
    celegans.add_panel_label(ax_network, "A", x=-0.10)

    ax_post_scatter = fig.add_subplot(row["dca_scatter"])
    draw_dca_scatter_panel(ax_post_scatter, "C. elegans", "PostDCA", "Post-DCA vs FCV")

    ax_fcv = fig.add_subplot(row["fcv"])
    draw_fcv_anatomy_panel(ax_fcv, "C. elegans", show_ylabel=True)

    ax_post = fig.add_subplot(row["post_dca"])
    draw_sc_anatomy_panel(ax_post, "C. elegans", "PostDCA", DCA_POST_LABEL, show_ylabel=True)

    ax_pre = fig.add_subplot(row["pre_dca"])
    draw_sc_anatomy_panel(ax_pre, "C. elegans", "PreDCA", DCA_PRE_LABEL, show_ylabel=False)

    axes = [ax_network, ax_post_scatter, ax_fcv, ax_post, ax_pre]
    remove_inherited_panel_labels(axes)
    return axes, axes


def draw_drosophila_row(fig: plt.Figure, row: dict) -> tuple[list[plt.Axes], list[plt.Axes]]:
    drosophila.set_style()

    ax_network = fig.add_subplot(row["network"])
    draw_drosophila_ncomms_network(ax_network)

    ax_post_scatter = fig.add_subplot(row["dca_scatter"])
    draw_dca_scatter_panel(ax_post_scatter, "Drosophila", "PreDCA", "Pre-DCA vs FCV")

    ax_fcv = fig.add_subplot(row["fcv"])
    draw_fcv_anatomy_panel(ax_fcv, "Drosophila", show_ylabel=True)

    ax_post = fig.add_subplot(row["post_dca"])
    draw_sc_anatomy_panel(ax_post, "Drosophila", "PostDCA", DCA_POST_LABEL, show_ylabel=True)

    ax_pre = fig.add_subplot(row["pre_dca"])
    draw_sc_anatomy_panel(ax_pre, "Drosophila", "PreDCA", DCA_PRE_LABEL, show_ylabel=False)

    axes = [ax_network, ax_post_scatter, ax_fcv, ax_post, ax_pre]
    remove_inherited_panel_labels(axes)
    return axes, axes


def prepare_layout() -> tuple[plt.Figure, list[dict]]:
    fs.apply_main_figure_style()
    fig = plt.figure(figsize=(16, 6.6))
    grid = fig.add_gridspec(
        2,
        5,
        left=0.075,
        right=0.985,
        top=0.935,
        bottom=0.105,
        width_ratios=[1.45, 0.92, 1.0, 1.0, 1.0],
        height_ratios=[1.0, 1.0],
        hspace=0.88,
        wspace=0.40,
    )
    rows = [
        {
            "network": grid[row, 0],
            "dca_scatter": grid[row, 1],
            "fcv": grid[row, 2],
            "post_dca": grid[row, 3],
            "pre_dca": grid[row, 4],
        }
        for row in range(2)
    ]
    return fig, rows


def save_figure(fig: plt.Figure) -> None:
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight", pad_inches=0.04, transparent=False)


def main() -> None:
    apply_local_group_colors()
    fig, rows = prepare_layout()
    top_axes, top_label_axes = draw_celegans_row(fig, rows[0])
    bottom_axes, bottom_label_axes = draw_drosophila_row(fig, rows[1])
    fig.canvas.draw()
    shift_axes_y(bottom_axes, BOTTOM_ROW_SHIFT_Y)
    shift_axes_x([top_axes[1], bottom_axes[1]], DCA_SCATTER_SHIFT_X)
    shift_axes_x([top_axes[2], bottom_axes[2]], FCV_PANEL_SHIFT_X)
    style_measure_ticks(top_axes, bottom_axes)
    fig.canvas.draw()
    add_aligned_panel_labels(fig, [top_label_axes, bottom_label_axes])
    add_species_title(fig, top_axes, "C. elegans")
    add_species_title(fig, bottom_axes, "Drosophila")
    save_figure(fig)
    plt.close(fig)
    print(f"Saved PNG: {OUT_PNG}")


if __name__ == "__main__":
    main()
