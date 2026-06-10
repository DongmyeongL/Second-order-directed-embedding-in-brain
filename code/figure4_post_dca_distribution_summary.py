"""Composite edge target-DCA distributions and Post-DCA versus FCV panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr, theilslopes

import figure_style as fs

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "final_summary_tables"
FIG = ROOT / "figures"

TARGET_SAMPLE = TAB / "edge_target_dca_distribution_sample_by_unit.csv"
TARGET_SUMMARY = TAB / "edge_target_dca_distribution_summary_by_unit.csv"
CLASS_TABLE = TAB / "current_plot_functional_classes.csv"
ZSCORE_FCV_TABLE = TAB / "highpass_ce_zf_plot_measures_recording_zscore_node_summary.csv"
OUT = FIG / "edge_target_dca_distribution_and_postdca_fcv_summary.png"
OUT_STATS = TAB / "edge_target_dca_distribution_and_postdca_fcv_summary_stats.csv"

SPECIES = ["C. elegans", "Drosophila", "Zebrafish"]
COLORS = fs.SPECIES_COLORS.copy()
YRANGE_PERCENTILES = {
    "C. elegans": (5, 95),
    "Drosophila": (8, 92),
    "Zebrafish": (8, 92),
}
Y_LIMITS = {
    "Drosophila": (-0.08, 0.02),
}
CLASS_COLORS = fs.FUNCTIONAL_GROUP_COLORS.copy()
CLASS_LABELS = {
    0: "olfactory / chemosensory",
    1: "visual / optic",
    2: "other sensory",
    3: "sensorimotor output",
    4: "integrative relay",
    5: "associative / learning",
    6: "state-dependent / modulatory",
    -1: "unclassified",
}
CLASS_MARKERS = {
    0: "o",
    1: "^",
    2: "v",
    3: "s",
    4: "D",
    5: "P",
    6: "X",
    -1: "o",
}


def p_text(p: float) -> str:
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}"


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.012, y: float = 1.02) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11,
        fontweight="bold",
    )


def add_functional_classes(df: pd.DataFrame) -> pd.DataFrame:
    classes = pd.read_csv(CLASS_TABLE)
    keep = ["species", "node", "class_order", "class_label"]
    out = df.merge(classes[keep], on=["species", "node"], how="left")
    out["class_order"] = out["class_order"].fillna(-1).astype(int)
    out["class_label"] = out["class_label"].fillna("unclassified")
    return out


def replace_fcv_with_recording_zscore(df: pd.DataFrame) -> pd.DataFrame:
    z = pd.read_csv(ZSCORE_FCV_TABLE)[["species", "node", "EdgeStdFCV"]].rename(
        columns={"EdgeStdFCV": "EdgeStdFCV_zscore"}
    )
    out = df.drop(columns=["EdgeStdFCV_zscore"], errors="ignore").merge(z, on=["species", "node"], how="left")
    if "EdgeStdFCV" in out.columns:
        out = out.rename(columns={"EdgeStdFCV": "EdgeStdFCV_raw"})
    out["EdgeStdFCV"] = out["EdgeStdFCV_zscore"]
    out["fcv_normalization"] = "within_recording_zscore_then_node_mean"
    return out


def permutation_p(x: np.ndarray, y: np.ndarray, statistic: str, n_perm: int = 1000) -> float:
    rng = np.random.default_rng(20260525)
    if statistic == "pearson":
        observed = pearsonr(x, y).statistic
        stat_fn = lambda xx, yy: pearsonr(xx, yy).statistic
    elif statistic == "spearman":
        observed = spearmanr(x, y).statistic
        stat_fn = lambda xx, yy: spearmanr(xx, yy).statistic
    else:
        raise ValueError(f"Unknown statistic: {statistic}")
    null = np.empty(n_perm, dtype=float)
    for idx in range(n_perm):
        null[idx] = stat_fn(x, rng.permutation(y))
    return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1))


def trend_stats(values: pd.DataFrame, species: str, predictor: str) -> dict[str, float | int | str]:
    sub = (
        values.loc[values["species"].eq(species), [predictor, "EdgeStdFCV"]]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )
    x = sub[predictor].to_numpy(float)
    y = sub["EdgeStdFCV"].to_numpy(float)
    pr = pearsonr(x, y)
    sr = spearmanr(x, y)
    kt = kendalltau(x, y)
    ols_slope, ols_intercept = np.polyfit(x, y, 1)
    ts_slope, ts_intercept, ts_low, ts_high = theilslopes(y, x, alpha=0.95)
    return {
        "species": species,
        "predictor": predictor,
        "response": "zEdgeStdFCV",
        "n": int(len(x)),
        "pearson_r": float(pr.statistic),
        "pearson_p": float(pr.pvalue),
        "pearson_perm_p": permutation_p(x, y, "pearson"),
        "spearman_rho": float(sr.statistic),
        "spearman_p": float(sr.pvalue),
        "spearman_perm_p": permutation_p(x, y, "spearman"),
        "kendall_tau": float(kt.statistic),
        "kendall_p": float(kt.pvalue),
        "ols_slope": float(ols_slope),
        "ols_intercept": float(ols_intercept),
        "theilsen_slope": float(ts_slope),
        "theilsen_intercept": float(ts_intercept),
        "theilsen_slope_ci_low": float(ts_low),
        "theilsen_slope_ci_high": float(ts_high),
        "r_squared": float(pr.statistic * pr.statistic),
    }


def distribution_limits(sub: pd.DataFrame, samp: pd.DataFrame, species: str) -> tuple[float, float]:
    vals_all = samp["target_DCA"].to_numpy(float)
    vals_all = vals_all[np.isfinite(vals_all)]
    if species in Y_LIMITS:
        return Y_LIMITS[species]
    lo_p, hi_p = YRANGE_PERCENTILES[species]
    ymin, ymax = np.nanpercentile(vals_all, [lo_p, hi_p])
    anchors = sub[["PostDCA", "edge_target_dca_q25", "edge_target_dca_q75"]].to_numpy(float).ravel()
    anchors = anchors[np.isfinite(anchors)]
    if len(anchors):
        ymin = min(ymin, float(np.nanpercentile(anchors, 1)))
        ymax = max(ymax, float(np.nanpercentile(anchors, 99)))
    ypad = max((ymax - ymin) * 0.08, 0.0015)
    return ymin - ypad, ymax + ypad


def plot_distribution(
    ax_bar: plt.Axes,
    ax: plt.Axes,
    summary: pd.DataFrame,
    sample: pd.DataFrame,
    species: str,
    panel_label: str,
) -> None:
    sub = summary[summary["species"].eq(species) & summary["EdgeStdFCV"].notna()].copy()
    sub = add_functional_classes(sub)
    sub = sub.sort_values("EdgeStdFCV", ascending=True).reset_index(drop=True)
    sub["x"] = np.arange(len(sub), dtype=float)
    x_lookup = dict(zip(sub["node"], sub["x"]))
    class_lookup = dict(zip(sub["node"], sub["class_order"]))
    samp = sample[sample["species"].eq(species) & sample["node"].isin(x_lookup)].copy()
    samp["x"] = samp["node"].map(x_lookup)
    samp["class_order"] = samp["node"].map(class_lookup).fillna(-1).astype(int)

    ymin, ymax = distribution_limits(sub, samp, species)
    fcv = sub["EdgeStdFCV"].to_numpy(float)
    fcv_norm = (fcv - np.nanmin(fcv)) / max(np.nanmax(fcv) - np.nanmin(fcv), 1e-12)

    ax_bar.imshow(
        fcv_norm[np.newaxis, :],
        aspect="auto",
        extent=(-0.5, len(sub) - 0.5, 0, 1),
        cmap="viridis",
        interpolation="nearest",
    )
    ax_bar.set_yticks([])
    ax_bar.set_xticks([])
    ax_bar.set_ylabel("z-FCV", fontsize=7, rotation=0, ha="right", va="center", labelpad=12)
    ax_bar.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax_bar.text(-0.5, 1.22, "low", fontsize=6.5, ha="left", va="bottom")
    ax_bar.text(len(sub) - 0.5, 1.22, "high", fontsize=6.5, ha="right", va="bottom")
    add_panel_label(ax_bar, panel_label, x=-0.012, y=1.45)

    plot_samp = samp.dropna(subset=["x", "target_DCA"]).copy()
    plot_samp = plot_samp[(plot_samp["target_DCA"] >= ymin) & (plot_samp["target_DCA"] <= ymax)]
    if len(plot_samp) > 90000:
        plot_samp = plot_samp.sample(90000, random_state=20260524)
    jitter_sd = 0.055 if len(sub) <= 60 else 0.035
    jitter = np.random.default_rng(20260524).normal(0, jitter_sd, size=len(plot_samp))
    plot_samp["x_jitter"] = plot_samp["x"].to_numpy(float) + jitter
    for class_order, group in plot_samp.groupby("class_order", sort=True):
        key = int(class_order)
        ax.scatter(
            group["x_jitter"].to_numpy(float),
            group["target_DCA"].to_numpy(float),
            s=4.0 if species == "C. elegans" else 2.0,
            marker=CLASS_MARKERS.get(key, "o"),
            color=CLASS_COLORS.get(key, CLASS_COLORS[-1]),
            alpha=0.22 if species == "C. elegans" else 0.14,
            linewidth=0,
            zorder=3,
        )

    box_half_width = 0.065 if len(sub) <= 60 else 0.035
    mean_x: list[float] = []
    mean_y: list[float] = []
    for _, row in sub.iterrows():
        x0 = float(row["x"])
        q25 = np.clip(row["edge_target_dca_q25"], ymin, ymax)
        q75 = np.clip(row["edge_target_dca_q75"], ymin, ymax)
        ax.add_patch(
            plt.Rectangle(
                (x0 - box_half_width, q25),
                2 * box_half_width,
                max(q75 - q25, 1e-9),
                facecolor="white",
                edgecolor="#111111",
                linewidth=0.45,
                alpha=0.72,
                zorder=2,
            )
        )
        mean_value = row["PostDCA"]
        if np.isfinite(mean_value):
            mean_x.append(x0)
            mean_y.append(float(np.clip(mean_value, ymin, ymax)))
            ax.scatter(
                x0,
                np.clip(mean_value, ymin, ymax),
                s=26,
                marker="D",
                color=fs.MAIN_COLORS["post_out"],
                edgecolor="white",
                linewidth=0.45,
                zorder=6,
            )
    if mean_x:
        ax.plot(mean_x, mean_y, color=fs.MAIN_COLORS["post_out"], lw=1.15, alpha=0.9, zorder=5)

    unit = "neurons" if species == "C. elegans" else "regions"
    ax.axhline(0, color="#777777", lw=0.7, ls="--")
    ax.set_title(f"{species}: outgoing edge target-DCA distribution", fontsize=10, pad=2)
    ax.set_xlabel(f"Source {unit} sorted by z-FCV", fontsize=8.5)
    ax.set_ylabel("Target DCA", fontsize=8.5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-0.75, len(sub) - 0.25)
    font_size = 3.0 if len(sub) > 120 else 5.0 if len(sub) > 60 else 6.1
    ax.set_xticks(sub["x"])
    ax.set_xticklabels(sub["node"], fontsize=font_size, rotation=90)
    for tick, class_order in zip(ax.get_xticklabels(), sub["class_order"], strict=False):
        tick.set_color(CLASS_COLORS.get(int(class_order), CLASS_COLORS[-1]))
    ax.tick_params(axis="y", labelsize=7, width=0.7, length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_scatter(ax: plt.Axes, values: pd.DataFrame, stats: pd.DataFrame, species: str, panel_label: str) -> None:
    sub = add_functional_classes(values.loc[values["species"].eq(species), ["species", "node", "PostDCA", "EdgeStdFCV"]])
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["PostDCA", "EdgeStdFCV"])
    x = sub["PostDCA"].to_numpy(float)
    y = sub["EdgeStdFCV"].to_numpy(float)

    for class_order, group in sub.groupby("class_order", sort=True):
        key = int(class_order)
        ax.scatter(
            group["PostDCA"].to_numpy(float),
            group["EdgeStdFCV"].to_numpy(float),
            s=28,
            marker=CLASS_MARKERS.get(key, "o"),
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

    ax.text(
        0.05,
        0.95,
        "n={n}\nrho={rho:.2f}\nTS slope={slope:.3f}\n95% CI [{lo:.3f}, {hi:.3f}]".format(
            n=int(stat["n"]),
            rho=float(stat["spearman_rho"]),
            slope=float(stat["theilsen_slope"]),
            lo=float(stat["theilsen_slope_ci_low"]),
            hi=float(stat["theilsen_slope_ci_high"]),
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.set_title(species, fontsize=10, pad=3)
    ax.set_xlabel("Post-DCA", fontsize=9)
    ax.set_ylabel("z-FCV", fontsize=9)
    ax.tick_params(labelsize=8, width=0.7, length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    add_panel_label(ax, panel_label, x=-0.13, y=1.08)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    sample = pd.read_csv(TARGET_SAMPLE)
    summary = replace_fcv_with_recording_zscore(pd.read_csv(TARGET_SUMMARY))
    values = summary[summary["EdgeStdFCV"].notna()][["species", "node", "PostDCA", "EdgeStdFCV"]].copy()
    stats = pd.DataFrame([trend_stats(values, species, "PostDCA") for species in SPECIES])
    stats.to_csv(OUT_STATS, index=False)
    fig = plt.figure(figsize=(11.2, 13.0), constrained_layout=True)
    gs = fig.add_gridspec(
        4,
        3,
        height_ratios=[1.23, 1.23, 1.23, 1.0],
        hspace=0.10,
        wspace=0.28,
    )

    for row, (species, panel_label) in enumerate(zip(SPECIES, ["a", "b", "c"])):
        top_gs = gs[row, :].subgridspec(2, 1, height_ratios=[0.18, 1.0], hspace=0.02)
        ax_bar = fig.add_subplot(top_gs[0, 0])
        ax = fig.add_subplot(top_gs[1, 0], sharex=ax_bar)
        plot_distribution(ax_bar, ax, summary, sample, species, panel_label)

    for col, (species, panel_label) in enumerate(zip(SPECIES, ["d", "e", "f"])):
        ax = fig.add_subplot(gs[3, col])
        plot_scatter(ax, values, stats, species, panel_label)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker=CLASS_MARKERS[key],
            linestyle="none",
            markerfacecolor=CLASS_COLORS[key],
            markeredgecolor="#202020",
            markeredgewidth=0.45,
            markersize=4.8,
            label=CLASS_LABELS[key],
        )
        for key in [0, 1, 2, 3, 4, 5, 6, -1]
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.018),
        ncol=4,
        frameon=False,
        fontsize=6.2,
        handletextpad=0.4,
        columnspacing=1.0,
    )

    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {OUT}")
    print(f"Saved stats: {OUT_STATS}")


if __name__ == "__main__":
    main()
