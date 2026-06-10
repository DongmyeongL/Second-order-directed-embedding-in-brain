"""Figure 7 synthesis: functional-group structure-function signature.

This final-figure variant keeps only the functional-group signature panel:
for each species and functional group, it summarizes FCV, Post-DCA, and
Pre-DCA on a within-species relative scale.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

import figure_style as fs


ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "source_inputs" / "ncomms_tables"
BUILD = ROOT
BUILD_FIG = BUILD / "figures"
BUILD_TAB = BUILD / "data" / "final_summary_tables"
BUILD_FIG.mkdir(parents=True, exist_ok=True)
BUILD_TAB.mkdir(parents=True, exist_ok=True)

FCV_TABLE = TAB / "figure7_recording_functional_group_points.csv"
FCV_RECORDING_NODE_TABLE = TAB / "highpass_ce_zf_plot_measures_recording_node.csv"
FCV_NODE_GROUP_TABLE = TAB / "figure7_comparative_fcv_fine_class_node_map.csv"
POST_TABLE = TAB / "figure7_postdca_functional_group_node_values.csv"
PRE_TABLE = TAB / "figure7_predca_functional_group_node_values.csv"
ZEBRAFISH_SUBJECT_DCA_TABLE = (
    ROOT
    / "data"
    / "source_inputs"
    / "external_processed"
    / "fcv_postdca_raw_recompute"
    / "out_data"
    / "zebrafish"
    / "post_dca_rank1"
    / "zebrafish_rank1_subject_region_post_dca.csv"
)

OUT_FIG = BUILD_FIG / "figure7_structure_function_synthesis_effects.png"
OUT_EFFECTS = BUILD_TAB / "figure7_structure_function_synthesis_planned_effects.csv"
OUT_GROUP = BUILD_TAB / "figure7_structure_function_synthesis_group_signature.csv"
OUT_GROUP_PERM = BUILD_TAB / "figure7_structure_function_synthesis_group_permutation_tests.csv"

STAT_CACHE_COLUMNS = {
    "effects": {
        "species",
        "metric",
        "mean_diff",
        "cohens_d",
        "permutation_two_sided_p",
    },
    "signature": {
        "species",
        "shared_fine_order",
        "metric",
        "leave_one_group_out_z_delta",
    },
    "group_perm": {
        "species",
        "metric",
        "shared_fine_order",
        "p_abs_permutation",
        "q_fdr_within_species_metric",
    },
}

SPECIES = ["C. elegans", "Drosophila", "Zebrafish"]
SPECIES_COLORS = fs.SPECIES_COLORS.copy()
METRIC_COLORS = fs.METRIC_COLORS.copy()
GROUP_LABELS = {
    0: "olf./chemo",
    1: "visual",
    2: "other sensory",
    3: "sensorimotor",
    4: "integrative",
    5: "assoc./learning",
    6: "state-mod.",
}
GROUP_COLORS = fs.FUNCTIONAL_GROUP_COLORS.copy()
HIGHLIGHT_GROUPS = {0, 5}
DISPLAY_GROUP_ORDER = [0, 5, 1, 2, 3, 4, 6]

CONTRASTS = {
    "C. elegans": {
        "high": [0],
        "low": [4],
        "label": "olf./chemo > integrative",
        "interpretation": "sensory-biased compact routing",
    },
    "Drosophila": {
        "high": [5],
        "low": [1, 4],
        "label": "assoc./learning > visual+integrative",
        "interpretation": "input-side MB axis",
    },
    "Zebrafish": {
        "high": [5],
        "low": [1, 3, 4],
        "label": "assoc./learning > visual+sensorimotor+integrative",
        "interpretation": "output-side telencephalic cascade",
    },
}


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.linewidth": 0.45,
        "axes.labelsize": 6.5,
    "axes.titlesize": 7.0,
        "xtick.labelsize": 6.0,
        "ytick.labelsize": 6.0,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "legend.fontsize": 6.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 300,
    }
)


def mean_ci95(values: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    mean = float(np.nanmean(values))
    if len(values) < 2:
        return mean, mean, mean
    sem = float(np.nanstd(values, ddof=1) / np.sqrt(len(values)))
    half_width = 1.96 * sem
    return mean, mean - half_width, mean + half_width


def format_p(p: float) -> str:
    if not np.isfinite(p):
        return "n/a"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def metric_tables() -> dict[str, pd.DataFrame]:
    # Match the scatter points used in figure7_functional_group_fcv_synthesis.png:
    # one point per recording x functional group for FCV.
    fcv = pd.read_csv(FCV_TABLE).rename(columns={"mean_fcv": "value"})
    fcv["metric"] = "FCV"
    fcv = fcv[["species", "recording_id", "shared_fine_order", "shared_fine_label", "value", "metric"]]

    post = pd.read_csv(POST_TABLE).rename(columns={"postdca_for_plot": "value"})
    post["metric"] = "Post-DCA"
    post = post[["species", "node", "shared_fine_order", "shared_fine_label", "value", "metric"]]
    post = use_zebrafish_subject_dca_points(post, "Post-DCA")

    pre = pd.read_csv(PRE_TABLE).rename(columns={"predca_for_plot": "value"})
    pre["metric"] = "Pre-DCA"
    pre = pre[["species", "node", "shared_fine_order", "shared_fine_label", "value", "metric"]]
    pre = use_zebrafish_subject_dca_points(pre, "Pre-DCA")
    return {"FCV": fcv, "Post-DCA": post, "Pre-DCA": pre}


def use_zebrafish_subject_dca_points(table: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Use the same zebrafish DCA scatter points as the original Figure 7.

    For C. elegans and Drosophila, the plotted points are node/region-level
    rows from the node-value tables. For zebrafish, the original synthesis
    figure plotted subject x region DCA estimates; use those rows for both
    the bar summary and permutation tests.
    """
    if not ZEBRAFISH_SUBJECT_DCA_TABLE.exists():
        return table
    source_col = {"Post-DCA": "Rank1PostDCA", "Pre-DCA": "Rank1PreDCA"}.get(metric)
    if source_col is None:
        return table
    z_classes = (
        table.loc[table["species"].eq("Zebrafish"), ["node", "shared_fine_order", "shared_fine_label"]]
        .drop_duplicates("node")
        .copy()
    )
    if z_classes.empty:
        return table
    z_subject = pd.read_csv(ZEBRAFISH_SUBJECT_DCA_TABLE).rename(
        columns={"Region": "node", source_col: "value", "Subject": "subject_id"}
    )
    z_subject = z_subject.merge(z_classes, on="node", how="inner")
    z_subject["species"] = "Zebrafish"
    z_subject["metric"] = metric
    z_subject = z_subject[
        ["species", "subject_id", "node", "shared_fine_order", "shared_fine_label", "value", "metric"]
    ].dropna(subset=["value", "shared_fine_order"])
    non_z = table.loc[~table["species"].eq("Zebrafish")].copy()
    return pd.concat([non_z, z_subject], ignore_index=True, sort=False)


def two_group_arrays(df: pd.DataFrame, species: str, high_orders: list[int], low_orders: list[int]) -> tuple[np.ndarray, np.ndarray]:
    sub = df[df["species"].eq(species)].copy()
    high = sub[sub["shared_fine_order"].isin(high_orders)]["value"].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
    low = sub[sub["shared_fine_order"].isin(low_orders)]["value"].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
    return high, low


def mean_diff(high: np.ndarray, low: np.ndarray) -> float:
    if len(high) == 0 or len(low) == 0:
        return np.nan
    return float(np.nanmean(high) - np.nanmean(low))


def cohens_d(high: np.ndarray, low: np.ndarray) -> float:
    if len(high) < 2 or len(low) < 2:
        return np.nan
    pooled = np.sqrt(((len(high) - 1) * np.nanvar(high, ddof=1) + (len(low) - 1) * np.nanvar(low, ddof=1)) / (len(high) + len(low) - 2))
    if not np.isfinite(pooled) or pooled == 0:
        return np.nan
    return float((np.nanmean(high) - np.nanmean(low)) / pooled)


def bootstrap_ci(
    high: np.ndarray,
    low: np.ndarray,
    statistic=mean_diff,
    n_boot: int = 5000,
    seed: int = 0,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    if len(high) == 0 or len(low) == 0:
        return np.nan, np.nan
    boot = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        h = rng.choice(high, size=len(high), replace=True)
        l = rng.choice(low, size=len(low), replace=True)
        boot[i] = statistic(h, l)
    return float(np.nanpercentile(boot, 2.5)), float(np.nanpercentile(boot, 97.5))


def permutation_p(high: np.ndarray, low: np.ndarray, n_perm: int = 5000, seed: int = 0) -> tuple[float, float]:
    if len(high) == 0 or len(low) == 0:
        return np.nan, np.nan
    values = np.concatenate([high, low])
    labels = np.concatenate([np.ones(len(high), dtype=int), np.zeros(len(low), dtype=int)])
    observed = float(np.nanmean(high) - np.nanmean(low))
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        p = rng.permutation(labels)
        null[i] = np.nanmean(values[p == 1]) - np.nanmean(values[p == 0])
    p_greater = float((np.sum(null >= observed) + 1) / (n_perm + 1))
    p_two = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1))
    return p_greater, p_two


def compute_effects(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for species in SPECIES:
        contrast = CONTRASTS[species]
        for metric, df in tables.items():
            high, low = two_group_arrays(df, species, contrast["high"], contrast["low"])
            diff = mean_diff(high, low)
            ci_low, ci_high = bootstrap_ci(high, low, seed=20260530 + sum(ord(c) for c in f"{species}_{metric}_boot"))
            effect_d = cohens_d(high, low)
            d_ci_low, d_ci_high = bootstrap_ci(
                high,
                low,
                statistic=cohens_d,
                seed=20260530 + sum(ord(c) for c in f"{species}_{metric}_d_boot"),
            )
            p_greater, p_two = permutation_p(high, low, seed=20260530 + sum(ord(c) for c in f"{species}_{metric}_perm"))
            rows.append(
                {
                    "species": species,
                    "metric": metric,
                    "contrast": contrast["label"],
                    "interpretation": contrast["interpretation"],
                    "high_orders": ",".join(map(str, contrast["high"])),
                    "low_orders": ",".join(map(str, contrast["low"])),
                    "n_high": int(len(high)),
                    "n_low": int(len(low)),
                    "mean_high": float(np.nanmean(high)) if len(high) else np.nan,
                    "mean_low": float(np.nanmean(low)) if len(low) else np.nan,
                    "mean_diff": diff,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "cohens_d": effect_d,
                    "cohens_d_ci_low": d_ci_low,
                    "cohens_d_ci_high": d_ci_high,
                    "permutation_greater_p": p_greater,
                    "permutation_two_sided_p": p_two,
                    "n_perm": 5000,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_EFFECTS, index=False)
    return out


def group_signature(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for metric, df in tables.items():
        for (species, order), sub in df.groupby(["species", "shared_fine_order"], sort=False):
            species_df = df[df["species"].eq(species)].copy()
            species_values = species_df["value"].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
            group_values = sub["value"].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
            other_values = (
                species_df.loc[~species_df["shared_fine_order"].eq(order), "value"]
                .replace([np.inf, -np.inf], np.nan)
                .dropna()
                .to_numpy(float)
            )
            mean = float(np.nanmean(group_values)) if len(group_values) else np.nan
            other_mean = float(np.nanmean(other_values)) if len(other_values) else np.nan
            spread = float(np.nanstd(species_values, ddof=1)) if len(species_values) > 1 else np.nan
            relative_z_mean = (mean - np.nanmean(species_values)) / spread if np.isfinite(spread) and spread > 0 else np.nan
            leave_one_group_out_z_delta = (mean - other_mean) / spread if np.isfinite(spread) and spread > 0 else np.nan
            rows.append(
                {
                    "species": species,
                    "shared_fine_order": int(order),
                    "shared_fine_label": GROUP_LABELS.get(int(order), str(order)),
                    "metric": metric,
                    "mean": mean,
                    "other_group_mean": other_mean,
                    "mean_minus_other_group_mean": mean - other_mean if np.isfinite(mean) and np.isfinite(other_mean) else np.nan,
                    "relative_z_mean": relative_z_mean,
                    "leave_one_group_out_z_delta": leave_one_group_out_z_delta,
                    "sem": float(pd.Series(group_values).sem()),
                    "n": int(len(group_values)),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_GROUP, index=False)
    return out


def _zdelta_for_group(values: np.ndarray, groups: np.ndarray, group_order: int) -> float:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=int)
    finite = np.isfinite(values)
    values = values[finite]
    groups = groups[finite]
    in_group = groups == int(group_order)
    if in_group.sum() == 0 or (~in_group).sum() == 0:
        return np.nan
    spread = float(np.nanstd(values, ddof=1))
    if not np.isfinite(spread) or spread <= 0:
        return np.nan
    return float((np.nanmean(values[in_group]) - np.nanmean(values[~in_group])) / spread)


def group_permutation_tests(
    tables: dict[str, pd.DataFrame],
    signature: pd.DataFrame,
    n_perm: int = 1000,
) -> pd.DataFrame:
    rows = []
    for metric, df in tables.items():
        for species in SPECIES:
            sub = (
                df[df["species"].eq(species)]
                .replace([np.inf, -np.inf], np.nan)
                .dropna(subset=["value", "shared_fine_order"])
                .copy()
            )
            if sub.empty:
                continue
            values = sub["value"].to_numpy(float)
            groups = sub["shared_fine_order"].to_numpy(int)
            unique_groups = sorted(np.unique(groups).astype(int).tolist())
            observed = {
                group: _zdelta_for_group(values, groups, group)
                for group in unique_groups
            }
            null = {group: np.empty(n_perm, dtype=float) for group in unique_groups}
            rng = np.random.default_rng(
                20260601 + sum(ord(ch) for ch in f"{species}_{metric}_group_zdelta")
            )
            for perm_idx in range(n_perm):
                perm_groups = rng.permutation(groups)
                for group in unique_groups:
                    null[group][perm_idx] = _zdelta_for_group(values, perm_groups, group)

            p_values = []
            for group in unique_groups:
                obs = observed[group]
                if not np.isfinite(obs):
                    p = np.nan
                else:
                    perm_vals = null[group][np.isfinite(null[group])]
                    # Magnitude-based permutation test:
                    # treat strong enrichment or depletion as large |ZDelta|
                    # and compare against the absolute permutation null.
                    p = float((np.sum(np.abs(perm_vals) >= abs(obs)) + 1) / (len(perm_vals) + 1))
                p_values.append(p)

            finite_idx = [idx for idx, p in enumerate(p_values) if np.isfinite(p)]
            p_holm = np.full(len(p_values), np.nan, dtype=float)
            q_fdr = np.full(len(p_values), np.nan, dtype=float)
            if finite_idx:
                finite_p = np.asarray([p_values[idx] for idx in finite_idx], dtype=float)
                order = np.argsort(finite_p)
                adjusted = np.empty_like(finite_p)
                running = 0.0
                m = len(finite_p)
                for rank, sorted_idx in enumerate(order):
                    value = min(1.0, (m - rank) * finite_p[sorted_idx])
                    running = max(running, value)
                    adjusted[sorted_idx] = running
                for local_idx, global_idx in enumerate(finite_idx):
                    p_holm[global_idx] = adjusted[local_idx]

                # Benjamini-Hochberg FDR correction within each species x metric.
                bh_adjusted = np.empty_like(finite_p)
                running_q = 1.0
                for rank_from_end, sorted_idx in enumerate(order[::-1]):
                    rank = len(finite_p) - rank_from_end
                    value = finite_p[sorted_idx] * len(finite_p) / rank
                    running_q = min(running_q, value)
                    bh_adjusted[sorted_idx] = min(1.0, running_q)
                for local_idx, global_idx in enumerate(finite_idx):
                    q_fdr[global_idx] = bh_adjusted[local_idx]

            for idx, group in enumerate(unique_groups):
                label = GROUP_LABELS.get(group, str(group))
                rows.append(
                    {
                        "species": species,
                        "metric": metric,
                        "shared_fine_order": group,
                        "shared_fine_label": label,
                        "observed_zdelta": observed[group],
                        "p_abs_permutation": p_values[idx],
                        "p_holm_within_species_metric": p_holm[idx],
                        "q_fdr_within_species_metric": q_fdr[idx],
                        "significant_holm_0.05": bool(np.isfinite(p_holm[idx]) and p_holm[idx] < 0.05),
                        "significant_fdr_0.05": bool(np.isfinite(q_fdr[idx]) and q_fdr[idx] < 0.05),
                        "n_group": int(np.sum(groups == group)),
                        "n_total": int(len(groups)),
                        "n_perm": int(n_perm),
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_GROUP_PERM, index=False)
    return out


def _valid_cache(path: Path, required_columns: set[str]) -> bool:
    if not path.exists():
        return False
    try:
        cols = set(pd.read_csv(path, nrows=0).columns)
    except Exception:
        return False
    return required_columns.issubset(cols)


def load_or_compute_statistics(
    tables: dict[str, pd.DataFrame] | None = None,
    *,
    recompute: bool = False,
    effects_path: Path = OUT_EFFECTS,
    group_path: Path = OUT_GROUP,
    group_perm_path: Path = OUT_GROUP_PERM,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load cached Figure 7 statistics unless explicit recomputation is needed."""
    cache_ok = (
        _valid_cache(effects_path, STAT_CACHE_COLUMNS["effects"])
        and _valid_cache(group_path, STAT_CACHE_COLUMNS["signature"])
        and _valid_cache(group_perm_path, STAT_CACHE_COLUMNS["group_perm"])
    )
    if not recompute and cache_ok:
        if tables is None:
            tables = metric_tables()
        return (
            tables,
            pd.read_csv(effects_path),
            pd.read_csv(group_path),
            pd.read_csv(group_perm_path),
        )

    if tables is None:
        tables = metric_tables()
    effects = compute_effects(tables)
    signature = group_signature(tables)
    group_perm = group_permutation_tests(tables, signature)
    effects.to_csv(effects_path, index=False)
    signature.to_csv(group_path, index=False)
    group_perm.to_csv(group_perm_path, index=False)
    return tables, effects, signature, group_perm


def draw_concept(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    boxes = [
        (0.06, 0.62, 0.23, 0.20, "Functional\ngroup identity", "#F4F4F4"),
        (0.39, 0.62, 0.23, 0.20, "Directed SC\npolarity", "#F4F4F4"),
        (0.72, 0.62, 0.23, 0.20, "Dynamic FC\nvariability", "#F4F4F4"),
    ]
    for x, y, w, h, text, fc in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.012", fc=fc, ec="#333333", lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8, fontweight="bold")
    for x0, x1 in [(0.30, 0.39), (0.63, 0.72)]:
        ax.add_patch(FancyArrowPatch((x0, 0.72), (x1, 0.72), arrowstyle="-|>", mutation_scale=10, lw=1.0, color="#333333"))

    y0 = 0.30
    species_rows = [
        ("C. elegans", "compact sensory /\nintegrative routing", "weak DCA\npolarity"),
        ("Drosophila", "MB / olfactory\ninput-side axis", "Pre-DCA"),
        ("Zebrafish", "telencephalic\noutput cascade", "Post-DCA"),
    ]
    for i, (species, left, right) in enumerate(species_rows):
        y = y0 - i * 0.12
        ax.add_patch(Circle((0.08, y), 0.026, color=SPECIES_COLORS[species], ec="#222222", lw=0.4))
        ax.text(0.13, y, species, ha="left", va="center", fontsize=7.4, fontweight="bold", color=SPECIES_COLORS[species])
        ax.text(0.37, y, left, ha="center", va="center", fontsize=6.7)
        ax.add_patch(FancyArrowPatch((0.52, y), (0.62, y), arrowstyle="-|>", mutation_scale=8, lw=0.8, color="#555555"))
        ax.text(0.76, y, right, ha="center", va="center", fontsize=6.9, fontweight="bold")
    ax.text(0.02, 0.94, "A", fontsize=13, fontweight="bold")
    ax.text(0.06, 0.90, "Structure-function synthesis model", fontsize=10, fontweight="bold")


def draw_effect_forest(ax: plt.Axes, effects: pd.DataFrame) -> None:
    ax.axvline(0, color="#777777", lw=0.8, ls="--")
    y_positions = []
    labels = []
    y = 0
    for species in SPECIES:
        for metric in ["FCV", "Post-DCA", "Pre-DCA"]:
            row = effects[(effects["species"].eq(species)) & (effects["metric"].eq(metric))].iloc[0]
            y_positions.append(y)
            labels.append(f"{species}  {metric}")
            color = METRIC_COLORS[metric]
            ax.plot([row["cohens_d_ci_low"], row["cohens_d_ci_high"]], [y, y], color=color, lw=1.3, solid_capstyle="round")
            ax.scatter(row["cohens_d"], y, s=28, color=SPECIES_COLORS[species], edgecolors="#222222", linewidths=0.35, zorder=3)
            p = row["permutation_two_sided_p"]
            x_text = max(row["cohens_d_ci_low"], row["cohens_d_ci_high"]) + 0.16
            ax.text(x_text, y, f"p={format_p(p)}", va="center", fontsize=5.8)
            y += 1
        y += 0.55
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlabel("planned contrast standardized effect size\nCohen's d, high-FCV group minus comparison group")
    ax.set_title("B  Core planned contrasts", loc="left", fontsize=10, fontweight="bold", pad=4)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def significance_lookup(effects: pd.DataFrame, alpha: float = 0.05) -> dict[tuple[str, int, str], float]:
    lookup: dict[tuple[str, int, str], float] = {}
    for _, row in effects.iterrows():
        p = float(row["p_holm_within_species_metric"])
        if not np.isfinite(p) or p >= alpha:
            continue
        lookup[(row["species"], int(row["shared_fine_order"]), row["metric"])] = p
    return lookup


def raw_p_lookup(effects: pd.DataFrame) -> dict[tuple[str, int, str], float]:
    lookup: dict[tuple[str, int, str], float] = {}
    for _, row in effects.iterrows():
        p = float(row["p_abs_permutation"])
        if np.isfinite(p):
            lookup[(row["species"], int(row["shared_fine_order"]), row["metric"])] = p
    return lookup


def p_to_stars(p: float) -> str:
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    return "*"


def draw_signature_matrix(
    ax: plt.Axes,
    signature: pd.DataFrame,
    group_perm: pd.DataFrame,
) -> None:
    ax.set_title("Functional-group signature", loc="left", fontsize=12, fontweight="bold", pad=6)
    groups_by_species = {
        sp: [
            group
            for group in DISPLAY_GROUP_ORDER
            if group in set(signature[signature["species"].eq(sp)]["shared_fine_order"].unique().astype(int))
        ]
        for sp in SPECIES
    }
    row_step = 4.35
    y_base = 0.0
    yticks = []
    ylabels = []
    separators = []
    max_cols = max(len(v) for v in groups_by_species.values())
    metrics = ["FCV", "Post-DCA", "Pre-DCA"]
    sig = significance_lookup(group_perm)
    fdr_q = fdr_q_lookup(group_perm)
    for sp in SPECIES:
        orders = groups_by_species[sp]
        for ref_value in [-2, -1, 0, 1, 2]:
            y_ref = y_base + ref_value
            if ref_value == 0:
                continue
            ax.axhline(y_ref, color="#EEEEEE", lw=0.38, zorder=-2)
        for col, order in enumerate(orders):
            if order in HIGHLIGHT_GROUPS:
                ax.add_patch(
                    Rectangle(
                        (col - 0.52, y_base - 1.45),
                        1.00,
                        2.90,
                        fc="none",
                        ec=GROUP_COLORS.get(order, "#AAAAAA"),
                        lw=1.7,
                        alpha=0.95,
                        zorder=0,
                    )
                )
                ax.add_patch(
                    Rectangle(
                        (col - 0.52, y_base - 1.45),
                        1.00,
                        2.90,
                        fc=GROUP_COLORS.get(order, "#EEEEEE"),
                        ec="none",
                        alpha=0.055,
                        zorder=-1,
                    )
                )
            vals = []
            for metric in metrics:
                row = signature[
                    signature["species"].eq(sp)
                    & signature["shared_fine_order"].eq(order)
                    & signature["metric"].eq(metric)
                ]
                vals.append(float(row["leave_one_group_out_z_delta"].iloc[0]) if not row.empty else np.nan)
            for k, val in enumerate(vals):
                x = col + (k - 1) * 0.22
                height = np.clip(val, -2.0, 2.0)
                color = METRIC_COLORS[metrics[k]]
                ax.bar(x, height, width=0.18, bottom=y_base, color=color, alpha=0.85, zorder=2)
                p = sig.get((sp, order, metrics[k]))
                if p is not None and np.isfinite(height):
                    y_star = y_base + height + (0.22 if height >= 0 else -0.30)
                    ax.text(
                        x,
                        y_star,
                        p_to_stars(p),
                        ha="center",
                        va="bottom" if height >= 0 else "top",
                        fontsize=8.0,
                        fontweight="bold",
                        color="#111111",
                        zorder=4,
                    )
                elif (
                    (sp == "Zebrafish" and order == 0 and metrics[k] == "FCV")
                    or (sp == "Drosophila" and order == 0 and metrics[k] in {"Post-DCA", "Pre-DCA"})
                ):
                    q_value = fdr_q.get((sp, order, metrics[k]))
                    if q_value is not None and np.isfinite(height):
                        y_text = y_base + height + (0.20 if height >= 0 else -0.28)
                        ax.text(
                            x,
                            y_text,
                            f"q={q_value:.3f}",
                            ha="center",
                            va="bottom" if height >= 0 else "top",
                            fontsize=5.8,
                            fontweight="bold",
                            color="#111111",
                            zorder=4,
                        )
            label_weight = "bold" if order in HIGHLIGHT_GROUPS else "normal"
            label_size = 7.8 if order in HIGHLIGHT_GROUPS else 7.0
            ax.text(
                col,
                y_base - 0.36,
                GROUP_LABELS.get(order, str(order)),
                ha="center",
                va="top",
                fontsize=label_size,
                fontweight=label_weight,
                rotation=35,
            )
            ax.add_patch(
                Rectangle(
                    (col - 0.42, y_base - 0.025),
                    0.84,
                    0.05,
                    fc=GROUP_COLORS.get(order, "#BDBDBD"),
                    ec="#333333" if order in HIGHLIGHT_GROUPS else "none",
                    lw=0.5 if order in HIGHLIGHT_GROUPS else 0.0,
                    alpha=0.95,
                    zorder=1,
                )
            )
        yticks.append(y_base)
        ylabels.append(sp)
        separators.append(y_base + row_step / 2.0)
        y_base += row_step
    for sep in separators[:-1]:
        ax.axhline(sep, color="#DDDDDD", lw=0.55)
    ax.set_xlim(-0.7, max_cols - 0.3)
    ax.set_ylim(-0.85, y_base - row_step + 2.05)
    numeric_ticks = []
    numeric_labels = []
    for base in yticks:
        for ref_value in [-2, -1, 0, 1, 2]:
            y_ref = base + ref_value
            if ax.get_ylim()[0] <= y_ref <= ax.get_ylim()[1]:
                numeric_ticks.append(y_ref)
                numeric_labels.append(f"{ref_value:g}")
    ax.set_yticks(numeric_ticks)
    ax.set_yticklabels(numeric_labels, fontsize=5.8, color="#666666")
    ax.tick_params(axis="y", length=2.0, width=0.45, colors="#666666", pad=2)
    for base, label in zip(yticks, ylabels, strict=False):
        ax.text(
            -0.96,
            base,
            label,
            ha="right",
            va="center",
            fontsize=9.0,
            fontweight="bold",
            color="#111111",
            clip_on=False,
        )
    ax.set_xticks([])
    ax.set_ylabel("ZDelta\n(group mean - other groups, within-species SD)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def draw_signature_plot_panels(
    axes: np.ndarray,
    signature: pd.DataFrame,
    group_perm: pd.DataFrame,
) -> None:
    metrics = ["FCV", "Post-DCA", "Pre-DCA"]
    offsets = {"FCV": -0.22, "Post-DCA": 0.0, "Pre-DCA": 0.22}
    width = 0.18
    groups_by_species = {
        sp: [
            group
            for group in DISPLAY_GROUP_ORDER
            if group in set(signature[signature["species"].eq(sp)]["shared_fine_order"].unique().astype(int))
        ]
        for sp in SPECIES
    }
    sig = significance_lookup(group_perm)
    raw_p = raw_p_lookup(group_perm)
    for panel_idx, (ax, sp) in enumerate(zip(axes, SPECIES, strict=False)):
        orders = groups_by_species[sp]
        x_positions = np.arange(len(orders), dtype=float)
        sp_values = (
            signature.loc[signature["species"].eq(sp), "leave_one_group_out_z_delta"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
        )
        if len(sp_values):
            y_abs = float(np.nanmax(np.abs(sp_values)))
            y_lim = float(np.ceil((y_abs + 0.28) * 2) / 2)
            y_lim = max(1.25, min(3.7, y_lim))
        else:
            y_lim = 2.0
        y_min, y_max = -y_lim, y_lim
        for idx, order in enumerate(orders):
            if order in HIGHLIGHT_GROUPS:
                ax.axvspan(
                    idx - 0.48,
                    idx + 0.48,
                    color=GROUP_COLORS.get(order, "#EEEEEE"),
                    alpha=0.055,
                    zorder=-3,
                )
                ax.add_patch(
                    Rectangle(
                        (idx - 0.48, y_min),
                        0.96,
                        y_max - y_min,
                        fc="none",
                        ec=GROUP_COLORS.get(order, "#AAAAAA"),
                        lw=0.95,
                        alpha=0.90,
                        zorder=4,
                        clip_on=False,
                    )
                )
            ax.add_patch(
                Rectangle(
                    (idx - 0.40, -0.028),
                    0.80,
                    0.056,
                    fc=GROUP_COLORS.get(order, "#BDBDBD"),
                    ec="none",
                    alpha=0.78,
                    zorder=1,
                )
            )
            for metric in metrics:
                row = signature[
                    signature["species"].eq(sp)
                    & signature["shared_fine_order"].eq(order)
                    & signature["metric"].eq(metric)
                ]
                if row.empty:
                    continue
                val = float(row["leave_one_group_out_z_delta"].iloc[0])
                if not np.isfinite(val):
                    continue
                x = idx + offsets[metric]
                ax.bar(
                    x,
                    val,
                    width=width,
                    color=METRIC_COLORS[metric],
                    alpha=0.96,
                    zorder=3,
                )
                q = sig.get((sp, order, metric))
                if q is not None:
                    ax.text(
                        x,
                        val + (0.13 if val >= 0 else -0.16),
                        p_to_stars(q),
                        ha="center",
                        va="bottom" if val >= 0 else "top",
                        fontsize=7.0,
                        fontweight="bold",
                        color="#111111",
                        clip_on=False,
                        zorder=5,
                    )
                elif (
                    (sp == "Zebrafish" and order == 0 and metric == "FCV")
                    or (sp == "Drosophila" and order == 0 and metric in {"Post-DCA", "Pre-DCA"})
                ):
                    p_value = raw_p.get((sp, order, metric))
                    if p_value is not None:
                        ax.text(
                            x,
                            val + (0.13 if val >= 0 else -0.16),
                            f"p={p_value:.3f}",
                            ha="center",
                            va="bottom" if val >= 0 else "top",
                            fontsize=5.2,
                            fontweight="bold",
                            color="#111111",
                            clip_on=False,
                            zorder=5,
                        )
        ax.axhline(0, color="#2F2F2F", lw=0.55, zorder=2)
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(-0.65, len(orders) - 0.35)
        ax.set_ylabel("ZDelta", fontsize=6.4)
        ax.text(
            0.50,
            1.03,
            sp,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=7.5,
            fontweight="bold",
            color="#111111",
            clip_on=False,
        )
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 1.5, 2, 2.5, 3, 4]))
        ax.grid(axis="y", color="#EDEDED", lw=0.35, zorder=-4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_color("#BFBFBF")
        ax.spines["bottom"].set_linewidth(0.45)
        ax.spines["left"].set_linewidth(0.45)
        ax.tick_params(axis="y", labelsize=5.8, width=0.45, length=2.2, colors="#333333")
        ax.tick_params(axis="x", length=0)
        ax.set_xticks(x_positions)
        if panel_idx == len(axes) - 1:
            ax.set_xticklabels(
                [GROUP_LABELS.get(order, str(order)) for order in orders],
                rotation=28,
                ha="right",
                fontsize=6.2,
            )
            for tick_label, order in zip(ax.get_xticklabels(), orders, strict=False):
                tick_label.set_fontweight("bold" if order in HIGHLIGHT_GROUPS else "normal")
        else:
            ax.set_xticklabels([])


def draw_boxstrip_grid(
    axes: np.ndarray,
    tables: dict[str, pd.DataFrame],
    group_perm: pd.DataFrame,
) -> None:
    metrics = ["FCV", "Post-DCA", "Pre-DCA"]
    y_labels = {"FCV": "within-recording z", "Post-DCA": "SC measure value", "Pre-DCA": "SC measure value"}
    sig = significance_lookup(group_perm)
    raw_p = raw_p_lookup(group_perm)

    for row_idx, sp in enumerate(SPECIES):
        for col_idx, metric in enumerate(metrics):
            ax = axes[row_idx, col_idx]
            data = (
                tables[metric]
                .loc[tables[metric]["species"].eq(sp)]
                .replace([np.inf, -np.inf], np.nan)
                .dropna(subset=["value", "shared_fine_order"])
                .copy()
            )
            orders = [
                group
                for group in DISPLAY_GROUP_ORDER
                if group in set(data["shared_fine_order"].astype(int).unique())
            ]
            values_by_group: list[np.ndarray] = []
            rng = np.random.default_rng(20260601 + sum(ord(ch) for ch in f"{sp}-{metric}-strip"))
            for idx, order in enumerate(orders):
                vals = (
                    data.loc[data["shared_fine_order"].astype(int).eq(order), "value"]
                    .dropna()
                    .to_numpy(float)
                )
                values_by_group.append(vals)
                if len(vals):
                    ax.scatter(
                        np.full(len(vals), idx) + rng.normal(0, 0.060, size=len(vals)),
                        vals,
                        s=9.0,
                        marker="o",
                        facecolor=GROUP_COLORS.get(order, "#BBBBBB"),
                        edgecolor="#222222",
                        alpha=0.55,
                        linewidth=0.22,
                        rasterized=True,
                        zorder=3,
                    )
            box = ax.boxplot(
                values_by_group,
                positions=np.arange(len(orders)),
                widths=0.48,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "#222222", "linewidth": 0.75},
                boxprops={"linewidth": 0.65, "edgecolor": "#222222"},
                whiskerprops={"linewidth": 0.55, "color": "#222222"},
                capprops={"linewidth": 0.55, "color": "#222222"},
            )
            for patch, order in zip(box["boxes"], orders, strict=False):
                color = GROUP_COLORS.get(order, "#BBBBBB")
                patch.set_facecolor(color)
                patch.set_alpha(0.28)
                patch.set_edgecolor(color)
                patch.set_linewidth(0.90)
                patch.set_zorder(1)

            finite_vals = data["value"].dropna().to_numpy(float)
            if len(finite_vals):
                ymin, ymax = float(np.nanmin(finite_vals)), float(np.nanmax(finite_vals))
                yr = max(ymax - ymin, 1e-6)
                ax.set_ylim(ymin - 0.15 * yr, ymax + 0.22 * yr)
            y0, y1 = ax.get_ylim()
            yr = y1 - y0
            for idx, order in enumerate(orders):
                vals = values_by_group[idx]
                if len(vals):
                    y_text = float(np.nanmax(vals)) + 0.045 * yr
                else:
                    y_text = y1 - 0.08 * yr
                q = sig.get((sp, int(order), metric))
                if q is not None:
                    ax.text(
                        idx,
                        y_text,
                        p_to_stars(q),
                        ha="center",
                        va="bottom",
                        fontsize=6.0,
                        fontweight="bold",
                        clip_on=False,
                        zorder=5,
                    )
                elif (
                    (sp == "Zebrafish" and order == 0 and metric == "FCV")
                    or (sp == "Drosophila" and order == 0 and metric in {"Post-DCA", "Pre-DCA"})
                ):
                    p_value = raw_p.get((sp, int(order), metric))
                    if p_value is not None:
                        ax.text(
                            idx,
                            y_text,
                            f"p={p_value:.3f}",
                            ha="center",
                            va="bottom",
                            fontsize=4.9,
                            fontweight="bold",
                            clip_on=False,
                            zorder=5,
                        )
            ax.axhline(0, color="#777777", lw=0.50, ls=":", zorder=0)
            if row_idx == 0:
                ax.set_title(metric, fontsize=7.8, pad=3)
            if col_idx == 0:
                ax.text(
                    -0.42,
                    0.5,
                    sp,
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=7.2,
                    fontweight="bold",
                    clip_on=False,
                )
            ax.set_ylabel(y_labels[metric] if col_idx == 0 else "", fontsize=6.2)
            ax.set_xticks(np.arange(len(orders)))
            if row_idx == len(SPECIES) - 1:
                ax.set_xticklabels(
                    [GROUP_LABELS.get(int(order), str(order)) for order in orders],
                    rotation=38,
                    ha="right",
                    fontsize=5.2,
                )
                for tick, order in zip(ax.get_xticklabels(), orders, strict=False):
                    tick.set_color(GROUP_COLORS.get(int(order), "#555555"))
                    tick.set_fontweight("bold")
            else:
                ax.set_xticklabels([])
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.grid(axis="y", color="#EDEDED", lw=0.33, zorder=-4)
            ax.tick_params(axis="both", direction="out", pad=1.2, length=2.0, width=0.45, labelsize=5.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_color("#BFBFBF")
            ax.spines["bottom"].set_linewidth(0.45)
            ax.spines["left"].set_linewidth(0.45)


def draw_combined_metric_boxstrip_panels(
    axes: np.ndarray,
    tables: dict[str, pd.DataFrame],
    group_perm: pd.DataFrame,
) -> None:
    """Draw one Figure-14C-like distribution panel per species.

    FCV and DCA measures live on different raw scales, so each metric is
    z-scored within species before plotting. This preserves the group-level
    ordering while allowing FCV, Post-DCA, and Pre-DCA to be read together.
    """
    metrics = ["FCV", "Post-DCA", "Pre-DCA"]
    offsets = {"FCV": -0.24, "Post-DCA": 0.0, "Pre-DCA": 0.24}
    metric_labels = {"FCV": "FCV", "Post-DCA": "Post-DCA", "Pre-DCA": "Pre-DCA"}
    metric_markers = {"FCV": "o", "Post-DCA": "^", "Pre-DCA": "s"}
    sig = significance_lookup(group_perm)
    raw_p = raw_p_lookup(group_perm)

    for panel_idx, (ax, sp) in enumerate(zip(axes, SPECIES, strict=False)):
        sp_tables = {}
        available_groups: set[int] = set()
        for metric in metrics:
            data = (
                tables[metric]
                .loc[tables[metric]["species"].eq(sp)]
                .replace([np.inf, -np.inf], np.nan)
                .dropna(subset=["value", "shared_fine_order"])
                .copy()
            )
            vals = data["value"].to_numpy(float)
            mean = float(np.nanmean(vals)) if len(vals) else np.nan
            sd = float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else np.nan
            if np.isfinite(sd) and sd > 0:
                data["plot_value"] = (data["value"] - mean) / sd
            else:
                data["plot_value"] = np.nan
            sp_tables[metric] = data
            available_groups.update(data["shared_fine_order"].astype(int).unique().tolist())

        orders = [group for group in DISPLAY_GROUP_ORDER if group in available_groups]
        x_positions = np.arange(len(orders), dtype=float)
        rng = np.random.default_rng(20260601 + sum(ord(ch) for ch in f"{sp}-combined-strip"))

        all_plot_values = []
        for idx, order in enumerate(orders):
            if order in HIGHLIGHT_GROUPS:
                ax.axvspan(
                    idx - 0.50,
                    idx + 0.50,
                    color=GROUP_COLORS.get(order, "#EEEEEE"),
                    alpha=0.050,
                    zorder=-4,
                )
                ax.add_patch(
                    Rectangle(
                        (idx - 0.50, -3.0),
                        1.0,
                        6.0,
                        fc="none",
                        ec=GROUP_COLORS.get(order, "#AAAAAA"),
                        lw=0.95,
                        alpha=0.90,
                        clip_on=False,
                        zorder=5,
                    )
                )

            for metric in metrics:
                data = sp_tables[metric]
                vals = (
                    data.loc[data["shared_fine_order"].astype(int).eq(order), "plot_value"]
                    .dropna()
                    .to_numpy(float)
                )
                if len(vals) == 0:
                    continue
                all_plot_values.extend(vals.tolist())
                x = idx + offsets[metric]
                group_color = GROUP_COLORS.get(order, "#BBBBBB")
                box = ax.boxplot(
                    [vals],
                    positions=[x],
                    widths=0.18,
                    patch_artist=True,
                    showfliers=False,
                    medianprops={"color": "#222222", "linewidth": 0.70},
                    boxprops={"linewidth": 0.70, "edgecolor": group_color},
                    whiskerprops={"linewidth": 0.55, "color": group_color},
                    capprops={"linewidth": 0.55, "color": group_color},
                )
                box["boxes"][0].set_facecolor(group_color)
                box["boxes"][0].set_alpha(0.16)
                box["boxes"][0].set_zorder(1)
                ax.scatter(
                    np.full(len(vals), x) + rng.normal(0, 0.024, size=len(vals)),
                    vals,
                    s=5.8,
                    marker=metric_markers[metric],
                    facecolor=group_color,
                    edgecolor="#222222",
                    alpha=0.36,
                    linewidth=0.16,
                    rasterized=True,
                    zorder=3,
                )

                q = sig.get((sp, int(order), metric))
                p_value = raw_p.get((sp, int(order), metric))
                if q is not None or (
                    (sp == "Zebrafish" and order == 0 and metric == "FCV")
                    or (sp == "Drosophila" and order == 0 and metric in {"Post-DCA", "Pre-DCA"})
                ):
                    y_text = float(np.nanmax(vals)) + 0.18
                    label = p_to_stars(q) if q is not None else f"p={p_value:.3f}"
                    ax.text(
                        x,
                        y_text,
                        label,
                        ha="center",
                        va="bottom",
                        fontsize=5.2 if q is None else 6.3,
                        fontweight="bold",
                        clip_on=False,
                        zorder=6,
                    )

        if all_plot_values:
            vals_for_scale = np.asarray(all_plot_values, dtype=float)
            ymin = float(np.nanpercentile(vals_for_scale, 1.0))
            ymax = float(np.nanpercentile(vals_for_scale, 99.0))
            yr = max(ymax - ymin, 1e-6)
            ymin = max(float(np.nanmin(vals_for_scale)) - 0.02 * yr, ymin - 0.12 * yr)
            ymax = min(float(np.nanmax(vals_for_scale)) + 0.02 * yr, ymax + 0.24 * yr)
            ax.set_ylim(ymin, ymax)
        ax.axhline(0, color="#777777", lw=0.50, ls=":", zorder=0)
        ax.set_xlim(-0.65, len(orders) - 0.35)
        ax.set_title(sp, fontsize=8.4, fontweight="bold", pad=3)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([GROUP_LABELS.get(int(order), str(order)) for order in orders], rotation=0, ha="center", fontsize=5.8)
        for tick, order in zip(ax.get_xticklabels(), orders, strict=False):
            tick.set_color(GROUP_COLORS.get(int(order), "#070707"))
            tick.set_fontweight("bold" if order in HIGHLIGHT_GROUPS else "normal")
        ax.set_ylabel("within-species metric z", fontsize=6.3)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.grid(axis="y", color="#EDEDED", lw=0.33, zorder=-5)
        ax.tick_params(axis="both", direction="out", pad=1.2, length=2.0, width=0.45, labelsize=5.6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_color("#BFBFBF")
        ax.spines["bottom"].set_linewidth(0.45)
        ax.spines["left"].set_linewidth(0.45)

        if panel_idx == 0:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker=metric_markers[metric],
                    color="none",
                    markerfacecolor="#777777",
                    markeredgecolor="#222222",
                    markersize=4.4,
                    label=metric_labels[metric],
                )
                for metric in metrics
            ]
            ax.legend(
                handles=handles,
                loc="upper right",
                bbox_to_anchor=(1.0, 1.03),
                frameon=False,
                ncol=3,
                columnspacing=0.9,
                handletextpad=0.3,
                borderaxespad=0.0,
            )


def draw_measure_comparison_panels(
    axes: np.ndarray,
    tables: dict[str, pd.DataFrame],
    group_perm: pd.DataFrame,
) -> None:
    """Draw one panel per measure so species can be compared directly."""
    metrics = ["FCV", "Post-DCA", "Pre-DCA"]
    species_offsets = {"C. elegans": -0.30, "Drosophila": 0.0, "Zebrafish": 0.30}
    species_markers = {"C. elegans": "o", "Drosophila": "^", "Zebrafish": "s"}
    sig = significance_lookup(group_perm)
    raw_p = raw_p_lookup(group_perm)

    all_groups = []
    for metric in metrics:
        for sp in SPECIES:
            data = tables[metric].loc[tables[metric]["species"].eq(sp)]
            all_groups.extend(data["shared_fine_order"].dropna().astype(int).unique().tolist())
    orders = [group for group in DISPLAY_GROUP_ORDER if group in set(all_groups)]
    x_positions = np.arange(len(orders), dtype=float)

    for panel_idx, (ax, metric) in enumerate(zip(axes, metrics, strict=False)):
        all_plot_values = []
        for idx, order in enumerate(orders):
            if order in HIGHLIGHT_GROUPS:
                ax.axvspan(
                    idx - 0.50,
                    idx + 0.50,
                    color=GROUP_COLORS.get(order, "#EEEEEE"),
                    alpha=0.052,
                    zorder=-4,
                )
                ax.add_patch(
                    Rectangle(
                        (idx - 0.50, -3.0),
                        1.0,
                        6.0,
                        fc="none",
                        ec=GROUP_COLORS.get(order, "#AAAAAA"),
                        lw=1.05,
                        alpha=0.82,
                        clip_on=False,
                        zorder=5,
                    )
                )

        for sp in SPECIES:
            data = (
                tables[metric]
                .loc[tables[metric]["species"].eq(sp)]
                .replace([np.inf, -np.inf], np.nan)
                .dropna(subset=["value", "shared_fine_order"])
                .copy()
            )
            vals_all = data["value"].to_numpy(float)
            mean = float(np.nanmean(vals_all)) if len(vals_all) else np.nan
            sd = float(np.nanstd(vals_all, ddof=1)) if len(vals_all) > 1 else np.nan
            if np.isfinite(sd) and sd > 0:
                data["plot_value"] = (data["value"] - mean) / sd
            else:
                data["plot_value"] = np.nan
            rng = np.random.default_rng(20260601 + sum(ord(ch) for ch in f"{metric}-{sp}-measure-strip"))

            for idx, order in enumerate(orders):
                vals = (
                    data.loc[data["shared_fine_order"].astype(int).eq(order), "plot_value"]
                    .dropna()
                    .to_numpy(float)
                )
                if len(vals) == 0:
                    continue
                all_plot_values.extend(vals.tolist())
                x = idx + species_offsets[sp]
                species_color = SPECIES_COLORS.get(sp, "#777777")
                ax.scatter(
                    np.full(len(vals), x) + rng.normal(0, 0.026, size=len(vals)),
                    vals,
                    s=6.2,
                    marker=species_markers[sp],
                    facecolor=species_color,
                    edgecolor="none",
                    alpha=0.26,
                    linewidth=0.0,
                    rasterized=True,
                    zorder=2,
                )
                center, ci_low, ci_high = mean_ci95(vals)
                if np.isfinite(center):
                    ax.plot(
                        [x, x],
                        [ci_low, ci_high],
                        color=species_color,
                        lw=1.25,
                        alpha=0.95,
                        solid_capstyle="round",
                        zorder=4,
                    )
                    ax.scatter(
                        [x],
                        [center],
                        s=34,
                        marker=species_markers[sp],
                        facecolor=species_color,
                        edgecolor="#111111",
                        alpha=0.98,
                        linewidth=0.55,
                        zorder=5,
                    )
                    ax.plot(
                        [x - 0.075, x + 0.075],
                        [ci_low, ci_low],
                        color=species_color,
                        lw=0.85,
                        alpha=0.95,
                        zorder=4,
                    )
                    ax.plot(
                        [x - 0.075, x + 0.075],
                        [ci_high, ci_high],
                        color=species_color,
                        lw=0.85,
                        alpha=0.95,
                        zorder=4,
                    )

                    ax.scatter(
                        [x],
                        [center],
                        s=70,
                        marker=species_markers[sp],
                        facecolor="none",
                        edgecolor=species_color,
                        alpha=0.18,
                        linewidth=2.3,
                        zorder=3,
                    )

                if len(vals) >= 2:
                    ax.plot(
                        [x - 0.11, x + 0.11],
                        [np.nanmedian(vals), np.nanmedian(vals)],
                        color="#111111",
                        lw=0.55,
                        alpha=0.65,
                        zorder=4,
                    )

                q = sig.get((sp, int(order), metric))
                p_value = raw_p.get((sp, int(order), metric))
                if q is not None or (
                    (sp == "Zebrafish" and order == 0 and metric == "FCV")
                    or (sp == "Drosophila" and order == 0 and metric in {"Post-DCA", "Pre-DCA"})
                ):
                    label = p_to_stars(q) if q is not None else f"p={p_value:.3f}".replace("p=0.", "p=.")
                    if "_stat_labels" not in ax.__dict__:
                        ax.__dict__["_stat_labels"] = []
                    ax.__dict__["_stat_labels"].append((x, label, q is None))

        if all_plot_values:
            vals_for_scale = np.asarray(all_plot_values, dtype=float)
            ymin = float(np.nanpercentile(vals_for_scale, 1.0))
            ymax = float(np.nanpercentile(vals_for_scale, 99.0))
            yr = max(ymax - ymin, 1e-6)
            ymin = max(float(np.nanmin(vals_for_scale)) - 0.02 * yr, ymin - 0.12 * yr)
            ymax = min(float(np.nanmax(vals_for_scale)) + 0.02 * yr, ymax + 0.24 * yr)
            ax.set_ylim(ymin, ymax)
            stat_y = ymax - 0.13 * (ymax - ymin)
            for x, label, is_p_label in ax.__dict__.get("_stat_labels", []):
                ax.text(
                    x,
                    stat_y,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=5.9 if is_p_label else 7.6,
                    fontweight="heavy",
                    color="#111111",
                    clip_on=False,
                    zorder=7,
                )
        ax.axhline(0, color="#777777", lw=0.50, ls=":", zorder=0)
        ax.set_xlim(-0.65, len(orders) - 0.35)
        ax.set_title(metric, fontsize=8.4, fontweight="bold", pad=3)
        ax.set_ylabel("within-species metric z", fontsize=6.3)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([GROUP_LABELS.get(int(order), str(order)) for order in orders], rotation=0, ha="center", fontsize=5.8)
        for tick, order in zip(ax.get_xticklabels(), orders, strict=False):
            tick.set_color(GROUP_COLORS.get(int(order), "#555555"))
            tick.set_fontweight("bold" if order in HIGHLIGHT_GROUPS else "normal")
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.grid(axis="y", color="#E9E9E9", lw=0.40, zorder=-5)
        ax.tick_params(axis="both", direction="out", pad=1.2, length=2.0, width=0.45, labelsize=5.6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_color("#BFBFBF")
        ax.spines["bottom"].set_linewidth(0.45)
        ax.spines["left"].set_linewidth(0.45)

        if panel_idx == 0:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker=species_markers[sp],
                    color="none",
                    markerfacecolor=SPECIES_COLORS.get(sp, "#777777"),
                    markeredgecolor="#222222",
                    markersize=5.0,
                    label=sp,
                )
                for sp in SPECIES
            ]
            ax.legend(
                handles=handles,
                loc="upper right",
                bbox_to_anchor=(1.0, 1.22),
                frameon=False,
                ncol=3,
                columnspacing=0.9,
                handletextpad=0.3,
                borderaxespad=0.0,
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Draw Figure 7 structure-function synthesis.")
    parser.add_argument(
        "--recompute-stats",
        action="store_true",
        help="Recompute permutation/bootstrap statistics instead of reading cached CSV files.",
    )
    args = parser.parse_args()

    tables, effects, signature, group_perm = load_or_compute_statistics(recompute=args.recompute_stats)

    fig, axes = plt.subplots(3, 1, figsize=(6.9, 6.35), sharex=False)
    fig.subplots_adjust(left=0.095, right=0.99, bottom=0.095, top=0.93, hspace=0.50)
    draw_measure_comparison_panels(axes, tables, group_perm)

    fig.text(0.985, 0.022, "* q<0.05, ** q<0.01, *** q<0.001 (BH-FDR)", ha="right", va="center", fontsize=5.7)
    fig.text(0.02, 0.985, "a", ha="left", va="top", fontsize=8.2, fontweight="bold")
    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_EFFECTS}")
    print(f"wrote {OUT_GROUP}")
    print(f"wrote {OUT_GROUP_PERM}")


if __name__ == "__main__":
    main()
