"""Summarize four SC measures against FCV across species.

Measures:
- Post-DCA
- Pre-DCA
- Modularity
- log(out/in)
- OO fraction

The figure mirrors the geometric FCV-variant figure: species are rows, the
left panel is a node-wise feature heat strip sorted by FCV, and the remaining
columns show FCV-vs-SC scatter plots.
"""

from __future__ import annotations

from pathlib import Path
import itertools

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage

from figure1_cross_species_fc_dynamics_measures import (
    ANATOMY_COLORS,
    ANATOMY_MARKERS,
    ANATOMY_ORDER,
    add_functional_classes,
)


ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "data" / "final_summary_tables"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

ZSCORE_FCV = TAB / "highpass_ce_zf_plot_measures_recording_zscore_node_summary.csv"
RECORDING_FCV = TAB / "highpass_ce_zf_plot_measures_recording_node.csv"
TARGET_SUMMARY = TAB / "edge_target_dca_distribution_summary_by_unit.csv"
TARGET_SAMPLE = TAB / "edge_target_dca_distribution_sample_by_unit.csv"
PRE_SUMMARY = TAB / "pre_edge_source_dca_distribution_summary_by_unit.csv"
SUBSC_MOD = TAB / "region_subsc_modularity_fcv_values.csv"
CE_LOG = TAB / "celegans_out_in_log_ratio_neuron_values.csv"
FLY_LOG = TAB / "drosophila_clustering_out_in_ratio_edge_fcv_values.csv"
ZF_LOG = TAB / "zebrafish_region_clustering_inter_out_in_ratio_edge_fcv_values.csv"
CE_SC = ROOT / "data" / "source_inputs" / "external_processed" / "fcv_postdca_raw_recompute" / "out_data" / "celegans" / "post_dca" / "celegans_sc_matrix_full297_no_diagonal.npz"

OUT_FIG = FIG / "sc_four_measures_vs_fcv_all_species.png"
OUT_VALUES = TAB / "sc_four_measures_vs_fcv_all_species_values.csv"
OUT_REGION_VALUES = TAB / "sc_four_measures_vs_fcv_all_species_region_mean_values.csv"
OUT_CORR = TAB / "sc_four_measures_vs_fcv_all_species_correlations.csv"
ZF_SUBJECT_POSTDCA = ROOT / "data" / "source_inputs" / "external_processed" / "fcv_postdca_raw_recompute" / "out_data" / "zebrafish" / "post_dca_rank1" / "zebrafish_rank1_subject_region_post_dca.csv"
ZF_SUBJECT_MOD = TAB / "region_subsc_modularity_subject_values.csv"
ZF_SUBJECT_LOG = TAB / "zebrafish_subject_region_clustering_inter_out_in_ratio.csv"
OO_II_VALUES = TAB / "oo_fraction_recomputed_values_by_species.csv"

SPECIES = ["C. elegans", "Drosophila", "Zebrafish"]
SPECIES_COLORS = {
    "C. elegans": "#4F6D8A",
    "Drosophila": "#B85C5A",
    "Zebrafish": "#5E8C61",
}
MEASURES = [
    ("PostDCA", "Post-DCA"),
    ("PreDCA", "Pre-DCA"),
    ("Modularity", "Modularity"),
    ("LogOutIn", "log(out/in)"),
    ("OO_fraction", "OO fraction"),
]
HEAT_MEASURES = MEASURES
YRANGE_PERCENTILES = {
    "C. elegans": (5, 95),
    "Drosophila": (8, 92),
    "Zebrafish": (8, 92),
}
TARGET_Y_LIMITS = {
    "Drosophila": (-0.08, 0.02),
}
GROUP_Y_LIMITS = {
    ("Zebrafish", "PostDCA"): (-0.060, 0.035),
    ("Zebrafish", "PreDCA"): (-0.040, 0.020),
}


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.6,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.2,
        "ytick.major.size": 2.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 300,
    }
)


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return out
    sd = np.nanstd(values[finite])
    if sd <= 1e-12:
        out[finite] = 0.0
    else:
        out[finite] = (values[finite] - np.nanmean(values[finite])) / sd
    return out


def ce_ego_modularity() -> pd.DataFrame:
    """Compute 1-hop structural neighborhood modularity for C. elegans neurons."""
    raw = np.load(CE_SC, allow_pickle=True)
    neurons = raw["neuron"].astype(str)
    sc = np.asarray(raw["sc"], dtype=float)
    und = sc + sc.T
    np.fill_diagonal(und, 0.0)
    rows = []
    for i, node in enumerate(neurons):
        nbr = np.flatnonzero(und[i] > 0)
        if len(nbr) < 4:
            q = np.nan
            n_edges = 0
            n_modules = 0
        else:
            sub = und[np.ix_(nbr, nbr)].copy()
            graph = nx.Graph()
            graph.add_nodes_from(range(len(nbr)))
            for a in range(len(nbr)):
                for b in range(a + 1, len(nbr)):
                    if sub[a, b] > 0:
                        graph.add_edge(a, b, weight=float(sub[a, b]))
            n_edges = graph.number_of_edges()
            if n_edges == 0:
                q = np.nan
                n_modules = 0
            else:
                communities = list(nx.algorithms.community.greedy_modularity_communities(graph, weight="weight"))
                q = nx.algorithms.community.quality.modularity(graph, communities, weight="weight")
                n_modules = len(communities)
        rows.append(
            {
                "species": "C. elegans",
                "node": node,
                "Modularity": q,
                "modularity_definition": "1-hop induced structural neighborhood modularity Q",
                "modularity_edges": n_edges,
                "modularity_n_modules": n_modules,
            }
        )
    return pd.DataFrame(rows)


def load_values() -> pd.DataFrame:
    fcv = pd.read_csv(ZSCORE_FCV)[["species", "node", "EdgeStdFCV"]].copy()
    fcv = fcv[fcv["species"].isin(SPECIES)]

    # Match the exact unit set used by the existing Post-DCA and Pre-DCA
    # distribution summary figures.  This avoids adding atlas nodes that have
    # FCV but were not present in those DCA summary panels.
    post_units = pd.read_csv(TARGET_SUMMARY)[["species", "node"]].drop_duplicates()
    pre_units = pd.read_csv(PRE_SUMMARY)[["species", "node"]].drop_duplicates()
    dca_units = post_units.merge(pre_units, on=["species", "node"], how="inner")
    dca_units = dca_units[dca_units["species"].isin(SPECIES)]
    fcv = fcv.merge(dca_units, on=["species", "node"], how="inner")

    bundled_features = pd.read_csv(OO_II_VALUES).replace([np.inf, -np.inf], np.nan)
    bundled_features = bundled_features[
        (
            bundled_features["species"].isin(["C. elegans", "Drosophila"])
            & bundled_features["oo_level"].isin(["neuron", "side-aware region"])
        )
        | (bundled_features["species"].eq("Zebrafish") & bundled_features["oo_level"].eq("region mean"))
    ].drop_duplicates(["species", "node"])

    post = pd.read_csv(TARGET_SUMMARY)[["species", "node", "PostDCA"]].drop_duplicates(["species", "node"])
    pre = bundled_features[["species", "node", "PreDCA"]].drop_duplicates(["species", "node"])
    dca = post.merge(pre, on=["species", "node"], how="left")

    mod = bundled_features[
        ["species", "node", "Modularity", "modularity_definition", "modularity_edges", "modularity_n_modules"]
    ].drop_duplicates(["species", "node"])
    logout = bundled_features[["species", "node", "LogOutIn"]].drop_duplicates(["species", "node"])
    oo_region = bundled_features[["species", "node", "OO_fraction"]].drop_duplicates(["species", "node"])

    values = fcv.merge(dca, on=["species", "node"], how="left")
    values = values.merge(mod, on=["species", "node"], how="left")
    values = values.merge(logout, on=["species", "node"], how="left")
    values = values.merge(oo_region, on=["species", "node"], how="left")
    values = add_functional_classes(values)
    values = values.replace([np.inf, -np.inf], np.nan)
    values.to_csv(OUT_REGION_VALUES, index=False)
    return values


def load_target_summary_with_zscore() -> pd.DataFrame:
    summary = pd.read_csv(TARGET_SUMMARY)
    z = pd.read_csv(ZSCORE_FCV)[["species", "node", "EdgeStdFCV"]].rename(
        columns={"EdgeStdFCV": "EdgeStdFCV_zscore"}
    )
    summary = summary.drop(columns=["EdgeStdFCV_zscore"], errors="ignore").merge(
        z,
        on=["species", "node"],
        how="left",
    )
    if "EdgeStdFCV" in summary.columns:
        summary = summary.rename(columns={"EdgeStdFCV": "EdgeStdFCV_raw"})
    summary["EdgeStdFCV"] = summary["EdgeStdFCV_zscore"]
    summary["fcv_normalization"] = "within_recording_zscore_then_node_mean"
    return summary


def trend_stats(values: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species in SPECIES:
        sdf = values[values["species"].eq(species)]
        for measure, _ in MEASURES:
            sub = sdf[[measure, "EdgeStdFCV"]].dropna()
            if len(sub) >= 4 and sub[measure].nunique() > 1 and sub["EdgeStdFCV"].nunique() > 1:
                pr = stats.pearsonr(sub[measure], sub["EdgeStdFCV"])
                sr = stats.spearmanr(sub[measure], sub["EdgeStdFCV"])
                slope, intercept = np.polyfit(sub[measure].to_numpy(float), sub["EdgeStdFCV"].to_numpy(float), 1)
            else:
                pr = sr = None
                slope = intercept = np.nan
            rows.append(
                {
                    "species": species,
                    "measure": measure,
                    "n": int(len(sub)),
                    "pearson_r": float(pr.statistic) if pr else np.nan,
                    "pearson_p": float(pr.pvalue) if pr else np.nan,
                    "spearman_rho": float(sr.statistic) if sr else np.nan,
                    "spearman_p": float(sr.pvalue) if sr else np.nan,
                    "ols_slope": float(slope),
                    "ols_intercept": float(intercept),
                }
            )
    corr = pd.DataFrame(rows)
    corr.to_csv(OUT_CORR, index=False)
    return corr


def expand_zebrafish_subject_points(values: pd.DataFrame) -> pd.DataFrame:
    """Use subject x region rows for Zebrafish FCV and SC predictors."""
    non_zf = values[~values["species"].eq("Zebrafish")].copy()
    rec = pd.read_csv(RECORDING_FCV).replace([np.inf, -np.inf], np.nan)
    rec = rec[rec["species"].eq("Zebrafish")][["species", "recording_id", "node", "EdgeStdFCV"]].copy()
    rec["Subject"] = rec["recording_id"].astype(str).str.extract(r"subject_(\d+)").astype(float).astype("Int64")
    zframes = []
    for _, group in rec.groupby(["species", "recording_id"], sort=False):
        out = group.copy()
        out["EdgeStdFCV"] = zscore(out["EdgeStdFCV"].to_numpy(float))
        zframes.append(out)
    zf_rec = pd.concat(zframes, ignore_index=True) if zframes else rec

    zf_dca = pd.read_csv(ZF_SUBJECT_POSTDCA).replace([np.inf, -np.inf], np.nan)
    zf_dca = zf_dca.rename(
        columns={
            "Region": "node",
            "Rank1PostDCA": "PostDCA",
            "Rank1PreDCA": "PreDCA",
        }
    )[["Subject", "node", "PostDCA", "PreDCA"]]

    zf_mod = pd.read_csv(ZF_SUBJECT_MOD).replace([np.inf, -np.inf], np.nan)
    zf_mod = zf_mod[zf_mod["species"].eq("Zebrafish")].rename(
        columns={
            "subject": "Subject",
            "subsc_modularity_q": "Modularity",
            "subsc_edges": "modularity_edges",
            "subsc_n_modules": "modularity_n_modules",
        }
    )[["Subject", "node", "Modularity", "modularity_edges", "modularity_n_modules"]]
    zf_mod["Subject"] = zf_mod["Subject"].astype(float).astype("Int64")
    zf_mod["modularity_definition"] = "subject-specific within-region cell-level induced sub-SC modularity Q"

    zf_log = pd.read_csv(ZF_SUBJECT_LOG).replace([np.inf, -np.inf], np.nan)
    zf_log = zf_log.rename(
        columns={
            "Region": "node",
            "inter_out_in_log_ratio": "LogOutIn",
        }
    )[["Subject", "node", "LogOutIn"]]

    zf_oo = pd.read_csv(OO_II_VALUES).replace([np.inf, -np.inf], np.nan)
    zf_oo = zf_oo[
        zf_oo["species"].eq("Zebrafish")
        & zf_oo["oo_level"].eq("subject-region")
    ][["Subject", "node", "OO_fraction"]].drop_duplicates(["Subject", "node"])

    zf = zf_rec.merge(zf_dca, on=["Subject", "node"], how="inner")
    zf = zf.merge(zf_mod, on=["Subject", "node"], how="left")
    zf = zf.merge(zf_log, on=["Subject", "node"], how="left")
    zf = zf.merge(zf_oo, on=["Subject", "node"], how="left")
    zf["recording_id"] = "subject_" + zf["Subject"].astype(str)
    zf["species"] = "Zebrafish"
    non_zf["recording_id"] = "node_mean"
    out = pd.concat([non_zf, zf], ignore_index=True, sort=False)
    out = out.drop(
        columns=["class_order", "class_label", "class_short", "fine_class", "anatomy_group"],
        errors="ignore",
    )
    out = add_functional_classes(out)
    return out


def p_text(p: float) -> str:
    if not np.isfinite(p):
        return "P=nan"
    if p < 1e-3:
        return "P<0.001"
    return f"P={p:.3f}"


def holm_adjust(pvals: list[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[finite])]
    m = len(order)
    running = 0.0
    for rank, original_idx in enumerate(order):
        adjusted = (m - rank) * p[original_idx]
        running = max(running, adjusted)
        out[original_idx] = min(running, 1.0)
    return out


def permutation_group_tests(
    values_by_group: list[np.ndarray],
    n_perm: int = 20000,
    seed: int = 20260526,
) -> tuple[float, list[dict[str, object]]]:
    """Global label permutation preserving group sizes, with Holm-corrected pairwise bars."""
    clean_groups = [np.asarray(vals, dtype=float)[np.isfinite(vals)] for vals in values_by_group]
    sizes = np.asarray([len(vals) for vals in clean_groups], dtype=int)
    if np.sum(sizes > 0) < 2:
        return np.nan, []
    pooled = np.concatenate([vals for vals in clean_groups if len(vals)])
    labels = np.concatenate([np.full(len(vals), idx, dtype=int) for idx, vals in enumerate(clean_groups) if len(vals)])
    valid_group_ids = np.asarray([idx for idx, vals in enumerate(clean_groups) if len(vals)], dtype=int)
    n = len(pooled)
    grand = float(np.mean(pooled))

    def global_stat(vals: np.ndarray, lab: np.ndarray) -> float:
        ss_between = 0.0
        ss_total = float(np.sum((vals - np.mean(vals)) ** 2))
        if ss_total <= 1e-15:
            return 0.0
        for group_id in valid_group_ids:
            group_vals = vals[lab == group_id]
            if len(group_vals):
                ss_between += len(group_vals) * float((np.mean(group_vals) - np.mean(vals)) ** 2)
        return ss_between / ss_total

    pairs = list(itertools.combinations(valid_group_ids.tolist(), 2))
    obs_global = global_stat(pooled, labels)
    obs_diffs = []
    means = {idx: float(np.mean(clean_groups[idx])) for idx in valid_group_ids}
    for i, j in pairs:
        obs_diffs.append(abs(means[i] - means[j]))
    obs_diffs_arr = np.asarray(obs_diffs, dtype=float)

    rng = np.random.default_rng(seed)
    null_global = np.empty(n_perm, dtype=float)
    shuffled = labels.copy()
    null_pair_diffs = np.empty((n_perm, len(pairs)), dtype=float)
    for perm_idx in range(n_perm):
        rng.shuffle(shuffled)
        null_global[perm_idx] = global_stat(pooled, shuffled)
        perm_means = {}
        for group_id in valid_group_ids:
            vals = pooled[shuffled == group_id]
            perm_means[group_id] = float(np.mean(vals)) if len(vals) else np.nan
        for pair_idx, (i, j) in enumerate(pairs):
            null_pair_diffs[perm_idx, pair_idx] = abs(perm_means[i] - perm_means[j])

    global_p = (1.0 + np.sum(null_global >= obs_global)) / (n_perm + 1.0)
    raw_pair_p = [
        (1.0 + np.sum(null_pair_diffs[:, pair_idx] >= obs)) / (n_perm + 1.0)
        for pair_idx, obs in enumerate(obs_diffs_arr)
    ]
    holm_pair_p = holm_adjust(raw_pair_p)
    pair_rows = []
    for (i, j), obs, raw_p, holm_p in zip(pairs, obs_diffs_arr, raw_pair_p, holm_pair_p, strict=False):
        pair_rows.append(
            {
                "group_i": int(i),
                "group_j": int(j),
                "observed_abs_mean_diff": float(obs),
                "permutation_p": float(raw_p),
                "holm_p": float(holm_p),
                "n_perm": int(n_perm),
            }
        )
    return float(global_p), pair_rows


def add_pairwise_sig_bars(
    ax: plt.Axes,
    values_by_group: list[np.ndarray],
    seed: int,
    alpha: float = 0.05,
) -> tuple[float, list[dict[str, object]]]:
    global_p, pair_rows = permutation_group_tests(values_by_group, seed=seed)
    sig = [
        ((row["group_i"], row["group_j"]), float(row["holm_p"]))
        for row in pair_rows
        if np.isfinite(row["holm_p"]) and row["holm_p"] < alpha
    ]
    sig = sorted(sig, key=lambda item: (item[0][1] - item[0][0], item[1]))
    if sig:
        y_min, y_max = ax.get_ylim()
        yr = y_max - y_min if y_max > y_min else 1.0
        step = yr * 0.050
        bar_h = step * 0.22
        for lvl, ((i, j), p) in enumerate(sig):
            y = y_max + yr * 0.035 + lvl * step
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*"
            ax.plot([i, i, j, j], [y, y + bar_h, y + bar_h, y], lw=0.55, c="#333333", clip_on=False)
            ax.text((i + j) / 2, y + bar_h, star, ha="center", va="bottom", fontsize=5.6, clip_on=False)
        ax.set_ylim(y_min, y_max)
    return global_p, pair_rows


def plot_group_boxstrip(
    ax: plt.Axes,
    df: pd.DataFrame,
    species: str,
    measure: str,
    label: str,
    show_ylabel: bool,
) -> list[dict[str, object]]:
    order = list(ANATOMY_ORDER[species])
    if any(group not in ANATOMY_ORDER[species] for group in df["anatomy_group"].dropna().unique()):
        order.append("other")

    rows: list[dict[str, object]] = []
    values_by_group = []
    rng_seed = sum(ord(ch) for ch in f"sc-{species}-{measure}") % (2**32)
    rng = np.random.default_rng(rng_seed)
    for idx, group in enumerate(order):
        vals = (
            df.loc[df["anatomy_group"].eq(group), measure]
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
                marker=ANATOMY_MARKERS.get(group, "o"),
                facecolor=ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]),
                edgecolor="#222222",
                alpha=0.58,
                linewidth=0.25,
                rasterized=True,
                zorder=3,
            )
        rows.append(
            {
                "species": species,
                "measure": measure,
                "anatomy_group": group,
                "n": int(len(vals)),
                "mean": float(np.nanmean(vals)) if len(vals) else np.nan,
                "median": float(np.nanmedian(vals)) if len(vals) else np.nan,
                "std": float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else np.nan,
            }
        )

    box_values = [vals if len(vals) else np.asarray([np.nan]) for vals in values_by_group]
    box = ax.boxplot(
        box_values,
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
    for idx, vals in enumerate(values_by_group):
        if len(vals) == 0:
            y_min, y_max = ax.get_ylim()
            ax.text(
                idx,
                y_min + (y_max - y_min) * 0.04,
                "n=0",
                ha="center",
                va="bottom",
                fontsize=4.9,
                color="#777777",
                rotation=90,
            )

    valid_groups = [vals for vals in values_by_group if len(vals) > 0]
    if len(valid_groups) >= 2:
        try:
            stat, pval = stats.kruskal(*valid_groups)
        except ValueError:
            stat, pval = np.nan, np.nan
    else:
        stat, pval = np.nan, np.nan
    if (species, measure) in GROUP_Y_LIMITS:
        ax.set_ylim(*GROUP_Y_LIMITS[(species, measure)])
    perm_seed = sum(ord(ch) for ch in f"perm-{species}-{measure}") % (2**32)
    global_perm_p, pair_rows = add_pairwise_sig_bars(ax, values_by_group, seed=perm_seed, alpha=0.05)
    ax.text(
        0.02,
        0.98,
        f"perm {p_text(global_perm_p)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.1,
    )
    rows.append(
        {
            "species": species,
            "measure": measure,
            "anatomy_group": "__global_permutation__",
            "n": int(sum(len(vals) for vals in values_by_group)),
            "kruskal_h": float(stat) if np.isfinite(stat) else np.nan,
            "kruskal_p": float(pval) if np.isfinite(pval) else np.nan,
            "global_permutation_p": global_perm_p,
            "n_perm": 20000,
        }
    )

    for pair_row in pair_rows:
        gi = int(pair_row["group_i"])
        gj = int(pair_row["group_j"])
        rows.append(
            {
                "species": species,
                "measure": measure,
                "anatomy_group": "__pairwise_permutation_holm__",
                "group_i": order[gi],
                "group_j": order[gj],
                "n_i": int(len(values_by_group[gi])),
                "n_j": int(len(values_by_group[gj])),
                "observed_abs_mean_diff": pair_row["observed_abs_mean_diff"],
                "permutation_p": pair_row["permutation_p"],
                "holm_p": pair_row["holm_p"],
                "n_perm": pair_row["n_perm"],
            }
        )

    ax.axhline(0, color="#777777", lw=0.55, ls=":", zorder=0)
    ax.set_title(label, fontsize=7.8, pad=3)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=38, ha="right", fontsize=5.4)
    for tick, group in zip(ax.get_xticklabels(), order, strict=False):
        tick.set_color(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
        tick.set_fontweight("bold")
    if show_ylabel:
        ax.set_ylabel("SC measure value", fontsize=7.0)
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="both", direction="out", pad=1.3, length=2.2, width=0.55)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return rows


def draw_species_boxgrid(fig: plt.Figure, slot, df: pd.DataFrame, species: str) -> tuple[plt.Axes, list[dict[str, object]]]:
    panel_ax = fig.add_subplot(slot)
    panel_ax.set_axis_off()
    panel_ax.text(
        -0.05,
        0.50,
        species,
        transform=panel_ax.transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold",
        color=SPECIES_COLORS[species],
    )
    inner = gridspec.GridSpecFromSubplotSpec(1, len(MEASURES), subplot_spec=slot, wspace=0.16)
    rows: list[dict[str, object]] = []
    sdf = df[df["species"].eq(species)].copy()
    for idx, (measure, label) in enumerate(MEASURES):
        ax = fig.add_subplot(inner[0, idx])
        rows.extend(plot_group_boxstrip(ax, sdf, species, measure, label, show_ylabel=idx == 0))
    return panel_ax, rows


def target_distribution_limits(sub: pd.DataFrame, samp: pd.DataFrame, species: str) -> tuple[float, float]:
    vals = samp["target_DCA"].to_numpy(float)
    vals = vals[np.isfinite(vals)]
    if species in TARGET_Y_LIMITS:
        return TARGET_Y_LIMITS[species]
    lo_p, hi_p = YRANGE_PERCENTILES[species]
    ymin, ymax = np.nanpercentile(vals, [lo_p, hi_p])
    anchors = sub[["PostDCA", "edge_target_dca_q25", "edge_target_dca_q75"]].to_numpy(float).ravel()
    anchors = anchors[np.isfinite(anchors)]
    if len(anchors):
        ymin = min(ymin, float(np.nanpercentile(anchors, 1)))
        ymax = max(ymax, float(np.nanpercentile(anchors, 99)))
    pad = max((ymax - ymin) * 0.08, 0.0015)
    return ymin - pad, ymax + pad


def draw_post_distribution(
    ax_bar: plt.Axes,
    ax: plt.Axes,
    target_summary: pd.DataFrame,
    target_sample: pd.DataFrame,
    species: str,
) -> None:
    sub = target_summary[target_summary["species"].eq(species) & target_summary["EdgeStdFCV"].notna()].copy()
    sub = add_functional_classes(sub)
    sub = sub.sort_values("EdgeStdFCV", ascending=True).reset_index(drop=True)
    sub["x"] = np.arange(len(sub), dtype=float)
    x_lookup = dict(zip(sub["node"], sub["x"]))
    anatomy_lookup = dict(zip(sub["node"], sub["anatomy_group"]))

    samp = target_sample[target_sample["species"].eq(species) & target_sample["node"].isin(x_lookup)].copy()
    samp["x"] = samp["node"].map(x_lookup)
    samp["anatomy_group"] = samp["node"].map(anatomy_lookup).fillna("other")

    ymin, ymax = target_distribution_limits(sub, samp, species)
    fcv = sub["EdgeStdFCV"].to_numpy(float)
    fcv_norm = (fcv - np.nanmin(fcv)) / max(np.nanmax(fcv) - np.nanmin(fcv), 1e-12)
    im = ax_bar.imshow(
        fcv_norm[np.newaxis, :],
        aspect="auto",
        extent=(-0.5, len(sub) - 0.5, 0, 1),
        cmap="viridis",
        interpolation="nearest",
    )
    ax_bar.set_yticks([])
    ax_bar.set_xticks([])
    ax_bar.set_ylabel("z-FCV", fontsize=6, rotation=0, ha="right", va="center", labelpad=11)
    ax_bar.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax_bar.text(-0.5, 1.08, "low", fontsize=5.4, ha="left", va="bottom")
    ax_bar.text(len(sub) - 0.5, 1.08, "high", fontsize=5.4, ha="right", va="bottom")

    plot_samp = samp.dropna(subset=["x", "target_DCA"]).copy()
    plot_samp = plot_samp[(plot_samp["target_DCA"] >= ymin) & (plot_samp["target_DCA"] <= ymax)]
    if len(plot_samp) > 60000:
        plot_samp = plot_samp.sample(60000, random_state=20260524)
    jitter_sd = 0.055 if len(sub) <= 60 else 0.035
    jitter = np.random.default_rng(20260524).normal(0, jitter_sd, size=len(plot_samp))
    plot_samp["x_jitter"] = plot_samp["x"].to_numpy(float) + jitter
    for anatomy_group, group in plot_samp.groupby("anatomy_group", sort=True):
        anatomy_group = str(anatomy_group)
        ax.scatter(
            group["x_jitter"].to_numpy(float),
            group["target_DCA"].to_numpy(float),
            s=3.0 if species == "C. elegans" else 1.8,
            marker=ANATOMY_MARKERS.get(anatomy_group, "o"),
            color=ANATOMY_COLORS.get(anatomy_group, ANATOMY_COLORS["other"]),
            alpha=0.20 if species == "C. elegans" else 0.14,
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
                linewidth=0.42,
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
                s=20,
                marker="D",
                color="#C44E52",
                edgecolor="white",
                linewidth=0.4,
                zorder=6,
            )
    if mean_x:
        ax.plot(mean_x, mean_y, color="#C44E52", lw=1.0, alpha=0.9, zorder=5)

    ax.axhline(0, color="#777777", lw=0.6, ls="--")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-0.5, len(sub) - 0.5)
    ax.set_xticks([])
    ax.set_ylabel("Target DCA", fontsize=6.5)
    ax.set_title(species, color=SPECIES_COLORS[species], fontsize=7.5, pad=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=5.5, width=0.6, length=2)
    return im


def draw_sc_heatmap(
    dend_ax: plt.Axes,
    bar_ax: plt.Axes,
    heat_ax: plt.Axes,
    df: pd.DataFrame,
    species: str,
) -> None:
    tmp = df[df["species"].eq(species)].copy()
    measure_cols = [col for col, _ in HEAT_MEASURES]
    tmp = tmp.dropna(subset=["EdgeStdFCV"]).reset_index(drop=True)
    matrix = np.vstack([zscore(tmp[col].to_numpy(float)) for col, _ in HEAT_MEASURES])
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)

    if len(tmp) >= 3:
        z = linkage(matrix.T, method="ward")
        info = dendrogram(
            z,
            ax=dend_ax,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#333333",
        )
        order = np.asarray(info["leaves"], dtype=int)
        tmp = tmp.iloc[order].reset_index(drop=True)
        matrix = matrix[:, order]
        dend_ax.axis("off")
    else:
        tmp = tmp.sort_values("EdgeStdFCV", ascending=False).reset_index(drop=True)
        matrix = np.vstack([zscore(tmp[col].to_numpy(float)) for col, _ in HEAT_MEASURES])
        matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
        dend_ax.axis("off")

    anatomy_order = ANATOMY_ORDER[species]
    group_to_i = {group: idx for idx, group in enumerate(anatomy_order)}
    class_idx = np.asarray([group_to_i.get(group, len(anatomy_order)) for group in tmp["anatomy_group"]])
    color_list = [ANATOMY_COLORS[group] for group in anatomy_order] + [ANATOMY_COLORS["other"]]
    bar_ax.imshow(
        class_idx[None, :],
        aspect="auto",
        cmap=ListedColormap(color_list),
        vmin=0,
        vmax=len(color_list) - 1,
        interpolation="nearest",
    )
    bar_ax.set_axis_off()

    vmax = max(1.0, float(np.nanpercentile(np.abs(matrix), 98)))
    heat_ax.imshow(matrix, aspect="auto", cmap="PiYG_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    heat_ax.set_yticks(np.arange(len(HEAT_MEASURES)))
    heat_ax.set_yticklabels([label for _, label in HEAT_MEASURES], fontsize=5.4)
    heat_ax.yaxis.tick_right()
    heat_ax.tick_params(axis="y", labelright=True, labelleft=False, length=0, pad=1.0)
    heat_ax.set_xticks(np.arange(len(tmp)))
    heat_ax.set_xticklabels(tmp["node"].astype(str).tolist(), rotation=90, ha="center", va="top")
    label_size = 2.3 if len(tmp) > 80 else 3.0 if len(tmp) > 45 else 4.0
    for tick, group in zip(heat_ax.get_xticklabels(), tmp["anatomy_group"], strict=False):
        tick.set_fontsize(label_size)
        tick.set_color(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
    for y in np.arange(-0.5, len(HEAT_MEASURES), 1):
        heat_ax.axhline(y, color="white", lw=0.35)
    for spine in heat_ax.spines.values():
        spine.set_linewidth(0.55)
        spine.set_color("#222222")
    heat_ax.set_ylabel(species, rotation=0, ha="right", va="center", labelpad=31, fontweight="bold")


def draw_scatter(ax: plt.Axes, df: pd.DataFrame, corr: pd.DataFrame, species: str, measure: str, xlabel: str) -> None:
    sub = df[df["species"].eq(species)].dropna(subset=[measure, "EdgeStdFCV"]).copy()
    for anatomy_group, group in sub.groupby("anatomy_group", dropna=False):
        anatomy_group = str(anatomy_group) if pd.notna(anatomy_group) else "other"
        ax.scatter(
            group[measure],
            group["EdgeStdFCV"],
            s=22,
            marker=ANATOMY_MARKERS.get(anatomy_group, "o"),
            facecolor=ANATOMY_COLORS.get(anatomy_group, ANATOMY_COLORS["other"]),
            edgecolor="#202020",
            linewidth=0.35,
            alpha=0.88,
            zorder=2,
        )
    stat = corr[(corr["species"].eq(species)) & (corr["measure"].eq(measure))]
    if not stat.empty and len(sub) >= 4 and sub[measure].nunique() > 1:
        row = stat.iloc[0]
        x = sub[measure].to_numpy(float)
        grid = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ax.plot(grid, row["ols_slope"] * grid + row["ols_intercept"], color="#202020", lw=0.9, zorder=3)
        ax.text(
            0.04,
            0.96,
            f"n={int(row['n'])}\nr={row['pearson_r']:.2f}\n{p_text(row['pearson_p'])}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=5.8,
        )
    else:
        ax.text(0.04, 0.96, f"n={len(sub)}", transform=ax.transAxes, ha="left", va="top", fontsize=5.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("FCV z" if measure == "PostDCA" else "")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#E7E7E7", lw=0.35)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))


def add_legend(fig: plt.Figure) -> None:
    handles = [
        Line2D(
            [0],
            [0],
            marker=ANATOMY_MARKERS.get(group, "o"),
            linestyle="none",
            markerfacecolor=ANATOMY_COLORS[group],
            markeredgecolor="#202020",
            markeredgewidth=0.35,
            markersize=4.5,
            label=group,
        )
        for group in [
            "olf./chemo",
            "sensory neurons",
            "interneurons",
            "motor / command",
            "state-modulatory",
            "olfactory neuropil",
            "optic neuropil",
            "mushroom body",
            "central complex",
            "protocerebrum",
            "premotor interface",
            "telencephalon",
            "diencephalon",
            "mesencephalon",
            "hindbrain / motor",
            "state system",
            "other",
        ]
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=8,
        frameon=False,
        fontsize=5.6,
        handletextpad=0.35,
        columnspacing=0.9,
    )


def main() -> None:
    values = load_values()
    plot_values = expand_zebrafish_subject_points(values)
    plot_values.to_csv(OUT_VALUES, index=False)
    corr = trend_stats(plot_values)

    fig = plt.figure(figsize=(17.8, 8.2), constrained_layout=True)
    gs = gridspec.GridSpec(
        3,
        2,
        figure=fig,
        width_ratios=[1.24, 2.76],
        wspace=0.04,
        hspace=0.18,
    )
    heat_axes = []
    right_axes = []
    anatomy_rows: list[dict[str, object]] = []
    for row, species in enumerate(SPECIES):
        left_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[row, 0], height_ratios=[0.35, 0.06, 0.74], hspace=0.03)
        dend_ax = fig.add_subplot(left_gs[0, 0])
        bar_ax = fig.add_subplot(left_gs[1, 0])
        hax = fig.add_subplot(left_gs[2, 0])
        heat_axes.append(hax)
        draw_sc_heatmap(dend_ax, bar_ax, hax, values, species)
        if row == 0:
            dend_ax.set_title("SC feature fingerprint", pad=2)
        panel_ax, rows = draw_species_boxgrid(fig, gs[row, 1], plot_values, species)
        right_axes.append(panel_ax)
        anatomy_rows.extend(rows)

    for label, ax in [("A", heat_axes[0]), ("B", right_axes[0]), ("C", right_axes[1]), ("D", right_axes[2])]:
        ax.text(-0.06, 1.04, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")
    add_legend(fig)

    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight")
    plt.close(fig)
    anatomy_path = TAB / "sc_four_measures_vs_fcv_all_species_anatomy_group_stats.csv"
    pd.DataFrame(anatomy_rows).to_csv(anatomy_path, index=False)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_VALUES}")
    print(f"wrote {OUT_REGION_VALUES}")
    print(f"wrote {OUT_CORR}")
    print(f"wrote {anatomy_path}")


if __name__ == "__main__":
    main()
