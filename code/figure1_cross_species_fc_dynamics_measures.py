"""Plot FCV-adjacent measures by anatomical group for all species."""

from __future__ import annotations

from pathlib import Path
import itertools

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

import figure_style as fs


ROOT = Path(__file__).resolve().parents[1]
TAB_DIR = ROOT / "data" / "final_summary_tables"
FIG_DIR = ROOT / "figures"
OUT_DATA = ROOT / "data" / "source_inputs" / "external_processed" / "fcv_postdca_raw_recompute" / "out_data"

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

MEASURES = [
    ("FCS", "FCS"),
    ("ProfileCorrDistFCV", "Profile corr-distance"),
    ("ObservedNetTE", "Observed NetTE"),
    ("NeighborNetTE", "Neighbor NetTE"),
]
HEATMAP_MEASURES = [
    ("EdgeStdFCV", "FCV"),
    *MEASURES,
]

GEOMETRIC_MEASURES = ["ProfileCorrDistFCV"]
PLOT_MEASURES_TABLE = TAB_DIR / "highpass_ce_zf_plot_measures_node_summary.csv"
PLOT_RECORDING_MEASURES_TABLE = TAB_DIR / "highpass_ce_zf_plot_measures_recording_node.csv"
PLOT_ZSCORE_MEASURES_TABLE = TAB_DIR / "highpass_ce_zf_plot_measures_recording_zscore_node_summary.csv"
OBSERVED_NETTE_RECORDING_TABLE = TAB_DIR / "observed_nette_no_p_recording_level.csv"
OBSERVED_NETTE_SUMMARY_TABLE = TAB_DIR / "observed_nette_no_p_node_summary.csv"
DROSOPHILA_FCS_RECORDING_TABLE = (
    ROOT
    / "figure15_drosophila_final/final_results/Branson999_full_FC_FCV/"
    / "Branson999_full_ROI_5measure_recording_level_w15_step5_hp030_FINAL.csv"
)
CLASS_TABLE = TAB_DIR / "current_plot_functional_classes.csv"

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

ANATOMY_COLORS = {
    "olf./chemo": fs.MAIN_COLORS["olfactory"],
    "sensory neurons": fs.MAIN_COLORS["neutral_light"],
    "interneurons": fs.MAIN_COLORS["neutral_light"],
    "motor / command": fs.MAIN_COLORS["neutral_light"],
    "state-modulatory": fs.MAIN_COLORS["neutral_light"],
    "olfactory neuropil": fs.MAIN_COLORS["olfactory"],
    "optic neuropil": fs.MAIN_COLORS["neutral_light"],
    "mushroom body": fs.MAIN_COLORS["associative"],
    "central complex": fs.MAIN_COLORS["neutral_light"],
    "protocerebrum": fs.MAIN_COLORS["neutral_light"],
    "premotor interface": fs.MAIN_COLORS["neutral_light"],
    "telencephalon": fs.ZEBRAFISH_DIVISION_COLORS["Tel"],
    "diencephalon": fs.ZEBRAFISH_DIVISION_COLORS["Di"],
    "mesencephalon": fs.ZEBRAFISH_DIVISION_COLORS["Mes"],
    "hindbrain / motor": fs.ZEBRAFISH_DIVISION_COLORS["Hind"],
    "state system": fs.MAIN_COLORS["neutral_light"],
    "other": fs.MAIN_COLORS["neutral_light"],
}

ANATOMY_ORDER = {
    "C. elegans": ["olf./chemo", "sensory neurons", "interneurons", "motor / command", "state-modulatory"],
    "Drosophila": ["olfactory neuropil", "optic neuropil", "mushroom body", "central complex", "protocerebrum", "premotor interface"],
    "Zebrafish": ["telencephalon", "diencephalon", "mesencephalon", "hindbrain / motor", "state system", "other"],
}
ANATOMY_MARKERS = {
    "olf./chemo": "o",
    "sensory neurons": "^",
    "interneurons": "D",
    "motor / command": "s",
    "state-modulatory": "X",
    "olfactory neuropil": "o",
    "optic neuropil": "^",
    "mushroom body": "P",
    "central complex": "X",
    "protocerebrum": "D",
    "premotor interface": "s",
    "telencephalon": "P",
    "diencephalon": "D",
    "mesencephalon": "^",
    "hindbrain / motor": "s",
    "state system": "X",
    "other": "o",
}

SPECIES = [
    {
        "species": "C. elegans",
        "level": "neuron",
        "path": OUT_DATA
        / "celegans/geometric_fcv/celegans_geometric_fcv_neuron_summary_w60_step15_minrec3_spontaneous.csv",
        "label_col": "neuron",
        "fcv_col": "EdgeStdFCV",
        "measure_prefix": "",
        "color": fs.SPECIES_COLORS["C. elegans"],
    },
    {
        "species": "Drosophila",
        "level": "side key",
        "path": TAB_DIR / "drosophila_roi_geometric_fcv_sidekey41_mean_summary.csv",
        "label_col": "side_key",
        "fcv_col": "EdgeStdFCV",
        "measure_prefix": "",
        "color": fs.SPECIES_COLORS["Drosophila"],
    },
    {
        "species": "Zebrafish",
        "level": "root area",
        "path": TAB_DIR / "zebrafish_raw_cluster_geometric_fcv_region_summary.csv",
        "label_col": "root_area_name",
        "fcv_col": "EdgeStdFCV",
        "measure_prefix": "",
        "color": fs.SPECIES_COLORS["Zebrafish"],
    },
]


def strip_side(name: str) -> str:
    raw = str(name)
    if raw.startswith("r") and raw not in {"RIR", "RID", "RIH", "RIS", "RIVL", "RIVR"}:
        return raw[1:]
    if raw.endswith("_L") or raw.endswith("_R"):
        return raw[:-2]
    return raw


def anatomy_group(species: str, node: str, fine_class: str) -> str:
    n = strip_side(str(node)).upper()
    fine = str(fine_class).lower()
    if species == "C. elegans":
        if any(k in fine for k in ["chemosensory", "olfactory"]):
            return "olf./chemo"
        if any(k in fine for k in ["sensory", "thermo", "mechanosensory"]):
            return "sensory neurons"
        if any(k in fine for k in ["motor", "command"]):
            return "motor / command"
        if "state" in fine:
            return "state-modulatory"
        return "interneurons"
    if species == "Drosophila":
        if n in {"AL", "LH"}:
            return "olfactory neuropil"
        if n in {"ME", "LO", "AME", "AOTU"}:
            return "optic neuropil"
        if n.startswith("MB"):
            return "mushroom body"
        if n in {"PB", "FB", "EB", "NO"}:
            return "central complex"
        if n in {"LAL", "VES", "SPS"}:
            return "premotor interface"
        return "protocerebrum"
    if species == "Zebrafish":
        if n in {"P", "SP", "OB", "OG", "OE"}:
            return "telencephalon"
        if n in {"TH", "PT", "PRT", "T"}:
            return "diencephalon"
        if n in {"TEO", "TL", "TS"}:
            return "mesencephalon"
        if n in {"HB", "HC", "HI", "HR", "IPN", "RA", "PO"}:
            return "state system"
        if n in {"MON", "CB", "MOS1", "MOS2", "MOS3", "MOS4", "MOS5", "IO", "ARF", "IMRF", "PRF", "TG", "VR", "NX"}:
            return "hindbrain / motor"
        return "other"
    return "other"


def load_species_table(spec: dict) -> pd.DataFrame:
    if PLOT_RECORDING_MEASURES_TABLE.exists():
        df = compute_recording_zscore_node_summary()
        df = df[df["species"].eq(spec["species"])]
        keep = ["species", "node", "EdgeStdFCV", *[name for name, _ in MEASURES]]
        df = df[[col for col in keep if col in df.columns]]
        return add_functional_classes(df)

    if PLOT_MEASURES_TABLE.exists():
        df = pd.read_csv(PLOT_MEASURES_TABLE)
        df = augment_observed_nette_summary(df)
        df = df[df["species"].eq(spec["species"])]
        keep = ["species", "node", "EdgeStdFCV", *[name for name, _ in MEASURES]]
        df = df[[col for col in keep if col in df.columns]]
        return add_functional_classes(df)

    df = pd.read_csv(spec["path"]).rename(columns={spec["label_col"]: "node"})
    df.insert(0, "species", spec["species"])
    df = augment_observed_nette_summary(df)
    keep = ["species", "node", spec["fcv_col"], *GEOMETRIC_MEASURES]
    df = df[[col for col in keep if col in df.columns]].rename(columns={spec["fcv_col"]: "EdgeStdFCV"})

    for path, cols in [
        (TAB_DIR / "fcs_autocorr_halflife_node_summary.csv", ["FCSHalfLifeWindows"]),
        (TAB_DIR / "fcs_fc_autocorr_node_summary.csv", ["FCS"]),
        (TAB_DIR / "participation_flex_node_summary.csv", ["ParticipationFlex"]),
    ]:
        extra = pd.read_csv(path)
        extra = extra[extra["species"].eq(spec["species"])][["species", "node", *cols]]
        df = df.merge(extra, on=["species", "node"], how="left")
    return add_functional_classes(df)


def add_functional_classes(df: pd.DataFrame) -> pd.DataFrame:
    classes = pd.read_csv(CLASS_TABLE)
    if "class_short" not in classes.columns:
        classes["class_short"] = classes["class_label"]
    classes = classes[["species", "node", "class_order", "class_label", "class_short", "fine_class"]]
    out = df.merge(classes, on=["species", "node"], how="left")
    out["class_order"] = out["class_order"].fillna(-1).astype(int)
    out["class_label"] = out["class_label"].fillna("unclassified")
    out["class_short"] = out["class_short"].fillna("unclassified")
    out["fine_class"] = out["fine_class"].fillna("unclassified")
    out["anatomy_group"] = [
        anatomy_group(species, node, fine)
        for species, node, fine in zip(out["species"], out["node"], out["fine_class"], strict=False)
    ]
    return out


def clean_xy(df: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    keep = [x_col, y_col, "class_order", "class_label"]
    out = df[keep].replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])
    return out.rename(columns={x_col: "x", y_col: "y"})


def add_fit(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str) -> None:
    if len(x) < 3:
        return
    slope, intercept, _, _, _ = stats.linregress(x, y)
    xs = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), 100)
    yhat = slope * xs + intercept
    fit = slope * x + intercept
    resid = y - fit
    dof = len(x) - 2
    if dof > 0:
        s_err = np.sqrt(np.sum(resid**2) / dof)
        x_mean = np.mean(x)
        ssx = np.sum((x - x_mean) ** 2)
        if ssx > 0:
            tcrit = stats.t.ppf(0.975, dof)
            ci = tcrit * s_err * np.sqrt(1 / len(x) + (xs - x_mean) ** 2 / ssx)
            ax.fill_between(xs, yhat - ci, yhat + ci, color="#222222", alpha=0.12, lw=0, zorder=1)
    ax.plot(xs, yhat, color="#222222", lw=0.8, zorder=3)


def corr_text(x: np.ndarray, y: np.ndarray) -> str:
    if len(x) < 3:
        return "n < 3"
    if np.nanstd(x) <= 1e-12 or np.nanstd(y) <= 1e-12:
        return f"n={len(x)}\nr=nan\nP=nan"
    pr, pp = stats.pearsonr(x, y)
    if pp < 1e-3:
        p_text = "P<0.001"
    else:
        p_text = f"P={pp:.3f}"
    return f"n={len(x)}\nr={pr:.2f}\n{p_text}"


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return out
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd <= 1e-12:
        out[finite] = 0.0
        return out
    out[finite] = (values[finite] - np.nanmean(values)) / sd
    return out


def augment_drosophila_recording_fcs(rec: pd.DataFrame) -> pd.DataFrame:
    fly_existing = rec[rec["species"].eq("Drosophila")]
    if "FCS" in fly_existing.columns and fly_existing["FCS"].notna().any():
        return rec
    if not DROSOPHILA_FCS_RECORDING_TABLE.exists():
        return rec
    fcs = pd.read_csv(DROSOPHILA_FCS_RECORDING_TABLE)
    fcs = (
        fcs.groupby(["recording_id", "atlas_region"], as_index=False)
        .agg(FCS=("FCS", "mean"))
        .rename(columns={"atlas_region": "node"})
    )
    fcs.insert(0, "species", "Drosophila")
    other = rec[~rec["species"].eq("Drosophila")].copy()
    fly = rec[rec["species"].eq("Drosophila")].drop(columns=["FCS"], errors="ignore").merge(
        fcs,
        on=["species", "recording_id", "node"],
        how="left",
    )
    return pd.concat([other, fly], ignore_index=True)


def augment_observed_nette_recording(rec: pd.DataFrame) -> pd.DataFrame:
    if not OBSERVED_NETTE_RECORDING_TABLE.exists():
        return rec
    nette = pd.read_csv(OBSERVED_NETTE_RECORDING_TABLE)
    nette = nette[["species", "recording_id", "node", "NetTE", "NeighborNetTE"]].rename(columns={"NetTE": "ObservedNetTE"})
    rec = rec.drop(columns=["ObservedNetTE", "NeighborNetTE"], errors="ignore")
    return rec.merge(nette, on=["species", "recording_id", "node"], how="left")


def augment_observed_nette_summary(df: pd.DataFrame) -> pd.DataFrame:
    if not OBSERVED_NETTE_SUMMARY_TABLE.exists():
        return df
    nette = pd.read_csv(OBSERVED_NETTE_SUMMARY_TABLE)
    nette = nette[["species", "node", "NetTE", "NeighborNetTE"]].rename(columns={"NetTE": "ObservedNetTE"})
    df = df.drop(columns=["ObservedNetTE", "NeighborNetTE"], errors="ignore")
    return df.merge(nette, on=["species", "node"], how="left")


def compute_recording_zscore_node_summary() -> pd.DataFrame:
    """Z-score each plotted measure within recording, then average by node."""
    measure_cols = [measure for measure, _ in HEATMAP_MEASURES]
    rec = pd.read_csv(PLOT_RECORDING_MEASURES_TABLE).replace([np.inf, -np.inf], np.nan)
    rec = augment_drosophila_recording_fcs(rec)
    rec = augment_observed_nette_recording(rec)
    zframes = []
    for _, group in rec.groupby(["species", "recording_id"], sort=False):
        out = group.copy()
        for col in measure_cols:
            if col in out.columns:
                out[col] = zscore(out[col].to_numpy(float))
        zframes.append(out)
    zrec = pd.concat(zframes, ignore_index=True)
    agg_cols = [col for col in measure_cols if col in zrec.columns]
    summary = (
        zrec.groupby(["species", "node"], as_index=False)
        .agg(
            **{col: (col, "mean") for col in agg_cols},
            n_recordings=("recording_id", "nunique"),
            level=("level", lambda s: s.dropna().iloc[0] if len(s.dropna()) else np.nan),
            window_config=("window_config", lambda s: s.dropna().iloc[0] if len(s.dropna()) else np.nan),
            root_area_id=("root_area_id", lambda s: s.dropna().iloc[0] if len(s.dropna()) else np.nan),
        )
    )
    summary["normalization"] = "within_recording_zscore_then_node_mean"
    summary.to_csv(PLOT_ZSCORE_MEASURES_TABLE, index=False)
    return summary


def load_species_recording_table(spec: dict) -> pd.DataFrame:
    """Return within-recording z-scored recording x node rows for box/strip panels."""
    measure_cols = [measure for measure, _ in HEATMAP_MEASURES]
    rec = pd.read_csv(PLOT_RECORDING_MEASURES_TABLE).replace([np.inf, -np.inf], np.nan)
    rec = augment_drosophila_recording_fcs(rec)
    rec = augment_observed_nette_recording(rec)
    rec = rec[rec["species"].eq(spec["species"])].copy()
    keep = ["species", "recording_id", "node", *measure_cols]
    rec = rec[[col for col in keep if col in rec.columns]]
    zframes = []
    for _, group in rec.groupby(["species", "recording_id"], sort=False):
        out = group.copy()
        for col in measure_cols:
            if col in out.columns:
                out[col] = zscore(out[col].to_numpy(float))
        zframes.append(out)
    if not zframes:
        return rec
    zrec = pd.concat(zframes, ignore_index=True)
    return add_functional_classes(zrec)


def draw_anatomy_panel(ax: plt.Axes, df: pd.DataFrame, species: str, add_colorbar: bool = False) -> pd.DataFrame:
    order = ANATOMY_ORDER[species]
    ax.set_axis_off()
    tmp = df.copy()
    measure_cols = [measure for measure, _ in HEATMAP_MEASURES]
    tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna(subset=measure_cols, how="all")
    tmp["anatomy_order"] = tmp["anatomy_group"].map({name: idx for idx, name in enumerate(order)}).fillna(len(order)).astype(int)
    tmp = tmp.sort_values(["anatomy_order", "anatomy_group", "node"]).reset_index(drop=True)
    matrix = np.vstack([zscore(tmp[measure].to_numpy(float)) for measure, _ in HEATMAP_MEASURES])
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)

    heat_ax = ax.inset_axes([0.0, 0.02, 1.0, 0.56], transform=ax.transAxes)
    bar_ax = ax.inset_axes([0.0, 0.62, 1.0, 0.055], transform=ax.transAxes)
    dend_ax = ax.inset_axes([0.0, 0.70, 1.0, 0.29], transform=ax.transAxes)
    legend_ax = ax.inset_axes([0.48, 0.835, 0.51, 0.145], transform=ax.transAxes)

    cluster_boundaries: list[float] = []
    if len(tmp) >= 3:
        z = linkage(matrix.T, method="ward")
        dend_info = dendrogram(
            z,
            ax=dend_ax,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#333333",
        )
        leaf_order = np.asarray(dend_info["leaves"], dtype=int)
        cluster_labels = fcluster(z, t=min(3, len(tmp)), criterion="maxclust")
        ordered_clusters = cluster_labels[leaf_order]
        cluster_boundaries = [
            idx - 0.5
            for idx in range(1, len(ordered_clusters))
            if ordered_clusters[idx] != ordered_clusters[idx - 1]
        ]
        max_height = max(max(coords) for coords in dend_info["dcoord"])
        dend_ax.set_xlim(0, len(tmp) * 10)
        dend_ax.set_ylim(0, max_height * 1.05 if max_height > 0 else 1)
        dend_ax.axis("off")
        tmp = tmp.iloc[leaf_order].reset_index(drop=True)
        tmp["heatmap_cluster_raw"] = cluster_labels[leaf_order]
        tmp["heatmap_cluster"] = np.r_[1, 1 + np.cumsum(ordered_clusters[1:] != ordered_clusters[:-1])]
        matrix = matrix[:, leaf_order]
    else:
        tmp["heatmap_cluster_raw"] = np.arange(len(tmp)) + 1
        tmp["heatmap_cluster"] = np.arange(len(tmp)) + 1
        dend_ax.axis("off")
    tmp["heatmap_order"] = np.arange(1, len(tmp) + 1)

    vmax = max(1.0, float(np.nanpercentile(np.abs(matrix), 98)))
    im = heat_ax.imshow(matrix, aspect="auto", cmap="PiYG_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    heat_ax.set_yticks(np.arange(len(HEATMAP_MEASURES)))
    heat_ax.set_yticklabels([label for _, label in HEATMAP_MEASURES], fontsize=5.6)
    heat_ax.yaxis.tick_right()
    heat_ax.tick_params(axis="y", labelright=True, labelleft=False)
    heat_ax.set_xticks(np.arange(len(tmp)))
    heat_ax.set_xticklabels(tmp["node"].astype(str).tolist(), rotation=90, ha="center", va="top")
    label_size = 2.4 if len(tmp) > 90 else 3.0 if len(tmp) > 55 else 4.0
    for tick, group in zip(heat_ax.get_xticklabels(), tmp["anatomy_group"], strict=False):
        tick.set_fontsize(label_size)
        tick.set_color(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
    dend_ax.set_title("Clustered\nmeasure profiles", fontsize=7, pad=1)
    heat_ax.tick_params(axis="x", length=0, pad=1.0)
    heat_ax.tick_params(axis="y", length=0, pad=1.0)
    for y in np.arange(-0.5, len(HEATMAP_MEASURES), 1):
        heat_ax.axhline(y=y, color="white", linewidth=0.35)
    for spine in heat_ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("black")

    color_indices = []
    group_to_idx = {name: idx for idx, name in enumerate(order)}
    for group in tmp["anatomy_group"]:
        color_indices.append(group_to_idx.get(group, len(order)))
    color_list = [ANATOMY_COLORS[name] for name in order] + [ANATOMY_COLORS["other"]]
    bar_ax.imshow(np.asarray(color_indices)[None, :], aspect="auto", cmap=ListedColormap(color_list), vmin=0, vmax=len(color_list) - 1)
    bar_ax.set_xlim(-0.5, len(tmp) - 0.5)
    bar_ax.set_axis_off()

    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.set_axis_off()
    legend_ax.patch.set_facecolor("white")
    legend_ax.patch.set_alpha(0.82)
    used_groups = [group for group in order if group in set(tmp["anatomy_group"])]
    if any(group not in order for group in tmp["anatomy_group"]):
        used_groups.append("other")
    n_groups = max(1, len(used_groups))
    n_cols = 2 if n_groups > 3 else n_groups
    n_rows = int(np.ceil(n_groups / n_cols))
    for idx, group in enumerate(used_groups):
        row = idx // n_cols
        col = idx % n_cols
        x0 = col / n_cols + 0.02
        y0 = 0.78 - row * (0.72 / max(1, n_rows - 1))
        legend_ax.add_patch(
            plt.Rectangle(
                (x0, y0 - 0.12),
                0.046,
                0.22,
                color=ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]),
                ec="#222222",
                lw=0.25,
                linewidth=0,
            )
        )
        legend_ax.text(
            x0 + 0.058,
            y0 - 0.02,
            group,
            ha="left",
            va="center",
            fontsize=4.5 if n_groups > 5 else 4.9,
            color="#222222",
            clip_on=False,
        )

    for boundary in cluster_boundaries:
        heat_ax.axvline(boundary, color="black", lw=0.8, zorder=5)
        bar_ax.axvline(boundary, color="black", lw=0.8, zorder=5)

    if add_colorbar:
        cax = ax.inset_axes([1.04, 0.02, 0.035, 0.56], transform=ax.transAxes)
        cbar = ax.figure.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=5, length=2, width=0.5)
        cbar.set_label("z-score", fontsize=5.5, labelpad=2)
    return tmp[["species", "node", "anatomy_group", "heatmap_order", "heatmap_cluster", "heatmap_cluster_raw"]]


def _p_text(p: float) -> str:
    if not np.isfinite(p):
        return "n.s."
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


def add_pairwise_sig_bars(
    ax: plt.Axes,
    values_by_group: list[np.ndarray],
    alpha: float = 0.05,
) -> list[dict[str, object]]:
    pairs = list(itertools.combinations(range(len(values_by_group)), 2))
    raw_p = []
    for i, j in pairs:
        a = values_by_group[i]
        b = values_by_group[j]
        if len(a) < 2 or len(b) < 2:
            raw_p.append(np.nan)
            continue
        try:
            raw_p.append(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except ValueError:
            raw_p.append(np.nan)
    adj_p = holm_adjust(raw_p)
    sig = [
        (pair, float(p), float(raw))
        for pair, p, raw in zip(pairs, adj_p, raw_p, strict=False)
        if np.isfinite(p) and p < alpha
    ]
    sig = sorted(sig, key=lambda item: (item[0][1] - item[0][0], item[1]))
    if not sig:
        return [
            {"group_i": i, "group_j": j, "mannwhitney_p": raw, "holm_p": adj}
            for (i, j), raw, adj in zip(pairs, raw_p, adj_p, strict=False)
        ]

    y_min, y_max = ax.get_ylim()
    yr = y_max - y_min if y_max > y_min else 1.0
    step = yr * 0.050
    bar_h = step * 0.22
    for lvl, ((i, j), p, _) in enumerate(sig):
        y = y_max + yr * 0.035 + lvl * step
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*"
        ax.plot([i, i, j, j], [y, y + bar_h, y + bar_h, y], lw=0.55, c="#333333", clip_on=False)
        ax.text((i + j) / 2, y + bar_h, star, ha="center", va="bottom", fontsize=5.6, clip_on=False)
    ax.set_ylim(y_min, y_max)
    return [
        {"group_i": i, "group_j": j, "mannwhitney_p": raw, "holm_p": adj}
        for (i, j), raw, adj in zip(pairs, raw_p, adj_p, strict=False)
    ]


def _plot_group_boxstrip(
    ax: plt.Axes,
    df: pd.DataFrame,
    species: str,
    measure: str,
    label: str,
    show_ylabel: bool,
) -> list[dict[str, object]]:
    order = [group for group in ANATOMY_ORDER[species] if group in set(df["anatomy_group"])]
    if any(group not in ANATOMY_ORDER[species] for group in df["anatomy_group"]):
        order.append("other")

    rows: list[dict[str, object]] = []
    values_by_group = []
    rng_seed = sum(ord(ch) for ch in f"{species}-{measure}") % (2**32)
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
            jitter = rng.normal(0, 0.070, size=len(vals))
            ax.scatter(
                np.full(len(vals), idx) + jitter,
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
        patch.set_facecolor(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
        patch.set_alpha(0.34)
        patch.set_edgecolor(ANATOMY_COLORS.get(group, ANATOMY_COLORS["other"]))
        patch.set_linewidth(1.05)
        patch.set_zorder(1)

    valid_groups = [vals for vals in values_by_group if len(vals) > 0]
    if len(valid_groups) >= 2:
        try:
            stat, pval = stats.kruskal(*valid_groups)
        except ValueError:
            stat, pval = np.nan, np.nan
    else:
        stat, pval = np.nan, np.nan
    ax.text(
        0.02,
        0.98,
        _p_text(pval),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.1,
        color="#222222",
    )
    rows.append(
        {
            "species": species,
            "measure": measure,
            "anatomy_group": "__global_kruskal__",
            "n": int(sum(len(vals) for vals in values_by_group)),
            "mean": np.nan,
            "median": np.nan,
            "std": np.nan,
            "kruskal_h": float(stat) if np.isfinite(stat) else np.nan,
            "kruskal_p": float(pval) if np.isfinite(pval) else np.nan,
        }
    )
    pair_rows = add_pairwise_sig_bars(ax, values_by_group, alpha=0.05)
    for pair_row in pair_rows:
        gi = int(pair_row["group_i"])
        gj = int(pair_row["group_j"])
        rows.append(
            {
                "species": species,
                "measure": measure,
                "anatomy_group": "__pairwise_mannwhitney_holm__",
                "group_i": order[gi],
                "group_j": order[gj],
                "n_i": int(len(values_by_group[gi])),
                "n_j": int(len(values_by_group[gj])),
                "mannwhitney_p": pair_row["mannwhitney_p"],
                "holm_p": pair_row["holm_p"],
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
        ax.set_ylabel("within-recording z", fontsize=7.0)
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="both", direction="out", pad=1.3, length=2.2, width=0.55)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    return rows


def draw_species_anatomy_boxgrid(fig: plt.Figure, slot, df: pd.DataFrame, spec: dict) -> tuple[plt.Axes, list[dict[str, object]]]:
    panel_ax = fig.add_subplot(slot)
    panel_ax.set_axis_off()
    panel_ax.text(
        -0.05,
        0.50,
        spec["species"],
        transform=panel_ax.transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold",
        color=spec["color"],
    )

    inner = gridspec.GridSpecFromSubplotSpec(
        1,
        len(HEATMAP_MEASURES),
        subplot_spec=slot,
        wspace=0.16,
    )
    rows: list[dict[str, object]] = []
    for idx, (measure, label) in enumerate(HEATMAP_MEASURES):
        ax = fig.add_subplot(inner[0, idx])
        rows.extend(_plot_group_boxstrip(ax, df, spec["species"], measure, label, show_ylabel=idx == 0))
    return panel_ax, rows


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15.6, 8.2), constrained_layout=True)
    outer = gridspec.GridSpec(
        len(SPECIES),
        2,
        figure=fig,
        width_ratios=[1.24, 2.76],
        wspace=0.04,
        hspace=0.18,
    )

    rows: list[dict[str, object]] = []
    cluster_rows = []
    heat_axes = []
    right_axes = []
    for row_idx, spec in enumerate(SPECIES):
        df = load_species_table(spec)
        rec_df = load_species_recording_table(spec)
        heat_ax = fig.add_subplot(outer[row_idx, 0])
        heat_axes.append(heat_ax)
        cluster_rows.append(draw_anatomy_panel(heat_ax, df, spec["species"], add_colorbar=row_idx == 0))
        heat_ax.text(
            -0.30,
            0.5,
            spec["species"],
            transform=heat_ax.transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
        panel_ax, panel_rows = draw_species_anatomy_boxgrid(fig, outer[row_idx, 1], rec_df, spec)
        right_axes.append(panel_ax)
        rows.extend(panel_rows)

    label_specs = [("A", heat_axes[0]), ("B", right_axes[0]), ("C", right_axes[1]), ("D", right_axes[2])]
    for label, ax in label_specs:
        ax.text(-0.06, 1.04, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")

    legend_handles = [
        Line2D([0], [0], color=ANATOMY_COLORS[group], lw=5, label=group)
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
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=8,
        frameon=False,
        fontsize=5.6,
        handletextpad=0.4,
        columnspacing=0.9,
    )

    png = FIG_DIR / "geometric_fcv_variants_vs_edge_fcv_all_species.png"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)

    corr_path = TAB_DIR / "geometric_fcv_variants_vs_edge_fcv_all_species_anatomy_group_stats.csv"
    pd.DataFrame(rows).to_csv(corr_path, index=False)
    cluster_path = TAB_DIR / "geometric_fcv_heatmap_measure_profile_clusters.csv"
    pd.concat(cluster_rows, ignore_index=True).to_csv(cluster_path, index=False)
    print(f"wrote {png}")
    print(f"wrote {corr_path}")
    print(f"wrote {cluster_path}")


if __name__ == "__main__":
    main()
