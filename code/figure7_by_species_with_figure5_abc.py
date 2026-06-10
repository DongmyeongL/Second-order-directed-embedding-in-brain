"""Combine Figure 7 by-species synthesis with Figure 5 A-C panels.

The combined figure is drawn from the underlying tables rather than pasting
existing PNG panels.
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
import pickle
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
logging.getLogger("fontTools").setLevel(logging.WARNING)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy import stats

import figure_style as fs
from figure5_stimulus_evoked_fcv_modulation import (
    format_p,
    load_celegans_evoked,
    load_tables,
)
from figure7_structure_function_synthesis_effects import (
    BUILD_FIG,
    BUILD_TAB,
    DISPLAY_GROUP_ORDER,
    GROUP_COLORS,
    GROUP_LABELS,
    HIGHLIGHT_GROUPS,
    SPECIES,
    SPECIES_COLORS,
    load_or_compute_statistics,
    metric_tables,
    p_to_stars,
    significance_lookup,
)


ROOT = Path(__file__).resolve().parents[1]
OUT_FIG = BUILD_FIG / "figure7_by_species_with_figure5_abc.png"
OUT_EFFECTS = BUILD_TAB / "figure7_by_species_with_figure5_abc_effects.csv"
OUT_GROUP = BUILD_TAB / "figure7_by_species_with_figure5_abc_group_signature.csv"
OUT_GROUP_PERM = BUILD_TAB / "figure7_by_species_with_figure5_abc_group_permutation_tests.csv"
NODE_GROUP_TABLE = (
    ROOT
    / "stats_tables"
    / "figure7_comparative_fcv_fine_class_node_map.csv"
)
ZEBRAFISH_LEGACY_FCV_DCA_TABLE = (
    ROOT
    / "data"
    / "source_inputs"
    / "external_processed"
    / "fcv_postdca_raw_recompute"
    / "out_data"
    / "zebrafish"
    / "post_dca_rank1"
    / "zebrafish_rank1_post_dca_fcv_legacy_merged.csv"
)
ZEBRAFISH_RAW_FC_TEMPLATE = (
    ROOT
    / "data"
    / "raw_not_bundled"
    / "subject_{subject}_data_cellular_synapse_sc_100_data.pkl"
)
ZEBRAFISH_ALL_STIMULUS_CACHE = BUILD_TAB / "figure7_zebrafish_all_region_stimulus_values.csv"
DEFAULT_GROUP_COLOR = "#9A9A9A"
FIG_SIZE = (8.2, 9.20)
AXIS_FS = fs.AXIS_LABEL_FS_2COL
TICK_FS = fs.TICK_FS_2COL
PANEL_FS = fs.PANEL_LABEL_FS_2COL
STAR_FS = fs.STAR_FS_2COL
STAT_FS = fs.STAT_FS_2COL
GROUP_LABEL_FS = fs.TICK_FS_2COL - 2
SCATTER_TITLE_FS = AXIS_FS
SCATTER_TICK_FS = TICK_FS
SCATTER_LABEL_FS = AXIS_FS
SCATTER_STAT_FS = STAT_FS
SCATTER_ANNOTATION_FS = STAT_FS
SCATTER_POINT_SIZE = fs.MAIN_SCATTER_SIZE * 1.45
SIGNATURE_GRID_Y_SHIFT = -0.018
SIGNATURE_YLABEL_X = -0.16
FOOTNOTE_FS = STAT_FS
SUMMARY_METRICS = ["FCV", "Post-DCA", "Pre-DCA"]
METRIC_LABELS = {
    "FCV": "FCV",
    "Post-DCA": r"$\mathrm{DCA}_{\mathrm{post}}$",
    "Pre-DCA": r"$\mathrm{DCA}_{\mathrm{pre}}$",
}
SPECIES_MARKERS = {"C. elegans": "o", "Drosophila": "s", "Zebrafish": "^"}
ASSOC_STATE_ORDER = 5
STATE_MOD_ORDER = 6
ZEBRAFISH_SUBJECTS = range(12, 19)
ZEBRAFISH_STIM_KEEP = (10, 11, 12)
STIM_WINDOW_SIZE = 20
STIM_OVERLAP = 15
PLOT_GROUP_ORDER = [group for group in DISPLAY_GROUP_ORDER if group != STATE_MOD_ORDER]
PLOT_GROUP_LABELS = GROUP_LABELS.copy()
PLOT_GROUP_LABELS[ASSOC_STATE_ORDER] = "assoc./state"
METRIC_YLIMS = {
    "FCV": (-2.6, 2.35),
    "Post-DCA": (-3.8, 2.35),
    "Pre-DCA": (-2.2, 2.05),
}


PANEL_LABEL_X = -0.3
PANEL_LABEL_Y = 1.08

def add_panel_label(ax: plt.Axes, label: str) -> None:
    if label=='C':
        dx = PANEL_LABEL_X -0.05
    else:
        dx = PANEL_LABEL_X
        
        
    ax.text(
        dx,
        PANEL_LABEL_Y,
        label,
        transform=ax.transAxes,
        fontsize=PANEL_FS,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    
    
#def add_panel_label(ax: plt.Axes, label: str) -> None:
#    fs.add_main_panel_label(ax, label, fontsize=PANEL_FS)


def functional_group_mapping(species: str) -> pd.DataFrame:
    mapping = pd.read_csv(NODE_GROUP_TABLE)
    return (
        mapping.loc[
            mapping["species"].eq(species),
            ["node", "shared_fine_order", "shared_fine_label"],
        ]
        .drop_duplicates("node")
        .copy()
    )


def add_functional_group(df: pd.DataFrame, species: str, key_col: str) -> pd.DataFrame:
    mapping = functional_group_mapping(species)
    out = df.copy()
    out["_fg_key"] = out[key_col].astype(str)
    if species == "Zebrafish":
        direct = dict(zip(mapping["node"].astype(str), mapping["shared_fine_order"]))
        labels = dict(zip(mapping["node"].astype(str), mapping["shared_fine_label"]))
        orders = []
        names = []
        for key in out["_fg_key"]:
            stripped = key[1:] if len(key) > 1 and key[0] in {"l", "r"} else key
            match = key if key in direct else stripped
            orders.append(direct.get(match, np.nan))
            names.append(labels.get(match, np.nan))
        out["shared_fine_order"] = orders
        out["shared_fine_label"] = names
        return out.drop(columns=["_fg_key"])
    out = out.merge(
        mapping.rename(columns={"node": "_fg_key"}),
        on="_fg_key",
        how="left",
    )
    return out.drop(columns=["_fg_key"])


def collapse_state_mod_into_assoc(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "shared_fine_order" not in out.columns:
        return out
    order = pd.to_numeric(out["shared_fine_order"], errors="coerce")
    state_mask = order.eq(STATE_MOD_ORDER)
    out.loc[state_mask, "shared_fine_order"] = ASSOC_STATE_ORDER
    if "shared_fine_label" in out.columns:
        out.loc[state_mask, "shared_fine_label"] = PLOT_GROUP_LABELS[ASSOC_STATE_ORDER]
        out.loc[order.eq(ASSOC_STATE_ORDER), "shared_fine_label"] = PLOT_GROUP_LABELS[ASSOC_STATE_ORDER]
    return out


def collapse_metric_tables(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    return {metric: collapse_state_mod_into_assoc(df) for metric, df in tables.items()}


def use_bilateral_zebrafish_ob_fcv(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Replace zebrafish olfactory FCV with subject means from OB and rOB."""
    out = {metric: df.copy() for metric, df in tables.items()}
    if "FCV" not in out or not ZEBRAFISH_LEGACY_FCV_DCA_TABLE.exists():
        return out

    zf = pd.read_csv(ZEBRAFISH_LEGACY_FCV_DCA_TABLE)
    zf = zf.loc[zf["Region"].isin(["OB", "rOB"]), ["Subject", "Region", "ZFCV"]].dropna(subset=["ZFCV"])
    if zf.empty:
        return out
    bilateral = (
        zf.groupby("Subject", as_index=False)
        .agg(value=("ZFCV", "mean"), n_nodes=("Region", "nunique"))
        .query("n_nodes >= 1")
        .copy()
    )
    bilateral["species"] = "Zebrafish"
    bilateral["recording_id"] = "subject_" + bilateral["Subject"].astype(int).astype(str)
    bilateral["shared_fine_order"] = 0
    bilateral["shared_fine_label"] = "olf./chemo"
    bilateral["metric"] = "FCV"
    bilateral = bilateral[
        ["species", "recording_id", "shared_fine_order", "shared_fine_label", "value", "metric"]
    ]

    fcv = out["FCV"]
    keep = ~(fcv["species"].eq("Zebrafish") & pd.to_numeric(fcv["shared_fine_order"], errors="coerce").eq(0))
    out["FCV"] = pd.concat([fcv.loc[keep].copy(), bilateral], ignore_index=True, sort=False)
    return out


def relabel_collapsed_group_tables(*frames: pd.DataFrame) -> tuple[pd.DataFrame, ...]:
    out_frames = []
    for frame in frames:
        out = frame.copy()
        if {"shared_fine_order", "shared_fine_label"}.issubset(out.columns):
            order = pd.to_numeric(out["shared_fine_order"], errors="coerce")
            out.loc[order.eq(ASSOC_STATE_ORDER), "shared_fine_label"] = PLOT_GROUP_LABELS[ASSOC_STATE_ORDER]
        out_frames.append(out)
    return tuple(out_frames)


def _zscore0(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd <= 1e-12:
        return arr * np.nan
    return (arr - np.nanmean(arr)) / sd


def _zscore_rows(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    mean = np.nanmean(matrix, axis=1, keepdims=True)
    sd = np.nanstd(matrix, axis=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd <= 1e-12)] = np.nan
    return (matrix - mean) / sd


def _region_mean_traces(raw: dict, data_key: str, region_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    traces = np.asarray(raw[data_key], dtype=np.float32)
    neuron_region_id = np.asarray(raw["neuron_region_id"], dtype=int)
    root_area = np.asarray(raw["root_area"], dtype=int)
    neuron_root = root_area[neuron_region_id]

    region_traces = np.full((len(region_ids), traces.shape[1]), np.nan, dtype=np.float32)
    neuron_counts = np.zeros(len(region_ids), dtype=int)
    for pos, region_id in enumerate(region_ids):
        idx = np.flatnonzero(neuron_root == int(region_id))
        neuron_counts[pos] = idx.size
        if idx.size:
            region_traces[pos] = np.nanmean(traces[idx], axis=0)
    return region_traces, neuron_counts


def _fcv_by_region(region_traces: np.ndarray) -> tuple[np.ndarray, int]:
    if region_traces.shape[1] < STIM_WINDOW_SIZE:
        return np.full(region_traces.shape[0], np.nan), 0
    traces = _zscore_rows(region_traces)
    step = STIM_WINDOW_SIZE - STIM_OVERLAP
    corr_list = []
    for start in range(0, traces.shape[1] - STIM_WINDOW_SIZE + 1, step):
        corr_list.append(np.corrcoef(traces[:, start : start + STIM_WINDOW_SIZE]))
    if not corr_list:
        return np.full(region_traces.shape[0], np.nan), 0
    corr_stack = np.asarray(corr_list)
    std_fc = np.nanstd(corr_stack, axis=0)
    np.fill_diagonal(std_fc, np.nan)
    return np.nanmean(std_fc, axis=1), corr_stack.shape[0]


def _all_region_stimulus_rows(region_values: pd.DataFrame) -> pd.DataFrame:
    """Return all-region OMR FCV rows, computing and caching them on first use."""
    if ZEBRAFISH_ALL_STIMULUS_CACHE.exists():
        cached = pd.read_csv(ZEBRAFISH_ALL_STIMULUS_CACHE)
        if (
            "normalization" in cached.columns
            and cached["normalization"].astype(str).eq("subject_stimulus_zscore").all()
        ):
            return cached

    meta = (
        region_values[["Region", "RegionID", "Division"]]
        .drop_duplicates()
        .loc[lambda df: ~df["Region"].astype(str).isin(["lOB", "OB", "rOB"])]
        .copy()
    )
    ob_meta = pd.DataFrame(
        [
            {"Region": "OB", "RegionID": 20, "Division": "Tel"},
            {"Region": "rOB", "RegionID": 56, "Division": "Tel"},
        ]
    )
    meta = pd.concat([ob_meta, meta], ignore_index=True, sort=False)
    meta["Division"] = pd.Categorical(meta["Division"], ["Tel", "Di", "Mes", "Hind"], ordered=True)
    meta = meta.sort_values(["Division", "RegionID"]).reset_index(drop=True)

    region_ids = pd.to_numeric(meta["RegionID"], errors="coerce").to_numpy(int)
    rows = []
    for subject in ZEBRAFISH_SUBJECTS:
        raw_path = Path(str(ZEBRAFISH_RAW_FC_TEMPLATE).format(subject=subject))
        with raw_path.open("rb") as f:
            raw = pickle.load(f)
        stim_traces, neuron_counts = _region_mean_traces(raw, "stim_data", region_ids)
        stim_array = np.asarray(raw["stim_array"]).ravel()
        for stim_index in ZEBRAFISH_STIM_KEEP:
            time_idx = np.flatnonzero(stim_array == stim_index)
            if time_idx.size >= STIM_WINDOW_SIZE:
                fcv_values, n_windows = _fcv_by_region(stim_traces[:, time_idx])
            else:
                fcv_values = np.full(len(meta), np.nan)
                n_windows = 0
            for pos, reg in enumerate(meta.itertuples(index=False)):
                rows.append(
                    {
                        "Subject": int(subject),
                        "StimulusIndex": int(stim_index),
                        "Region": str(reg.Region),
                        "RegionID": int(reg.RegionID),
                        "Division": str(reg.Division),
                        "FCV": float(fcv_values[pos]),
                        "NeuronCount": int(neuron_counts[pos]),
                        "NTimepoints": int(time_idx.size),
                        "NWindows": int(n_windows),
                    }
                )
        del raw, stim_traces
        gc.collect()

    stim = pd.DataFrame(rows)
    # For stimulus-condition comparisons, normalize regions within each
    # subject and stimulus condition so each OMR condition is centered.
    stim["SubjectZFCV"] = stim.groupby(["Subject", "StimulusIndex"])["FCV"].transform(_zscore0)
    stim["normalization"] = "subject_stimulus_zscore"
    ZEBRAFISH_ALL_STIMULUS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    stim.to_csv(ZEBRAFISH_ALL_STIMULUS_CACHE, index=False)
    return stim


def add_recomputed_zebrafish_stimulus_summary(region_values: pd.DataFrame) -> pd.DataFrame:
    """Replace zebrafish B/C stimulus readouts with all-region recomputed OMR FCV."""
    stim = _all_region_stimulus_rows(region_values)
    legacy = pd.read_csv(ZEBRAFISH_LEGACY_FCV_DCA_TABLE)

    def region_key(region: object) -> str:
        name = str(region)
        if name == "lOB":
            return "OB"
        if name.startswith("l") and len(name) > 1 and name[1].isupper():
            return name[1:]
        return name

    stim_summary = (
        stim.groupby("Region", as_index=False)
        .agg(
            MeanStimulusFCV=("SubjectZFCV", "mean"),
            SEMStimulusFCV=("SubjectZFCV", lambda x: np.nanstd(x, ddof=1) / np.sqrt(np.isfinite(x).sum())),
            StdStimulusFCV=("SubjectZFCV", "std"),
            RegionID=("RegionID", "first"),
            Division=("Division", "first"),
        )
    )
    stim_summary["RegionKey"] = stim_summary["Region"].map(region_key)
    legacy = legacy.copy()
    legacy["RegionKey"] = legacy["Region"].map(region_key)
    spont_summary = (
        legacy.groupby("RegionKey", as_index=False)
        .agg(
            SpontaneousFCV=("ZFCV", "mean"),
            SpontaneousFCVSEM=("ZFCV", lambda x: np.nanstd(x, ddof=1) / np.sqrt(np.isfinite(x).sum())),
            RawPostDCA=("LegacyPostDCA", "mean"),
            RawPreDCA=("LegacyPreDCA", "mean"),
        )
    )

    out = region_values.copy()
    out["Region"] = out["Region"].astype(str).replace({"lOB": "OB"})
    out = out.loc[~out["Region"].eq("rOB")].copy()
    out["RegionKey"] = out["Region"].map(region_key)

    missing_regions = sorted(set(stim_summary["Region"]) - set(out["Region"]))
    rows = []
    for region in missing_regions:
        out_row = {col: np.nan for col in region_values.columns}
        out_row["Region"] = region
        out_row["RegionKey"] = region_key(region)
        rows.append(out_row)
    if rows:
        missing_df = pd.DataFrame(rows).dropna(axis=1, how="all")
        out = pd.concat([out, missing_df], ignore_index=True, sort=False)

    update_cols = ["MeanStimulusFCV", "SEMStimulusFCV", "StdStimulusFCV", "RegionID", "Division"]
    out = out.drop(columns=[col for col in update_cols if col in out.columns], errors="ignore")
    out = out.merge(stim_summary[["Region", *update_cols]], on="Region", how="left")

    for row in spont_summary.itertuples(index=False):
        mask = out["RegionKey"].astype(str).eq(str(row.RegionKey))
        for col in ["SpontaneousFCV", "SpontaneousFCVSEM", "RawPostDCA", "RawPreDCA"]:
            if col not in out.columns:
                out[col] = np.nan
            out.loc[mask, col] = getattr(row, col)

    if "MeanStimulusFCVZ" in out.columns:
        out["MeanStimulusFCVZ"] = _zscore0(out["MeanStimulusFCV"])
    if "StdStimulusFCVZ" in out.columns:
        out["StdStimulusFCVZ"] = _zscore0(out["StdStimulusFCV"])
    out["Division"] = pd.Categorical(out["Division"], ["Tel", "Di", "Mes", "Hind"], ordered=True)
    return out.drop(columns=["RegionKey"], errors="ignore").sort_values(["Division", "RegionID"]).reset_index(drop=True)


def scatter_with_functional_group(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    xlabel: str,
    ylabel: str,
) -> None:
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(float)
    group = pd.to_numeric(df.get("shared_fine_order", np.nan), errors="coerce").to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    colors = [
        GROUP_COLORS.get(int(g), DEFAULT_GROUP_COLOR) if np.isfinite(g) else DEFAULT_GROUP_COLOR
        for g in group
    ]
    ax.scatter(
        x[ok],
        y[ok],
        s=SCATTER_POINT_SIZE,
        c=np.asarray(colors, dtype=object)[ok],
        alpha=0.82,
        linewidth=fs.MAIN_SCATTER_EDGE_LW,
        edgecolors="#222222",
    )
    if ok.sum() >= 4:
        r, p = stats.pearsonr(x[ok], y[ok])
        m, b = np.polyfit(x[ok], y[ok], 1)
        xx = np.linspace(np.nanpercentile(x[ok], 2), np.nanpercentile(x[ok], 98), 80)
        ax.plot(xx, m * xx + b, color="#333333", lw=fs.MAIN_LINE_W)
        ax.text(
            0.04,
            0.95,
            f"r={r:.2f}\np={format_p(p)}",
            transform=ax.transAxes,
            fontsize=SCATTER_STAT_FS,
            ha="left",
            va="top",
        )
    ax.axhline(0, color=fs.MAIN_REFERENCE_COLOR, lw=fs.MAIN_GRID_LW, ls=":")
    ax.axvline(0, color=fs.MAIN_REFERENCE_COLOR, lw=fs.MAIN_GRID_LW, ls=":")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fs.style_main_axis(ax, tick_fs=TICK_FS)


def load_stimulus_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    _, zf_region_values, _ = load_tables()
    zf_region_values = add_recomputed_zebrafish_stimulus_summary(zf_region_values)
    zf_region_values = collapse_state_mod_into_assoc(add_functional_group(zf_region_values, "Zebrafish", "Region"))
    ce = collapse_state_mod_into_assoc(add_functional_group(load_celegans_evoked(), "C. elegans", "neuron"))
    return ce, zf_region_values


def style_scatter_axes(axes: list[plt.Axes]) -> None:
    for ax, label in zip(axes, ["A", "B", "C"], strict=False):
        ax.tick_params(axis="both", labelsize=SCATTER_TICK_FS, width=0.55, length=2.2)
        ax.xaxis.label.set_size(SCATTER_LABEL_FS)
        ax.yaxis.label.set_size(SCATTER_LABEL_FS)
        add_panel_label(ax, label)


def shift_axes(axes: list[plt.Axes], dx: float = 0.0, dy: float = 0.0) -> None:
    for ax in axes:
        pos = ax.get_position()
        ax.set_position([pos.x0 + dx, pos.y0 + dy, pos.width, pos.height])


def group_tick_color(group_order: int) -> str:
    return GROUP_COLORS.get(int(group_order), "#030303") if group_order in HIGHLIGHT_GROUPS else "#030303"


def zscore_within(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values)
    sd = np.nanstd(values, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return np.full(values.shape, np.nan, dtype=float)
    return (values - mean) / sd


def add_violin_box(
    ax: plt.Axes,
    values: np.ndarray,
    position: float,
    color: str,
    marker: str,
    mean_size: float = 64,
) -> None:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return
    if len(values) >= 2:
        parts = ax.violinplot(
            [values],
            positions=[position],
            widths=0.22,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.42)
            body.set_edgecolor("none")
            body.set_zorder(2)

        q1, med, q3 = np.percentile(values, [25, 50, 75])
        iqr = q3 - q1
        wlo = max(np.min(values), q1 - 1.5 * iqr)
        whi = min(np.max(values), q3 + 1.5 * iqr)
        ax.plot([position, position], [wlo, whi], color="#222222", linewidth=0.60, zorder=4)
        ax.add_patch(
            plt.Rectangle(
                (position - 0.045, q1),
                0.09,
                max(q3 - q1, 1e-6),
                facecolor="white",
                edgecolor="#222222",
                linewidth=0.55,
                zorder=5,
            )
        )
    else:
        med = float(values[0])

    mean_value = float(np.nanmean(values))

    ax.scatter(
        [position],
        [mean_value],
        s=mean_size,
        marker=marker,
        facecolor=color,
        edgecolor="#222222",
        linewidth=0.95,
        alpha=0.98,
        zorder=7,
    )


def draw_metric_species_grid(
    axes: list[list[plt.Axes]],
    tables: dict,
    group_perm: pd.DataFrame,
) -> None:
    sig = significance_lookup(group_perm)
    species_group_orders: dict[str, list[int]] = {}
    for species in SPECIES:
        groups = sorted(
            {
                int(group)
                for metric in SUMMARY_METRICS
                for group in (
                    tables[metric]
                    .loc[tables[metric]["species"].eq(species), "shared_fine_order"]
                    .dropna()
                    .astype(int)
                    .unique()
                )
            }
        )
        species_group_orders[species] = [group for group in PLOT_GROUP_ORDER if group in groups]

    for row_idx, metric in enumerate(SUMMARY_METRICS):
        for col_idx, species in enumerate(SPECIES):
            ax = axes[row_idx][col_idx]
            data = (
                tables[metric]
                .loc[tables[metric]["species"].eq(species)]
                .replace([np.inf, -np.inf], np.nan)
                .dropna(subset=["value", "shared_fine_order"])
                .copy()
            )
            data["shared_fine_order"] = data["shared_fine_order"].astype(int)
            data["plot_value"] = zscore_within(data["value"].to_numpy(float))

            orders = species_group_orders[species]
            stat_labels: list[tuple[float, str]] = []
            trend_x: list[float] = []
            trend_y: list[float] = []
            for idx, group_order in enumerate(orders):
                vals = (
                    data.loc[data["shared_fine_order"].eq(group_order), "plot_value"]
                    .dropna()
                    .to_numpy(float)
                )
                if len(vals) == 0:
                    continue
                mean_value = float(np.nanmean(vals))
                trend_x.append(float(idx))
                trend_y.append(mean_value)
                add_violin_box(
                    ax,
                    vals,
                    idx,
                    GROUP_COLORS.get(int(group_order), "#BBBBBB"),
                    SPECIES_MARKERS.get(species, "o"),
                    mean_size=72 if group_order in {0, ASSOC_STATE_ORDER} else 58,
                )
                q = sig.get((species, int(group_order), metric))
                if q is not None:
                    stat_labels.append((idx, p_to_stars(q)))

            if len(trend_x) >= 2:
                ax.plot(
                    trend_x,
                    trend_y,
                    color="#242424",
                    lw=1.25,
                    alpha=0.82,
                    zorder=6,
                    solid_capstyle="round",
                )

            ymin, ymax = METRIC_YLIMS.get(metric, (-3.0, 3.0))
            ax.set_ylim(ymin, ymax)
            stat_y = ymax - 0.12 * (ymax - ymin)
            for x, label in stat_labels:
                ax.text(
                    x,
                    stat_y,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=STAR_FS,
                    fontweight="heavy",
                    color="#111111",
                    clip_on=False,
                    zorder=8,
                )

            ax.axhline(0, color="#777777", lw=0.55, ls="--", alpha=0.55, zorder=0)
            for boundary in np.arange(0.5, len(orders) - 0.5, 1.0):
                ax.axvline(boundary, color="#D2D2D2", lw=0.45, ls=(0, (2.0, 2.0)), alpha=0.75, zorder=-1)
            ax.set_xlim(-0.55, len(orders) - 0.45)
            if row_idx == 0:
                ax.set_title(
                    species,
                    fontsize=AXIS_FS,
                    fontweight="bold",
                    color=SPECIES_COLORS.get(species, "#111111"),
                    pad=3,
                )
            else:
                ax.set_title("")
            ax.set_ylabel(METRIC_LABELS[metric], fontsize=AXIS_FS)
            ax.yaxis.set_label_coords(SIGNATURE_YLABEL_X, 0.5)
            ax.set_xticks(np.arange(len(orders), dtype=float))
            if row_idx == len(SUMMARY_METRICS) - 1:
                ax.set_xticklabels(
                    [PLOT_GROUP_LABELS.get(int(group), str(group)) for group in orders],
                    fontsize=GROUP_LABEL_FS - 1,
                    rotation=32,
                    ha="right",
                )
                for tick, group_order in zip(ax.get_xticklabels(), orders, strict=False):
                    tick.set_color(group_tick_color(int(group_order)))
                    tick.set_fontweight("bold" if group_order in HIGHLIGHT_GROUPS else "normal")
            else:
                ax.set_xticklabels([])
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            #ax.grid(axis="y", color="#ECECEC", lw=0.42, zorder=-5)
            ax.tick_params(axis="both", direction="out", pad=1.1, length=2.0, width=0.55, labelsize=TICK_FS - 1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_color("#BFBFBF")
            ax.spines["bottom"].set_linewidth(0.55)
            ax.spines["left"].set_linewidth(0.60)


def style_signature_axes(axes: list[plt.Axes]) -> None:
    for ax in axes:
        ax.title.set_fontsize(AXIS_FS)
        ax.xaxis.label.set_size(AXIS_FS)
        ax.yaxis.label.set_size(AXIS_FS)
        ax.tick_params(axis="y", labelsize=TICK_FS)
        ax.tick_params(axis="x", labelsize=GROUP_LABEL_FS)
        for text in ax.texts:
            value = text.get_text()
            if value.startswith("p="):
                text.set_fontsize(STAT_FS)
            elif "*" in value:
                text.set_fontsize(STAR_FS)
            else:
                text.set_fontsize(STAT_FS)
        legend = ax.get_legend()
        if legend is not None:
            for text in legend.get_texts():
                text.set_fontsize(TICK_FS)


def main() -> None:
    parser = argparse.ArgumentParser(description="Draw combined Figure 7 plus stimulus FCV panels.")
    parser.add_argument(
        "--recompute-stats",
        action="store_true",
        help="Recompute Figure 7 permutation/bootstrap statistics instead of reading cached CSV files.",
    )
    args = parser.parse_args()

    fs.apply_main_figure_style()

    tables = use_bilateral_zebrafish_ob_fcv(collapse_metric_tables(metric_tables()))
    tables, effects, group_signature, group_perm = load_or_compute_statistics(
        tables,
        recompute=args.recompute_stats,
        effects_path=OUT_EFFECTS,
        group_path=OUT_GROUP,
        group_perm_path=OUT_GROUP_PERM,
    )
    effects, group_signature, group_perm = relabel_collapsed_group_tables(effects, group_signature, group_perm)
    effects.to_csv(OUT_EFFECTS, index=False)
    group_signature.to_csv(OUT_GROUP, index=False)
    group_perm.to_csv(OUT_GROUP_PERM, index=False)

    ce, zf_region_values = load_stimulus_tables()

    fig = plt.figure(figsize=FIG_SIZE)
    gs = fig.add_gridspec(
        4,
        3,
        height_ratios=[1.18, 0.82, 0.82, 0.92],
        left=0.055,
        right=0.985,
        bottom=0.125,
        top=0.955,
        wspace=0.34,
        hspace=0.44,
    )
    scatter_axes = [fig.add_subplot(gs[0, col]) for col in range(3)]
    signature_axes_grid = [
        [fig.add_subplot(gs[row, col]) for col in range(3)]
        for row in range(1, 4)
    ]
    signature_axes = [ax for row in signature_axes_grid for ax in row]
    shift_axes(signature_axes, dy=SIGNATURE_GRID_Y_SHIFT)

    scatter_with_functional_group(
        scatter_axes[0],
        ce,
        "SpontaneousFCV",
        "HeatFCV",
        "spontaneous FCV",
        "Heat FCV",
    )
    scatter_with_functional_group(
        scatter_axes[1],
        zf_region_values,
        "SpontaneousFCV",
        "StdStimulusFCVZ",
        "spontaneous FCV",
        "Std OMR FCV",
    )
    scatter_with_functional_group(
        scatter_axes[2],
        zf_region_values,
        "SpontaneousFCV",
        "MeanStimulusFCV",
        "spontaneous FCV",
        "Mean OMR FCV",
    )

    scatter_species = ["C. elegans", "Zebrafish", "Zebrafish"]
    for ax, species in zip(scatter_axes, scatter_species, strict=False):
        ax.set_title(
            species,
            fontsize=SCATTER_TITLE_FS,
            fontweight="bold",
            color=SPECIES_COLORS.get(species, "#111111"),
            pad=3,
        )

    style_scatter_axes(scatter_axes)

    draw_metric_species_grid(signature_axes_grid, tables, group_perm)
    style_signature_axes(signature_axes)
    for ax, label in zip([row[0] for row in signature_axes_grid], ["D", "E", "F"], strict=False):
        add_panel_label(ax, label)

    #fig.text(
    #    0.52,
    #    0.030,
    #    "* pHolm<0.05, ** pHolm<0.01, *** pHolm<0.001",
    #    ha="center",
    #    va="center",
    #    fontsize=FOOTNOTE_FS,
    #)
    fig.savefig(OUT_FIG, dpi=600, bbox_inches="tight", transparent=False)
    plt.close(fig)
    print(f"wrote {OUT_FIG}")
    print(f"wrote {OUT_EFFECTS}")
    print(f"wrote {OUT_GROUP}")
    print(f"wrote {OUT_GROUP_PERM}")


if __name__ == "__main__":
    main()
