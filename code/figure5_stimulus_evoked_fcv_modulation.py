#!/usr/bin/env python3
"""Build Figure 5 draft: zebrafish stimulus-evoked FCV hierarchy."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from PIL import Image
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT
FIG_DIR = OUT / "figures"
TAB_DIR = OUT / "data" / "final_summary_tables"
QC_DIR = OUT / "qc"

STIM_KEEP = [10, 11, 12]
STIM_LABELS = {10: "OMR forward", 11: "OMR right", 12: "OMR left"}
DIVISION_ORDER = ["Tel", "Di", "Mes", "Hind"]
DIVISION_COLORS = {
    "Tel": "#C44E52",
    "Di": "#8172B2",
    "Mes": "#55A868",
    "Hind": "#4C72B0",
}


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd <= 1e-12:
        return arr * np.nan
    return (arr - np.nanmean(arr)) / sd


def format_p(p: float) -> str:
    if not np.isfinite(p):
        return "nan"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="left", va="bottom")


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stim = pd.read_csv(TAB_DIR / "figure7_zebrafish_all_region_stimulus_values.csv")
    stim = stim[stim["StimulusIndex"].isin(STIM_KEEP)].copy()
    if "SubjectZFCV" not in stim.columns:
        stim["SubjectZFCV"] = stim.groupby(["Subject", "StimulusIndex"])["FCV"].transform(zscore)
    if "SubjectStimulusZFCV" not in stim.columns:
        stim["SubjectStimulusZFCV"] = stim["SubjectZFCV"]

    region_values = (
        stim[["Region", "RegionID", "Division"]]
        .drop_duplicates("Region")
        .copy()
    )
    stim_summary = (
        stim.groupby("Region", as_index=False)
        .agg(
            MeanStimulusFCV=("SubjectZFCV", "mean"),
            StdStimulusFCVZ=("SubjectZFCV", "std"),
        )
    )
    region_values = region_values.merge(stim_summary, on="Region", how="left")
    region_values["Division"] = pd.Categorical(region_values["Division"], DIVISION_ORDER, ordered=True)
    region_values = region_values.sort_values(["Division", "RegionID"]).reset_index(drop=True)

    heat = (
        stim.groupby(["StimulusIndex", "Region"], as_index=False)
        .agg(
            MeanStimulusFCV=("SubjectZFCV", "mean"),
            SEMStimulusFCV=("SubjectZFCV", lambda x: np.nanstd(x, ddof=1) / np.sqrt(np.isfinite(x).sum())),
            N=("SubjectZFCV", "count"),
        )
        .merge(region_values[["Region", "RegionID", "Division"]], on="Region", how="inner")
    )
    heat["Division"] = pd.Categorical(heat["Division"], DIVISION_ORDER, ordered=True)
    heat = heat.sort_values(["StimulusIndex", "Division", "RegionID"]).reset_index(drop=True)
    return stim, region_values, heat


def compute_stats(region_values: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("SpontaneousFCV", "MeanStimulusFCV", "spontaneous FCV", "mean OMR FCV"),
        ("SpontaneousFCV", "StdStimulusFCVZ", "spontaneous FCV", "OMR FCV modulation"),
        ("PostDCA", "MeanStimulusFCV", "Post-DCA", "mean OMR FCV"),
        ("PostDCA", "StdStimulusFCVZ", "Post-DCA", "OMR FCV modulation"),
    ]
    rows = []
    for x_col, y_col, x_label, y_label in pairs:
        x = pd.to_numeric(region_values[x_col], errors="coerce").to_numpy(float)
        y = pd.to_numeric(region_values[y_col], errors="coerce").to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() >= 4:
            pr, pp = stats.pearsonr(x[ok], y[ok])
            sr, sp = stats.spearmanr(x[ok], y[ok])
        else:
            pr = pp = sr = sp = np.nan
        rows.append(
            {
                "x": x_col,
                "y": y_col,
                "x_label": x_label,
                "y_label": y_label,
                "n": int(ok.sum()),
                "pearson_r": float(pr),
                "pearson_p": float(pp),
                "spearman_rho": float(sr),
                "spearman_p": float(sp),
            }
        )
    out = pd.DataFrame(rows)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(TAB_DIR / "figure5_stimulus_fcv_correlations.csv", index=False)
    return out


def load_celegans_evoked() -> pd.DataFrame:
    spont = pd.read_csv(
        ROOT
        / "data/source_inputs/external_processed/fcv_postdca_raw_recompute/out_data/celegans/geometric_fcv/"
        / "celegans_geometric_fcv_post_dca_merged_w60_step15_minrec3_spontaneous.csv"
    )
    heat = pd.read_csv(
        ROOT
        / "data/source_inputs/external_processed/fcv_postdca_raw_recompute/out_data/celegans/geometric_fcv/"
        / "celegans_geometric_fcv_post_dca_merged_w60_step15_minrec3_heat.csv"
    )
    keep = [
        "neuron",
        "cell_class",
        "PostDCA",
        "EdgeStdFCV_z",
        "ProfileCorrDistFCV_z",
        "ProfileTransitionCorrDistFCV_z",
    ]
    out = spont[keep].rename(
        columns={
            "EdgeStdFCV_z": "SpontaneousFCV",
            "ProfileCorrDistFCV_z": "SpontaneousPatternFCV",
            "ProfileTransitionCorrDistFCV_z": "SpontaneousTransitionFCV",
        }
    ).merge(
        heat[keep].rename(
            columns={
                "EdgeStdFCV_z": "HeatFCV",
                "ProfileCorrDistFCV_z": "HeatPatternFCV",
                "ProfileTransitionCorrDistFCV_z": "HeatTransitionFCV",
                "PostDCA": "PostDCA_heat_join",
                "cell_class": "cell_class_heat_join",
            }
        ),
        on="neuron",
        how="inner",
    )
    out = out.drop(columns=["PostDCA_heat_join", "cell_class_heat_join"], errors="ignore")
    return out


def load_drosophila_evoked() -> pd.DataFrame:
    return pd.read_csv(TAB_DIR / "positive_negative_prepost_dca_drosophila_roi.csv").rename(
        columns={"FCV_z_final": "EvokedFCV"}
    )


def load_drosophila_evoked_vs_original(fly_roi: pd.DataFrame) -> pd.DataFrame:
    original = pd.read_csv(TAB_DIR / "figure1_cross_species_fcv_node_table.csv")
    original = original[original["species"].eq("Drosophila")][["node", "node_class", "fcv", "fcv_z"]].rename(
        columns={"node": "side_key", "node_class": "big_group", "fcv": "SpontaneousFCV", "fcv_z": "SpontaneousFCV_z"}
    )
    evoked = (
        fly_roi.groupby(["side_key", "big_group"], as_index=False)
        .agg(
            EvokedFCV=("mean_FCV_raw_recording_mean", "mean"),
            EvokedFCV_z=("EvokedFCV", "mean"),
            n_rois=("roi", "count"),
        )
    )
    out = original.merge(evoked, on=["side_key", "big_group"], how="inner")
    return out


def append_cross_species_stats(zf_stats: pd.DataFrame, ce: pd.DataFrame, fly: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def add(species: str, x_col: str, y_col: str, x_label: str, y_label: str, df: pd.DataFrame) -> None:
        x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(float)
        y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() >= 4:
            pr, pp = stats.pearsonr(x[ok], y[ok])
            sr, sp = stats.spearmanr(x[ok], y[ok])
        else:
            pr = pp = sr = sp = np.nan
        rows.append(
            {
                "species": species,
                "x": x_col,
                "y": y_col,
                "x_label": x_label,
                "y_label": y_label,
                "n": int(ok.sum()),
                "pearson_r": float(pr),
                "pearson_p": float(pp),
                "spearman_rho": float(sr),
                "spearman_p": float(sp),
            }
        )

    zf_region = pd.read_csv(TAB_DIR / "figure5_stimulus_fcv_region_values.csv")
    add("Zebrafish", "SpontaneousFCV", "MeanStimulusFCV", "spontaneous FCV", "OMR FCV", zf_region)
    add("Zebrafish", "PostDCA", "MeanStimulusFCV", "Post-DCA", "OMR FCV", zf_region)
    add("C. elegans", "SpontaneousFCV", "HeatFCV", "spontaneous FCV", "heat FCV", ce)
    add("C. elegans", "PostDCA", "HeatFCV", "Post-DCA", "heat FCV", ce)
    add("Drosophila", "SpontaneousFCV_z", "EvokedFCV_z", "spontaneous FCV", "Branson evoked FCV", fly)
    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / "figure5_cross_species_stimulus_fcv_correlations.csv", index=False)
    return out


def plot_heatmap(ax: plt.Axes, heat: pd.DataFrame, region_values: pd.DataFrame, fig: plt.Figure) -> None:
    region_meta = region_values[["Region", "RegionID", "Division"]].drop_duplicates().sort_values(["Division", "RegionID"])
    regions = region_meta["Region"].tolist()
    divisions = region_meta["Division"].astype(str).tolist()
    lookup = {(int(r.StimulusIndex), r.Region): r.MeanStimulusFCV for r in heat.itertuples(index=False)}
    matrix = np.full((len(STIM_KEEP), len(regions)), np.nan)
    for i, stim in enumerate(STIM_KEEP):
        for j, region in enumerate(regions):
            matrix[i, j] = lookup.get((stim, region), np.nan)
    vmax = np.nanpercentile(np.abs(matrix), 98)
    im = ax.imshow(matrix, aspect="auto", cmap="PiYG_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    boundaries = [i - 0.5 for i in range(1, len(divisions)) if divisions[i] != divisions[i - 1]]
    for boundary in boundaries:
        ax.axvline(boundary, color="white", lw=0.8)
    ax.set_xticks(np.arange(len(regions)))
    ax.set_xticklabels(regions, rotation=90, fontsize=4.8)
    for label, div in zip(ax.get_xticklabels(), divisions, strict=False):
        label.set_color(DIVISION_COLORS.get(div, "#333333"))
    ax.set_yticks(np.arange(len(STIM_KEEP)))
    ax.set_yticklabels([STIM_LABELS[s] for s in STIM_KEEP], fontsize=6.2)
    ax.set_xlabel("Region")
    ax.set_ylabel("Stimulus")
    ax.set_title("OMR-evoked FCV across regions", fontsize=8, pad=3)
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.01)
    cbar.set_label("mean subject-z FCV", fontsize=6)
    cbar.ax.tick_params(labelsize=5)


def scatter_with_fit(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    xlabel: str,
    ylabel: str,
    annotate: bool = True,
    color: str = "#777777",
) -> None:
    if "Division" in df.columns:
        colors = [DIVISION_COLORS.get(str(d), "#777777") for d in df["Division"]]
    elif "cell_class" in df.columns:
        class_colors = {"sensory": "#3A6EA5", "interneuron": "#2A9D8F", "motor": "#E76F51", "other": "#777777"}
        colors = [class_colors.get(str(c), color) for c in df["cell_class"]]
    else:
        colors = [color] * len(df)
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(float)
    y = pd.to_numeric(df[y_col], errors="coerce").to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[ok], y[ok], s=28, c=np.asarray(colors, dtype=object)[ok], alpha=0.78, linewidth=0)
    if ok.sum() >= 4:
        r, p = stats.pearsonr(x[ok], y[ok])
        m, b = np.polyfit(x[ok], y[ok], 1)
        xx = np.linspace(np.nanpercentile(x[ok], 2), np.nanpercentile(x[ok], 98), 80)
        ax.plot(xx, m * xx + b, color="#333333", lw=1.2)
        ax.text(0.04, 0.95, f"r={r:.2f}, p={format_p(p)}, n={ok.sum()}", transform=ax.transAxes, fontsize=6.6, ha="left", va="top")
    if annotate:
        top = df.sort_values(y_col, ascending=False).head(4)
        for row in top.itertuples(index=False):
            label = getattr(row, "Region", None)
            if label is None:
                label = getattr(row, "neuron", None)
            if label is None:
                label = getattr(row, "label", None)
            if label is None:
                label = getattr(row, "side_key", None)
            if label is None:
                label = getattr(row, "roi", "")
            ax.text(getattr(row, x_col), getattr(row, y_col), str(label), fontsize=5.5, ha="left", va="bottom")
    ax.axhline(0, color="#AAAAAA", lw=0.6, ls=":")
    ax.axvline(0, color="#AAAAAA", lw=0.6, ls=":")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_division_summary(ax: plt.Axes, region_values: pd.DataFrame) -> None:
    data = [region_values.loc[region_values["Division"].astype(str) == d, "MeanStimulusFCV"].dropna().to_numpy(float) for d in DIVISION_ORDER]
    bp = ax.boxplot(data, positions=np.arange(len(DIVISION_ORDER)), widths=0.55, patch_artist=True, showfliers=False)
    for patch, div in zip(bp["boxes"], DIVISION_ORDER, strict=False):
        patch.set_facecolor(DIVISION_COLORS[div])
        patch.set_alpha(0.45)
        patch.set_edgecolor("#333333")
    rng = np.random.default_rng(3)
    for i, (div, vals) in enumerate(zip(DIVISION_ORDER, data, strict=False)):
        x = i + rng.normal(0, 0.045, size=len(vals))
        ax.scatter(x, vals, s=16, color=DIVISION_COLORS[div], alpha=0.75, linewidth=0)
    ax.axhline(0, color="#AAAAAA", lw=0.6, ls=":")
    ax.set_xticks(np.arange(len(DIVISION_ORDER)))
    ax.set_xticklabels(DIVISION_ORDER)
    ax.set_ylabel("mean OMR FCV")
    ax.set_title("Stimulus FCV by division", fontsize=8, pad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_condition_profiles(ax: plt.Axes, heat: pd.DataFrame, region_values: pd.DataFrame) -> None:
    top_regions = region_values.sort_values("SpontaneousFCV", ascending=False).head(4)["Region"].tolist()
    low_regions = region_values.sort_values("SpontaneousFCV", ascending=True).head(4)["Region"].tolist()
    xs = np.arange(len(STIM_KEEP))
    for regions, label, color, alpha in [
        (top_regions, "high spontaneous FCV", "#C44E52", 0.95),
        (low_regions, "low spontaneous FCV", "#4C72B0", 0.75),
    ]:
        profiles = []
        for region in regions:
            vals = []
            for stim in STIM_KEEP:
                sub = heat[(heat["Region"] == region) & (heat["StimulusIndex"] == stim)]
                vals.append(float(sub["MeanStimulusFCV"].iloc[0]) if len(sub) else np.nan)
            profiles.append(vals)
            ax.plot(xs, vals, color=color, alpha=0.25, lw=0.9)
        mean = np.nanmean(np.asarray(profiles, dtype=float), axis=0)
        ax.plot(xs, mean, color=color, lw=2.0, label=label)
    ax.axhline(0, color="#AAAAAA", lw=0.6, ls=":")
    ax.set_xticks(xs)
    ax.set_xticklabels([STIM_LABELS[s].replace("OMR ", "") for s in STIM_KEEP], rotation=20, ha="right")
    ax.set_ylabel("mean subject-z FCV")
    ax.set_title("OMR profile by intrinsic-FCV rank", fontsize=8, pad=3)
    ax.legend(frameon=False, fontsize=6, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_cross_species_bars(ax: plt.Axes, stats_df: pd.DataFrame) -> None:
    show = stats_df.copy()
    label_map = {
        ("Zebrafish", "SpontaneousFCV", "MeanStimulusFCV"): "ZF spont -> OMR",
        ("Zebrafish", "PostDCA", "MeanStimulusFCV"): "ZF Post-DCA -> OMR",
        ("C. elegans", "SpontaneousFCV", "HeatFCV"): "Ce spont -> heat",
        ("C. elegans", "PostDCA", "HeatFCV"): "Ce Post-DCA -> heat",
        ("Drosophila", "SpontaneousFCV_z", "EvokedFCV_z"): "Fly spont -> evoked",
    }
    show["label"] = [label_map.get((r.species, r.x, r.y), f"{r.species}: {r.x_label} -> {r.y_label}") for r in show.itertuples()]
    y = np.arange(len(show), dtype=float)
    colors = [
        {"Zebrafish": "#3A6EA5", "C. elegans": "#2A9D8F", "Drosophila": "#E76F51"}.get(s, "#777777")
        for s in show["species"]
    ]
    ax.barh(y, show["pearson_r"], color=colors, alpha=0.88)
    ax.axvline(0, color="#777777", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(show["label"], fontsize=6.1)
    ax.set_xlabel("Pearson r")
    ax.set_xlim(-0.65, 0.72)
    ax.set_title("Stimulus/evoked FCV associations", fontsize=8, pad=3)
    for yy, p in zip(y, show["pearson_p"], strict=False):
        ax.text(0.50, yy, f"p={format_p(p)}", fontsize=5.6, va="center", ha="left")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_stimulus_scope(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    rows = [
        ("Zebrafish", "OMR 10-12", "condition-resolved"),
        ("C. elegans", "heat phase", "phase-resolved"),
        ("Drosophila", "Branson responses", "evoked epoch only"),
    ]
    ax.text(0.02, 0.95, "Stimulus definitions", fontsize=8, fontweight="bold", ha="left", va="top")
    y = 0.78
    for species, condition, scope in rows:
        color = {"Zebrafish": "#3A6EA5", "C. elegans": "#2A9D8F", "Drosophila": "#E76F51"}[species]
        ax.scatter(0.05, y, s=70, color=color)
        ax.text(0.12, y + 0.035, species, fontsize=7.2, fontweight="bold", ha="left", va="center")
        ax.text(0.12, y - 0.025, f"{condition}; {scope}", fontsize=6.5, ha="left", va="center")
        y -= 0.22
    ax.text(
        0.02,
        0.08,
        "Drosophila condition labels were not present locally;\n"
        "therefore this panel uses whole Branson\n"
        "stimulus-response FCV.",
        fontsize=6.0,
        ha="left",
        va="bottom",
    )


def plot_figure() -> Path:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7})
    _stim, region_values, heat = load_tables()
    zf_stats = compute_stats(region_values)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    region_values.to_csv(TAB_DIR / "figure5_stimulus_fcv_region_values.csv", index=False)
    heat.to_csv(TAB_DIR / "figure5_stimulus_fcv_heatmap_values.csv", index=False)
    ce = load_celegans_evoked()
    fly = load_drosophila_evoked()
    fly_region = load_drosophila_evoked_vs_original(fly)
    ce.to_csv(TAB_DIR / "figure5_celegans_heat_fcv_table.csv", index=False)
    fly.to_csv(TAB_DIR / "figure5_drosophila_branson_evoked_fcv_table.csv", index=False)
    fly_region.to_csv(TAB_DIR / "figure5_drosophila_spontaneous_evoked_fcv_table.csv", index=False)
    cross_stats = append_cross_species_stats(zf_stats, ce, fly_region)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes_arr = plt.subplots(2, 2, figsize=(7.4, 6.6), constrained_layout=True)
    axes = axes_arr.ravel().tolist()
    for ax, label in zip(axes, list("ABCD"), strict=False):
        panel_label(ax, label)

    scatter_with_fit(axes[0], region_values, "SpontaneousFCV", "MeanStimulusFCV", "spontaneous FCV", "mean OMR FCV")
    scatter_with_fit(axes[1], ce, "SpontaneousFCV", "HeatFCV", "spontaneous FCV", "heat FCV", color="#2A9D8F")
    scatter_with_fit(
        axes[2],
        fly_region,
        "SpontaneousFCV_z",
        "EvokedFCV_z",
        "spontaneous FCV",
        "evoked FCV",
        color="#E76F51",
    )
    scatter_with_fit(
        axes[3],
        region_values,
        "SpontaneousFCV",
        "StdStimulusFCVZ",
        "spontaneous FCV",
        "std mean OMR FCV",
    )

    fig.suptitle("Figure 5. Stimulus-evoked FCV preserves the intrinsic FCV hierarchy", fontsize=11, fontweight="bold")
    out = FIG_DIR / "figure5_stimulus_fcv_draft.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def run_qc(path: Path) -> dict:
    img = Image.open(path)
    arr = np.asarray(img)
    stats_df = pd.read_csv(TAB_DIR / "figure5_cross_species_stimulus_fcv_correlations.csv")
    qc = {
        "figure": str(path.relative_to(ROOT)),
        "figure_exists": path.exists(),
        "image_width_px": img.width,
        "image_height_px": img.height,
        "pixel_std": float(arr.std()),
        "nonblank": bool(arr.std() > 1.0),
        "correlation_rows": int(len(stats_df)),
        "min_correlation_n": int(stats_df["n"].min()),
        "species": sorted(stats_df["species"].unique().tolist()),
    }
    qc["pass"] = bool(qc["figure_exists"] and qc["nonblank"] and qc["correlation_rows"] >= 5 and len(qc["species"]) == 3)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    (QC_DIR / "figure5_stimulus_fcv_qc.json").write_text(json.dumps(qc, indent=2))
    return qc


def main() -> None:
    path = plot_figure()
    print(json.dumps(run_qc(path), indent=2))


if __name__ == "__main__":
    main()
