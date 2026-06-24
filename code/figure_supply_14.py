import os
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

import figure7_by_species_with_figure5_abc as fig7abc
import figure_style as fs

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

fs.set_paper_style()
plt.rcParams.update({
    "font.size": fs.AXIS_LABEL_FS_2COL,
    "axes.labelsize": fs.AXIS_LABEL_FS_2COL,
    "axes.titlesize": fs.AXIS_LABEL_FS_2COL,
    "xtick.labelsize": max(5, fs.TICK_FS_2COL - 2),
    "ytick.labelsize": fs.TICK_FS_2COL,
})

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(PROJECT_ROOT, "data")
IN_CSV = os.path.join(
    DATA,
    "final_summary_tables",
    "figure7_zebrafish_all_region_stimulus_values.csv",
)
OUT_PNG = os.path.join(PROJECT_ROOT, "output", "png", "figure_supply_14.png")
OUT_CSV = os.path.join(PROJECT_ROOT, "output", "stats", "figure_supply_14_stimulus_fcv_heatmap_values.csv")
FIG_PNG = os.path.join(PROJECT_ROOT, "figures", "figure_supply_14.png")
SUPP_PNG = os.path.join(PROJECT_ROOT, "outputs", "supplementary", "figure_supply_14.png")

DIVISION_ORDER = ["Tel", "Di", "Mes", "Hind"]
DIVISION_COLORS = {
    "Tel": fs.division_colors[2],
    "Di": fs.division_colors[1],
    "Mes": fs.division_colors[3],
    "Hind": fs.division_colors[0],
}
STIM_KEEP = [10, 11, 12]
STIM_LABELS = {
    10: "10: OMR forward",
    11: "11: OMR rightward",
    12: "12: OMR leftward",
}


def zscore(values):
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd == 0:
        return values * np.nan
    return (values - np.nanmean(values)) / sd


def load_heatmap_table():
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df = pd.read_csv(IN_CSV)
    df = df[df["StimulusIndex"].isin(STIM_KEEP)].copy()
    if "SubjectZFCV" in df.columns:
        df["SubjectStimulusZFCV"] = df["SubjectZFCV"]
    else:
        df["SubjectStimulusZFCV"] = df.groupby(["Subject", "StimulusIndex"])["FCV"].transform(zscore)
    region_meta = (
        df[["Region", "RegionID", "Division"]]
        .drop_duplicates("Region")
        .copy()
    )
    region_meta["Division"] = pd.Categorical(
        region_meta["Division"],
        categories=DIVISION_ORDER,
        ordered=True,
    )
    region_meta = region_meta.sort_values(["Division", "RegionID"]).reset_index(drop=True)
    regions = region_meta["Region"].tolist()
    stim_values = STIM_KEEP

    mean_df = (
        df.groupby(["StimulusIndex", "Region"], as_index=False)
        .agg(
            MeanStimulusFCV=("SubjectStimulusZFCV", "mean"),
            SEMStimulusFCV=("SubjectStimulusZFCV", lambda x: np.nanstd(x, ddof=1) / np.sqrt(np.isfinite(x).sum())),
            N=("SubjectStimulusZFCV", "count"),
        )
    )
    mean_df = mean_df.merge(region_meta, on="Region", how="left")
    mean_df = mean_df.sort_values(["StimulusIndex", "Division", "RegionID"])
    mean_df.to_csv(OUT_CSV, index=False)

    lookup = {
        (int(row.StimulusIndex), row.Region): row.MeanStimulusFCV
        for row in mean_df.itertuples(index=False)
    }
    matrix = np.full((len(stim_values), len(regions)), np.nan)
    for i, stim in enumerate(stim_values):
        for j, region in enumerate(regions):
            matrix[i, j] = lookup.get((stim, region), np.nan)
    return matrix, stim_values, region_meta


def add_division_bar(ax, region_meta, label_inside=False):
    divisions = region_meta["Division"].astype(str).tolist()
    div_to_idx = {div: i for i, div in enumerate(DIVISION_ORDER)}
    div_arr = np.array([[div_to_idx[d] for d in divisions]])
    div_cmap = ListedColormap([DIVISION_COLORS[d] for d in DIVISION_ORDER])
    ax.imshow(div_arr, aspect="auto", cmap=div_cmap, vmin=-0.5, vmax=len(DIVISION_ORDER) - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    boundaries = division_boundaries(divisions)
    for boundary in boundaries:
        ax.axvline(boundary, color="white", lw=0.8)

    x0 = 0
    for div in DIVISION_ORDER:
        n = sum(d == div for d in divisions)
        if n:
            y_pos = 0.0 if label_inside else -0.65
            ax.text(
                x0 + (n - 1) / 2,
                y_pos,
                div,
                ha="center",
                va="center" if label_inside else "bottom",
                fontsize=fs.TICK_FS_2COL,
                color="white" if label_inside else DIVISION_COLORS[div],
                fontweight="bold",
            )
        x0 += n


def division_boundaries(divisions):
    return [
        idx - 0.5
        for idx in range(1, len(divisions))
        if divisions[idx] != divisions[idx - 1]
    ]


def attach_bar_to_heatmap(ax_bar, ax_heat, gap=0.006, height=0.026):
    heat_pos = ax_heat.get_position()
    ax_bar.set_position([
        heat_pos.x0,
        heat_pos.y1 + gap,
        heat_pos.width,
        height,
    ])


def make_figure():
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    os.makedirs(os.path.dirname(FIG_PNG), exist_ok=True)
    os.makedirs(os.path.dirname(SUPP_PNG), exist_ok=True)
    matrix, stim_values, region_meta = load_heatmap_table()
    regions = region_meta["Region"].tolist()
    divisions = region_meta["Division"].astype(str).tolist()
    ce, zf_region_values = fig7abc.load_stimulus_tables()

    fig = plt.figure(figsize=(8.4, 7.1))
    gs = fig.add_gridspec(
        4,
        4,
        height_ratios=[1.08, 0.12, 0.94, 0.04],
        width_ratios=[1.0, 1.0, 1.0, 0.040],
        left=0.08,
        right=0.95,
        top=0.94,
        bottom=0.14,
        hspace=0.78,
        wspace=0.34,
    )
    scatter_axes = [fig.add_subplot(gs[0, col]) for col in range(3)]
    ax_bar = fig.add_subplot(gs[1, 0:3])
    ax = fig.add_subplot(gs[2, 0:3])
    cax = fig.add_subplot(gs[2, 3])

    fig7abc.scatter_with_functional_group(
        scatter_axes[0],
        ce,
        "SpontaneousFCV",
        "HeatFCV",
        "spontaneous FCV",
        "Heat FCV",
    )
    fig7abc.scatter_with_functional_group(
        scatter_axes[1],
        zf_region_values,
        "SpontaneousFCV",
        "StdStimulusFCVZ",
        "spontaneous FCV",
        "Std OMR FCV",
    )
    fig7abc.scatter_with_functional_group(
        scatter_axes[2],
        zf_region_values,
        "SpontaneousFCV",
        "MeanStimulusFCV",
        "spontaneous FCV",
        "Mean OMR FCV",
    )
    for scatter_ax, species in zip(scatter_axes, ["C. elegans", "Zebrafish", "Zebrafish"], strict=False):
        scatter_ax.set_title(
            species,
            fontsize=fig7abc.SCATTER_TITLE_FS,
            fontweight="bold",
            color=fig7abc.SPECIES_COLORS.get(species, "#111111"),
            pad=3,
        )
    fig7abc.style_scatter_axes(scatter_axes)

    add_division_bar(ax_bar, region_meta, label_inside=True)
    attach_bar_to_heatmap(ax_bar, ax, gap=0.008, height=0.028)
    vmax = np.nanpercentile(np.abs(matrix), 98)
    im = ax.imshow(
        matrix,
        aspect="auto",
        cmap="PiYG_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    for boundary in division_boundaries(divisions):
        ax.axvline(boundary, color="white", lw=0.8)

    ax.set_xlabel("Region")
    ax.set_ylabel("Stimulus")
    ax.set_xticks(np.arange(len(regions)))
    ax.set_xticklabels(
        regions,
        rotation=90,
        ha="center",
        fontsize=fs.TICK_FS_2COL - 2,
        fontweight="bold",
    )
    for tick_label, division in zip(ax.get_xticklabels(), divisions):
        tick_label.set_color(DIVISION_COLORS[division])
    ax.set_yticks(np.arange(len(stim_values)))
    ax.set_yticklabels([STIM_LABELS[stim] for stim in stim_values])

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Mean subject z-scored FCV")
    cbar.ax.tick_params(labelsize=fs.TICK_FS_2COL)

    ax.text(
        -0.065,
        1.12,
        "D",
        transform=ax.transAxes,
        fontsize=fs.PANEL_LABEL_FS_2COL,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(FIG_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(SUPP_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    make_figure()
    print(f"Saved {OUT_PNG}")
    print(f"Saved {FIG_PNG}")
    print(f"Saved {SUPP_PNG}")
    print(f"Saved {OUT_CSV}")


if __name__ == "__main__":
    main()
