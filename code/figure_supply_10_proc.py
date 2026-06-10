import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import figure_style as fs
from figure_style import add_panel_label_fig, darw_region_bar


fs.set_paper_style()

MAIN_AXIS_FS = 11
MAIN_TICK_FS = 10
MAIN_PANEL_FS = 12
MAIN_AXIS_LW = 1.0
MAIN_ERROR_COLOR = "#2F3437"
MAIN_DIVISION_COLORS = {
    0: "#6FA8C9",  # Hind
    1: "#DDA15E",  # Di
    2: "#6FAF6B",  # Tel
    3: "#E76F51",  # Mes
}
MAIN_DIVISION_ORDER = [2, 1, 3, 0]

PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PROJECT_ROOT.parent
DATA = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "figures"
SUBJECT_IDS = range(12, 19)
BASE_NET = DATA / "region_community_io"

DIVISION_COLUMNS = ["Tel", "Di", "Mes", "Hind"]
DIVISION_ORDER = {"Tel": 0, "Di": 1, "Mes": 2, "Hind": 3}

TEL_REGIONS = {"lOB", "lSP", "lP", "lPO", "rOB", "rSP", "rP", "rPO"}
DI_REGIONS = {"lHb", "lTh", "lPT", "lPrT", "rHb", "rTh", "rPT", "rPrT"}
MES_REGIONS = {"lTeO", "lTL", "lTS", "rTeO", "rTL", "rTS"}
HIND_REGIONS = {
    "lCb", "lT", "laRF", "lMOS1", "lMOS2", "lMOS3", "lMOS4", "lMOS5",
    "limRF", "lpRF", "lMON", "lNX", "lRa",
    "rCb", "rT", "raRF", "rMOS1", "rMOS2", "rMOS3", "rMOS4", "rMOS5",
    "rimRF", "rpRF", "rMON", "rNX", "rRa",
}


def first_existing_path(*paths):
    for path in paths:
        candidate = Path(path)
        if candidate.exists():
            return candidate
        if not candidate.is_absolute():
            project_candidate = PROJECT_ROOT / candidate
            if project_candidate.exists():
                return project_candidate
            workspace_candidate = WORKSPACE_ROOT / candidate
            if workspace_candidate.exists():
                return workspace_candidate
    return Path(paths[0])


def region_to_division(region_name):
    if region_name in TEL_REGIONS:
        return "Tel"
    if region_name in DI_REGIONS:
        return "Di"
    if region_name in MES_REGIONS:
        return "Mes"
    if region_name in HIND_REGIONS:
        return "Hind"
    return None


def node_to_region_index(node):
    node = str(node)
    if node in fs.region:
        return fs.region.index(node)
    left_node = f"l{node}"
    if left_node in fs.region:
        return fs.region.index(left_node)
    return None


def display_region_name(region_name):
    region_name = str(region_name)
    if len(region_name) > 1 and region_name[0] == "l" and region_name[1].isupper():
        return region_name[1:]
    return region_name


def zebrafish_division(node):
    base = display_region_name(node).upper()
    if base in {"P", "SP", "OB", "OG", "OE", "PO"}:
        return "Tel"
    if base in {"HB", "HI", "HR", "TH", "PT", "PRT"}:
        return "Di"
    if base in {"TEO", "TL", "TS"}:
        return "Mes"
    return "Hind"


def zscore(values):
    values = np.asarray(values, dtype=float)
    std = np.nanstd(values)
    if not np.isfinite(std) or std == 0:
        return np.zeros_like(values, dtype=float)
    return (values - np.nanmean(values)) / std


def load_subject_causality(subject_id):
    return np.load(BASE_NET / f"subject_{subject_id}" / f"subject_{subject_id}_causality.npz")


def load_subject_fc_neighbors(subject_id):
    return np.load(
        BASE_NET
        / f"subject_{subject_id}"
        / f"subject_{subject_id}_net_te_drive_fc_neighbors.npz"
    )


def aggregate_subject_region_lists(values_by_subject, region_order):
    region_values = [[] for _ in region_order]
    order_lookup = {region_idx: pos for pos, region_idx in enumerate(region_order)}

    for region_num, values in values_by_subject:
        region_num = np.asarray(region_num)
        values = np.asarray(values, dtype=float)
        for reg_idx in region_order:
            mask = region_num == reg_idx
            if not np.any(mask):
                continue
            region_subject_values = values[mask]
            region_subject_values = region_subject_values[np.isfinite(region_subject_values)]
            if region_subject_values.size:
                region_values[order_lookup[reg_idx]].append(float(np.mean(region_subject_values)))

    return region_values


def ordered_lists_to_region_lists(region_order, ordered_lists):
    region_lists = [[] for _ in range(len(fs.region))]
    selected = []
    for region_idx, values in zip(region_order, ordered_lists):
        arr = np.asarray(values, dtype=float)
        arr = arr[np.isfinite(arr)]
        if len(arr) > 0:
            region_lists[int(region_idx)] = arr.tolist()
            selected.append(int(region_idx))
    return region_lists, selected


def replace_region_with_division_values(ordered_lists, region_order, target_region):
    target_idx = fs.region.index(target_region)
    if target_idx not in region_order:
        return ordered_lists

    target_pos = region_order.index(target_idx)
    target_division = fs.brain_division_list[target_idx]
    division_values = []

    for region_idx, values in zip(region_order, ordered_lists):
        if region_idx == target_idx:
            continue
        if fs.brain_division_list[region_idx] != target_division:
            continue
        arr = np.asarray(values, dtype=float)
        division_values.extend(arr[np.isfinite(arr)].tolist())

    if division_values:
        ordered_lists[target_pos] = division_values
    return ordered_lists


def draw_region_bar_main(ax, region_lists, selected_regions):
    means = []
    sems = []
    ordered_names = []
    ordered_divisions = []

    for idx in selected_regions:
        arr = np.asarray(region_lists[idx], dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            continue
        means.append(float(np.mean(arr)))
        sems.append(float(np.std(arr, ddof=1) / np.sqrt(arr.size)) if arr.size > 1 else 0.0)
        ordered_names.append(display_region_name(fs.region[idx]))
        ordered_divisions.append(int(fs.brain_division_list[idx]))

    x = 0
    xticks = []
    xticklabels = []
    xtick_divisions = []
    last_division = None
    for division in MAIN_DIVISION_ORDER:
        for i, node_division in enumerate(ordered_divisions):
            if node_division != division:
                continue
            if last_division is not None and node_division != last_division:
                ax.axvline(x - 0.5, color="#B5B5B5", lw=0.8, zorder=0)
            ax.bar(
                x,
                means[i],
                yerr=sems[i],
                color=MAIN_DIVISION_COLORS[division],
                alpha=0.86,
                width=0.62,
                error_kw=dict(
                    ecolor=MAIN_ERROR_COLOR,
                    elinewidth=0.8,
                    capsize=1.8,
                    capthick=0.8,
                ),
                zorder=2,
            )
            xticks.append(x)
            xticklabels.append(ordered_names[i])
            xtick_divisions.append(node_division)
            last_division = node_division
            x += 1

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=60, ha="right", fontsize=MAIN_TICK_FS)
    for lbl, division in zip(ax.get_xticklabels(), xtick_divisions):
        lbl.set_color(MAIN_DIVISION_COLORS.get(division, "#333333"))
        if division == 2:
            lbl.set_fontweight("bold")
            lbl.set_fontstyle("italic")
    ax.tick_params(axis="both", which="both", direction="out", bottom=True, left=True, length=3, width=0.8, labelsize=MAIN_TICK_FS)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(MAIN_AXIS_LW)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ============================================================
# Load Figure 9 node order
# ============================================================

fig9_summary_path = first_existing_path(
    DATA / "final_summary_tables" / "figure1_dynamic_fc_fingerprint_overview_values.csv",
)
fig9_summary = pd.read_csv(fig9_summary_path)
fig9_summary = fig9_summary.loc[fig9_summary["species"].eq("Zebrafish")].copy()
fig9_summary["display_node"] = fig9_summary["node"].astype(str).map(display_region_name)
fig9_summary["division"] = fig9_summary["display_node"].map(zebrafish_division)
fig9_summary["_division_order"] = fig9_summary["division"].map(DIVISION_ORDER)
fig9_summary["region_index"] = fig9_summary["display_node"].map(node_to_region_index)
fig9_summary = (
    fig9_summary.dropna(subset=["_division_order", "region_index"])
    .sort_values(["_division_order", "display_node"])
    .drop_duplicates("region_index")
    .reset_index(drop=True)
)
regions_for_heatmap = fig9_summary["display_node"].to_numpy()
region_order = fig9_summary["region_index"].astype(int).tolist()

def metric_lists_from_recording_table(value_col):
    table_path = first_existing_path(
        DATA / "source_inputs" / "ncomms_tables" / "highpass_ce_zf_plot_measures_recording_node.csv",
    )
    df = pd.read_csv(table_path)
    df = df.loc[df["species"].eq("Zebrafish")].replace([np.inf, -np.inf], np.nan)
    region_lists = []
    for region_name in regions_for_heatmap:
        values = df.loc[df["node"].astype(str).eq(region_name), value_col].to_numpy(float)
        values = values[np.isfinite(values)]
        region_lists.append(values.tolist())
    return region_lists


fcs_lists = metric_lists_from_recording_table("FCS")
fcv_lists = metric_lists_from_recording_table("EdgeStdFCV")
fc_partner_lists = metric_lists_from_recording_table("ProfileCorrDistFCV")

net_te_subject_data = []
neighbor_subject_data = []
for subject_id in SUBJECT_IDS:
    causality = load_subject_causality(subject_id)
    net_te_subject_data.append(
        (causality["region_num"], np.nanmean(causality["net_te_matrix"], axis=1))
    )

    neighbor = load_subject_fc_neighbors(subject_id)
    neighbor_subject_data.append(
        (causality["region_num"], neighbor["fc_neighbor_mean_drive"])
    )

net_te_lists = aggregate_subject_region_lists(net_te_subject_data, region_order)
neighbor_net_te_lists = aggregate_subject_region_lists(neighbor_subject_data, region_order)
net_te_lists = replace_region_with_division_values(net_te_lists, region_order, "rTS")
neighbor_net_te_lists = replace_region_with_division_values(
    neighbor_net_te_lists, region_order, "rTS"
)

feature_panels = [
    ("A", "FCV", fcv_lists),
    ("B", "FCS", fcs_lists),
    ("C", "FC Partner", fc_partner_lists),
    ("D", r"$\mathrm{TE}_{\mathrm{net}}$", net_te_lists),
    ("E", "Neighbor " + r"$\mathrm{TE}_{\mathrm{net}}$", neighbor_net_te_lists),
]


# ============================================================
# Create Figure Supply 10
# ============================================================

fig = plt.figure(figsize=(16, 16))
axes = [
    plt.subplot2grid((5, 6), (row, 0), colspan=6)
    for row in range(5)
]

for ax, (panel_label, ylabel, values) in zip(axes, feature_panels):
    region_lists, selected_regions = ordered_lists_to_region_lists(region_order, values)
    draw_region_bar_main(ax, region_lists, selected_regions)
    ax.set_ylabel(ylabel, fontsize=MAIN_AXIS_FS)
    ax.set_xlim(-1.1, len(selected_regions))
    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        bottom=True,
        left=True,
        length=3,
        width=0.8,
        labelsize=MAIN_TICK_FS,
    )
    ax.yaxis.set_label_coords(-0.055, 0.5)
    add_panel_label_fig(fig, ax, panel_label, dx=-0.075, dy=0.008, fontsize=MAIN_PANEL_FS)

for ax in (axes[0], axes[2]):
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0.2, ymax)

for ax in axes:
    pos = ax.get_position()
    ax.set_position([pos.x0 + 0.02, pos.y0, pos.width * 0.95, pos.height * 0.90])

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT_DIR / "figure_supply_10_proc.png", dpi=600, bbox_inches="tight")
