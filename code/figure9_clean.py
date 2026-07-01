import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

# Figure 9 workflow:
# 1. Load data and prepare derived quantities for plotting.
# 2. Prepare the figure layout and axes.
# 3. Draw each panel.
# 4. Adjust panel positions and add panel labels.
# 5. Save figure files and statistics.

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu, rankdata
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.colors as mcolors
import itertools

import figure_style as fs

def _strip_side(name):
    name = str(name)
    if len(name) > 1 and name[0] in {"l", "r"} and name[1].isupper():
        return name[1:]
    return name

def _panel_b_side_label(name):
    """Keep right-side root areas separate while matching left labels to root names."""
    name = str(name)
    if len(name) > 1 and name[0] == "l" and name[1].isupper():
        return name[1:]
    return name

def _zebrafish_anatomy_group(node):
    n = _strip_side(node).upper()
    if n in {"P", "SP", "OB", "OG", "OE", "PO"}:
        return "Tel"
    if n in {"HB", "HI", "HR", "TH", "PT", "PRT"}:
        return "Di"
    if n in {"TEO", "TL", "TS"}:
        return "Mes"
    if n in {"MON", "CB", "MOS1", "MOS2", "MOS3", "MOS4", "MOS5", "IO",
             "ARF", "IMRF", "PRF", "TG", "VR", "NX", "HC", "IPN", "RA", "T"}:
        return "Hind"
    return "Hind"


lregion = ['lMON','lCb','lMOS1','lMOS2','lMOS3','lMOS4','lMOS5','lIPN','lIO','lHc','lRa','lT',
           'laRF','limRF','lpRF','lGG','lHb','lHi','lHR','lOG','lOB','lOE','lP','lPi','lPT',
           'lPO','lPrT','lR','lSP','lTeO','lTh','lTL','lTS','lTG','lVR','lNX']

rregion = ['rMON','rCb','rMOS1','rMOS2','rMOS3','rMOS4','rMOS5','rIPN','rIO','rHc',
           'rRa','rT','raRF','rimRF','rpRF','rGG','rHb','rHi','rHR','rOG','rOB',
           'rOE','rP','rPi','rPT','rPO','rPrT','rR','rSP','rTeO','rTh','rTL',
           'rTS','rTG','rVR','rNX']

region = lregion + rregion
SUBJECT_IDS = range(12, 19)
FIG1_FC_FINGERPRINT_TABLE = "final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv"
FIG1_RECORDING_MEASURE_TABLE = "source_inputs/ncomms_tables/highpass_ce_zf_plot_measures_recording_node.csv"
OBSERVED_NETTE_RECORDING_TABLE = "final_summary_tables/observed_nette_no_p_recording_level.csv"
ZF_GROUP_ORDER = ["Tel", "Di", "Mes", "Hind"]
ZF_GROUP_COLORS = fs.ZEBRAFISH_DIVISION_COLORS.copy()
fs.apply_main_figure_style()
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(PROJECT_ROOT, "data")
OUTPUT_PNG = os.path.join(PROJECT_ROOT, "figures", "figure9_final.png")
STATS_DIR = os.path.join(PROJECT_ROOT, "data", "final_summary_tables")
STATS_CSV = os.path.join(STATS_DIR, "figure9_stats.csv")
BASE_NET = os.path.join(DATA, "region_community_io")
STATS_ROWS = []


# ============================================================
# Helpers
# ============================================================
def _zscore(values):
    values = np.asarray(values, dtype=float)
    return (values - np.nanmean(values)) / np.nanstd(values)

def _reorder_linkage_tel_left_hind_right(z_linkage, divisions):
    """Flip dendrogram branches while preserving clustering topology.

    Branches enriched for Tel are placed left when possible, and branches
    enriched for Hind are placed right. This does not force division blocks;
    it only chooses the left/right orientation of existing dendrogram branches.
    """
    reordered = np.array(z_linkage, copy=True)
    n_leaves = len(divisions)
    divisions = np.asarray(divisions).astype(str)
    stats_by_id = {}
    display_rank = {"Tel": 0.0, "Di": 1.0, "Mes": 2.0, "Hind": 5.0}
    for idx, div in enumerate(divisions):
        stats_by_id[idx] = {
            "n": 1,
            "tel": int(div == "Tel"),
            "hind": int(div == "Hind"),
            "rank_sum": float(display_rank.get(div, max(display_rank.values()) + 1.0)),
        }

    def _score(stats):
        n = max(int(stats["n"]), 1)
        mean_rank = stats["rank_sum"] / n
        # Smaller scores are placed to the left. Mean division rank encourages
        # the left-to-right progression Tel -> Di -> Mes -> Hind, so Mes-rich
        # branches tend to stay between forebrain and hindbrain clusters while
        # preserving the original Ward clustering topology.
        hind_minus_tel = (stats["hind"] - stats["tel"]) / n
        return (mean_rank, hind_minus_tel)

    for row_idx, row in enumerate(reordered):
        left = int(row[0])
        right = int(row[1])
        parent = n_leaves + row_idx

        left_stats = stats_by_id[left]
        right_stats = stats_by_id[right]
        if _score(right_stats) < _score(left_stats):
            reordered[row_idx, 0], reordered[row_idx, 1] = reordered[row_idx, 1], reordered[row_idx, 0]
            left, right = right, left
            left_stats, right_stats = right_stats, left_stats

        stats_by_id[parent] = {
            "n": left_stats["n"] + right_stats["n"],
            "tel": left_stats["tel"] + right_stats["tel"],
            "hind": left_stats["hind"] + right_stats["hind"],
            "rank_sum": left_stats["rank_sum"] + right_stats["rank_sum"],
        }

    return reordered

def _move_last_cluster_run_to_middle(z_linkage, leaf_order, n_clusters=3):
    """Cut leaves into contiguous cluster runs and move the last run to middle."""
    leaf_order = np.asarray(leaf_order, dtype=int)
    cluster_labels = fcluster(z_linkage, t=n_clusters, criterion="maxclust")
    runs = []
    start = 0
    current_label = cluster_labels[leaf_order[0]]
    for pos, leaf_idx in enumerate(leaf_order[1:], start=1):
        label = cluster_labels[leaf_idx]
        if label != current_label:
            runs.append(leaf_order[start:pos])
            start = pos
            current_label = label
    runs.append(leaf_order[start:])
    if len(runs) != n_clusters:
        return leaf_order
    return np.concatenate([runs[0], runs[-1], *runs[1:-1]])

def _load_subject_causality(subject_id):
    return np.load(f"{BASE_NET}/subject_{subject_id}/subject_{subject_id}_causality.npz")

def _add_figure9_panel_labels(ax_dendro, ax_hier, ax_c, ax_d, ax_e, ax_f, ax_g):
    ax_dendro.text(
        -0.05, 1.5, 'A',
        transform=ax_dendro.transAxes,
        fontsize=fs.PANEL_LABEL_FS_2COL,
        fontweight='bold',
        va='bottom',
        ha='right',
    )
    
    ax_hier.text(
        -0.02, 0.98, 'B',
        transform=ax_hier.transAxes,
        fontsize=fs.PANEL_LABEL_FS_2COL,
        fontweight='bold',
        va='bottom',
        ha='right',
    )
    
    for ax, label in [
        (ax_c, 'C'),
        (ax_d, 'D'),
        (ax_e, 'E'),
        (ax_f, 'F'),
        (ax_g, 'G'),
    ]:
        ax.text(
            -0.31, 1.05, label,
            transform=ax.transAxes,
            fontsize=fs.PANEL_LABEL_FS_2COL,
            fontweight='bold',
            va='bottom',
        )


def _compute_panel_te_from_causality():
    rows = []
    for subject_id in SUBJECT_IDS:
        subject_data = _load_subject_causality(subject_id)
        rn = subject_data["region_num"]
        nte = subject_data["net_te_matrix"]
        fc_sig = subject_data["fc_neighbor_mask_fdr"]
        comm_nodes = np.array([
            _panel_b_side_label(region[int(reg_idx)]) if int(reg_idx) < len(region) else ""
            for reg_idx in rn
        ])
        unique_nodes = [name for name in np.unique(comm_nodes) if name]
        for src_name in unique_nodes:
            src_comm = np.where(comm_nodes == src_name)[0]
            all_targets = []
            fc_targets = []
            for tgt_name in unique_nodes:
                if src_name == tgt_name:
                    continue
                tgt_comm = np.where(comm_nodes == tgt_name)[0]
                sub_nte = nte[np.ix_(src_comm, tgt_comm)]
                sub_nte = sub_nte[np.isfinite(sub_nte)]
                if len(sub_nte):
                    target_value = float(np.mean(sub_nte))
                    all_targets.append(target_value)
                    if np.any(fc_sig[np.ix_(src_comm, tgt_comm)]):
                        fc_targets.append(target_value)
            rows.append({
                "species": "Zebrafish",
                "recording_id": f"subject_{subject_id}",
                "node": src_name,
                "ObservedNetTE_computed": float(np.mean(all_targets)) if all_targets else np.nan,
                "NeighborNetTE_computed": float(np.mean(fc_targets)) if fc_targets else np.nan,
            })
    return pd.DataFrame(rows)


# ============================================================
# 1. Data loading and plotting calculations
# ============================================================
# Figure 1 zebrafish FC-dynamics table.  These values are already the
# finalized, within-species z-scored node summaries used in Figure 1.
fig1_fc_df = pd.read_csv(os.path.join(DATA, FIG1_FC_FINGERPRINT_TABLE))
fig1_zf = (
    fig1_fc_df.loc[fig1_fc_df["species"].eq("Zebrafish")]
    .replace([np.inf, -np.inf], np.nan)
    .dropna(subset=["EdgeStdFCV", "FCS", "ProfileCorrDistFCV", "ObservedNetTE", "NeighborNetTE"])
    .copy()
)
fig1_zf["anatomy_group"] = fig1_zf["node"].map(_zebrafish_anatomy_group)
fig1_zf["_group_order"] = fig1_zf["anatomy_group"].map(
    {group: idx for idx, group in enumerate(ZF_GROUP_ORDER)}
).fillna(len(ZF_GROUP_ORDER)).astype(int)
fig1_zf = fig1_zf.sort_values(["_group_order", "anatomy_group", "node"]).reset_index(drop=True)

fig1_rec = pd.read_csv(os.path.join(DATA, FIG1_RECORDING_MEASURE_TABLE))
fig1_rec = fig1_rec.loc[fig1_rec["species"].eq("Zebrafish")].replace([np.inf, -np.inf], np.nan).copy()
fig1_rec = fig1_rec[["species", "recording_id", "node", "EdgeStdFCV", "FCS", "ProfileCorrDistFCV"]]
nette_rec = pd.read_csv(os.path.join(DATA, OBSERVED_NETTE_RECORDING_TABLE))
nette_rec = (
    nette_rec.loc[nette_rec["species"].eq("Zebrafish"),
                  ["species", "recording_id", "node", "NetTE", "NeighborNetTE"]]
    .rename(columns={"NetTE": "ObservedNetTE"})
)
fig1_zf_rec = fig1_rec.merge(nette_rec, on=["species", "recording_id", "node"], how="left")
_panel_te_rec = _compute_panel_te_from_causality()
fig1_zf_rec = fig1_zf_rec.merge(_panel_te_rec, on=["species", "recording_id", "node"], how="left")
for _col in ["ObservedNetTE", "NeighborNetTE"]:
    fig1_zf_rec[_col] = fig1_zf_rec[_col].fillna(fig1_zf_rec[f"{_col}_computed"])
fig1_zf_rec = fig1_zf_rec.drop(columns=["ObservedNetTE_computed", "NeighborNetTE_computed"])
_measure_cols = ["EdgeStdFCV", "FCS", "ProfileCorrDistFCV", "ObservedNetTE", "NeighborNetTE"]
_zrec_frames = []
for _, _group in fig1_zf_rec.groupby(["species", "recording_id"], sort=False):
    _out = _group.copy()
    for _col in _measure_cols:
        _out[_col] = _zscore(_out[_col].to_numpy(float))
    _zrec_frames.append(_out)
fig1_zf_rec = pd.concat(_zrec_frames, ignore_index=True)
fig1_zf_rec["anatomy_group"] = fig1_zf_rec["node"].map(_zebrafish_anatomy_group)
fig1_zf_rec["_group_order"] = fig1_zf_rec["anatomy_group"].map(
    {group: idx for idx, group in enumerate(ZF_GROUP_ORDER)}
).fillna(len(ZF_GROUP_ORDER)).astype(int)
fig1_zf_rec = fig1_zf_rec.sort_values(
    ["_group_order", "anatomy_group", "recording_id", "node"]
).reset_index(drop=True)

def _add_panel_a_missing_root_nodes(summary_df, recording_df, nodes):
    """Add root nodes missing from the final summary, filling absent measures by division mean."""
    rows = []
    existing = set(summary_df["node"].astype(str))
    for node in nodes:
        if node in existing:
            continue
        node_division = _zebrafish_anatomy_group(node)
        node_rec = recording_df.loc[recording_df["node"].astype(str).eq(node)]
        division_summary = summary_df.loc[summary_df["anatomy_group"].astype(str).eq(node_division)]

        row = {col: np.nan for col in summary_df.columns}
        row.update({
            "species": "Zebrafish",
            "node": node,
            "anatomy_group": node_division,
            "level": "root_area_from_louvain_communities",
        })
        if "n_recordings" in row:
            row["n_recordings"] = int(node_rec["recording_id"].nunique()) if not node_rec.empty else 0

        for col in _measure_cols:
            value = float(node_rec[col].mean()) if (not node_rec.empty and col in node_rec) else np.nan
            if not np.isfinite(value):
                value = float(division_summary[col].mean()) if col in division_summary else np.nan
            row[col] = value
        rows.append(row)

    if not rows:
        return summary_df

    out = pd.concat([summary_df, pd.DataFrame(rows)], ignore_index=True)
    out["_group_order"] = out["anatomy_group"].map(
        {group: idx for idx, group in enumerate(ZF_GROUP_ORDER)}
    ).fillna(len(ZF_GROUP_ORDER)).astype(int)
    return out.sort_values(["_group_order", "anatomy_group", "node"]).reset_index(drop=True)

# The final node-summary table collapses OB to a single node. For Panel A and
# Panel B, restore rOB from the recording-level table where available; measures
# missing for rOB, such as TE summaries, are filled with the corresponding
# division mean so the node can be displayed without dropping complete rows.
fig1_zf = _add_panel_a_missing_root_nodes(fig1_zf, fig1_zf_rec, ["rOB"])

fig1_zf_panel_b = fig1_zf.copy()
_panel_summary_missing_nodes = ["rOB"]
_panel_summary_rows = []
for _node in _panel_summary_missing_nodes:
    if _node in set(fig1_zf_panel_b["node"].astype(str)):
        continue
    _node_rec = fig1_zf_rec.loc[fig1_zf_rec["node"].astype(str).eq(_node)]
    if _node_rec.empty:
        continue
    _row = {col: np.nan for col in fig1_zf.columns}
    _row.update({
        "species": "Zebrafish",
        "node": _node,
        "anatomy_group": _zebrafish_anatomy_group(_node),
    })
    for _col in _measure_cols:
        _value = float(_node_rec[_col].mean())
        _row[_col] = _value if np.isfinite(_value) else 0.0
    _panel_summary_rows.append(_row)

if _panel_summary_rows:
    fig1_zf_panel_b = pd.concat([fig1_zf_panel_b, pd.DataFrame(_panel_summary_rows)], ignore_index=True)
    fig1_zf_panel_b["_group_order"] = fig1_zf_panel_b["anatomy_group"].map(
        {group: idx for idx, group in enumerate(ZF_GROUP_ORDER)}
    ).fillna(len(ZF_GROUP_ORDER)).astype(int)
    fig1_zf_panel_b = fig1_zf_panel_b.sort_values(["_group_order", "anatomy_group", "node"]).reset_index(drop=True)

regions_d = fig1_zf["node"].astype(str).to_numpy()
divs_d = fig1_zf["anatomy_group"].astype(str).to_numpy()
fcv_vals = fig1_zf["EdgeStdFCV"].to_numpy(float)
fcs_vals = fig1_zf["FCS"].to_numpy(float)
profile_corr_dist_vals = fig1_zf["ProfileCorrDistFCV"].to_numpy(float)
net_transfer_data = fig1_zf["ObservedNetTE"].to_numpy(float)
fc_neighbor_mean_drive_data = fig1_zf["NeighborNetTE"].to_numpy(float)

out_data = np.array([
    fcv_vals,
    fcs_vals,
    profile_corr_dist_vals,
    net_transfer_data,
    fc_neighbor_mean_drive_data,
])

div_color_map = {
    **ZF_GROUP_COLORS,
}

# ── 클러스터링: out_data.T = (n_regions × 5 features) ──
# Panel A preserves Ward clustering of the measured FC features, but flips
# dendrogram branches so clusters enriched for Tel tend to appear on the left
# and Hind-enriched clusters on the right when this is possible within the
# existing hierarchy.
Z = linkage(out_data.T, method='ward')
Z = _reorder_linkage_tel_left_hind_right(Z, divs_d)

# 덴드로그램 leaf 순서 추출
dend_info  = dendrogram(Z, no_plot=True)
leaf_order = np.array(dend_info['leaves'])
leaf_order = _move_last_cluster_run_to_middle(Z, leaf_order, n_clusters=3)

# 클러스터링 순서로 데이터 재정렬
out_data_c = out_data[:, leaf_order]
regions_c  = regions_d[leaf_order]
divs_c     = divs_d[leaf_order]
n_regions = out_data_c.shape[1]
y_labels  = [
    'FCV',
    'FCS',
    'FC Recon. \nDeg.',
    r'$\mathbf{TE}_{\mathbf{net}}$',
    'Neighbor\n' + r'$\mathbf{TE}_{\mathbf{net}}$',
]

# ============================================================
# 2. Layout preparation
# ============================================================
# Left=[A (dendro+divbar+heat)] | Right=[B (network)] | Bottom=[C D E F G]
_fig_w   = 16;#fs.TWO_COL_IN
_fig_h   = 9;#fs.TWO_COL_IN * 0.65
fig = plt.figure(figsize=(_fig_w, _fig_h))
gs  = GridSpec(5, 10,
               figure=fig,
               height_ratios=[1.0, 0.2, 2.5, 2.0, 2.0],
               width_ratios=[1]*10,
               left=0.08, right=0.97,
               top=0.94, bottom=0.20,
               hspace=0.20, wspace=0.30)

# ──── LEFT COLUMN: Panel A (Dendrogram, Division bar, Heatmap) ────
# Row 0: Dendrogram (left, cols 0-5)
ax_dendro = fig.add_subplot(gs[0, 0:4])

# Row 1: Division bar (left, cols 0-5)
ax_divbar_a = fig.add_subplot(gs[1, 0:4])

# Row 2: Heatmap (left, cols 0-4), Colorbar (col 4)
ax_heat   = fig.add_subplot(gs[2, 0:4])
ax_cbar   = fig.add_subplot(gs[2, 4:5])

# ──── RIGHT TOP: Panel B (Hierarchical network at top) ────
ax_hier = fig.add_subplot(gs[0:3, 5:10])

# ──── BOTTOM ROW: Panels C, D, E, F, G (각 2 columns, 동일 크기) ────
ax_c = fig.add_subplot(gs[3:5, 0:2])    # Panel C
ax_d = fig.add_subplot(gs[3:5, 2:4])    # Panel D
ax_e = fig.add_subplot(gs[3:5, 4:6])    # Panel E
ax_f = fig.add_subplot(gs[3:5, 6:8])    # Panel F
ax_g = fig.add_subplot(gs[3:5, 8:10])   # Panel G

# Panel C-G plotting helpers and source tables

# ── figure2 helper functions ──
_div_colors_f2 = ZF_GROUP_COLORS
_div_order_f2  = [group for group in ZF_GROUP_ORDER if group in set(fig1_zf["anatomy_group"])]
_div_short_labels_f2 = {
    "Tel": "Tel",
    "Di": "Di",
    "Mes": "Mes",
    "Hind": "Hind",
}
_STAR_FS = fs.STAR_FS_2COL

def _record_division_stats(panel, df, order):
    groups = [df[k].dropna().values for k in order]
    pairs  = list(itertools.combinations(range(len(order)), 2))
    pvals  = [mannwhitneyu(groups[i], groups[j], alternative='two-sided')[1]
              for i, j in pairs]
    reject, corr, _, _ = multipletests(pvals, method='holm')
    for idx, (i, j) in enumerate(pairs):
        STATS_ROWS.append({
            'figure': 'figure9',
            'panel': panel,
            'test': 'Mann-Whitney U',
            'alternative': 'two-sided',
            'group_1': order[i],
            'group_2': order[j],
            'n_group_1': len(groups[i]),
            'n_group_2': len(groups[j]),
            'p_uncorrected': pvals[idx],
            'p_holm': corr[idx],
            'reject_holm_0.05': bool(reject[idx]),
        })
    return pairs, reject, corr

def _add_sig_bars(ax, df, order, panel):
    pairs, reject, corr = _record_division_stats(panel, df, order)
    sig = sorted([(pairs[k], corr[k]) for k in range(len(pairs)) if reject[k]],
                 key=lambda x: x[0][1] - x[0][0])
    if not sig:
        return
    y_min, y_max = ax.get_ylim()
    yr    = y_max - y_min
    step  = yr * 0.090
    bar_h = yr * 0.018
    for lvl, ((i, j), p) in enumerate(sig):
        y    = y_max + yr * 0.025 + lvl * step
        star = '***' if p < 0.001 else '**' if p < 0.01 else '*'
        ax.plot(
            [i, i, j, j],
            [y, y + bar_h, y + bar_h, y],
            lw=0.65,
            c='#333',
            clip_on=False,
        )
        ax.text((i + j) / 2, y + bar_h, star,
                ha='center', va='center', fontsize=_STAR_FS,
                clip_on=False)
    ax.set_ylim(y_min, y_max)

def _boxplot_panel(ax, df, ylabel):
    vals, labels = [], []
    for col in _div_order_f2:
        v = df[col].dropna().values
        vals.extend(v); labels.extend([col] * len(v))
    fs.draw_main_box_strip(
        ax,
        labels,
        vals,
        _div_order_f2,
        palette=[_div_colors_f2[d] for d in _div_order_f2],
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel('')
    ax.set_xticks(range(len(_div_order_f2)))
    ax.set_xticklabels([_div_short_labels_f2.get(group, group) for group in _div_order_f2],
                       rotation=0, ha='center')
    #_add_sig_bars(ax, {c: df[c] for c in _div_order_f2}, _div_order_f2, ylabel)

def _fig1_measure_division_df(plot_df, measure):
    plot_df = plot_df[['anatomy_group', measure]].dropna().copy()
    grouped = {
        div: plot_df.loc[plot_df['anatomy_group'] == div, measure].to_numpy()
        for div in _div_order_f2
    }
    max_len = max((len(v) for v in grouped.values()), default=0)
    return pd.DataFrame({
        div: pd.Series(vals, dtype=float).reindex(range(max_len))
        for div, vals in grouped.items()
    })

# ── Figure 1 zebrafish data (Panels C-G) ──
_fcv = _fig1_measure_division_df(fig1_zf_rec, 'EdgeStdFCV')
_fcs = _fig1_measure_division_df(fig1_zf_rec, 'FCS')
_profile_corr = _fig1_measure_division_df(fig1_zf_rec, 'ProfileCorrDistFCV')
out_net_te_data_df = _fig1_measure_division_df(fig1_zf_rec, 'ObservedNetTE')
out_neigh_net_te_data_df = _fig1_measure_division_df(fig1_zf_rec, 'NeighborNetTE')

# ============================================================
# 3. Draw each panel
# ============================================================
# Panels C-G: division-level boxplots
_boxplot_panel(ax_c, _fcv, 'FCV')
_boxplot_panel(ax_d, _fcs, 'FCS')
_boxplot_panel(ax_e, _profile_corr, 'FC Reconfig. Deg')
_boxplot_panel(ax_f, out_net_te_data_df, r'$\mathrm{TE}_{\mathrm{net}}$')
_boxplot_panel(ax_g, out_neigh_net_te_data_df, 'Neighbor ' + r'$\mathrm{TE}_{\mathrm{net}}$')


_add_sig_bars(ax_c, {c:  _fcv[c] for c in _div_order_f2}, _div_order_f2, 'FCV')
_add_sig_bars(ax_d, {c:  _fcs[c] for c in _div_order_f2}, _div_order_f2, 'FCS')
_add_sig_bars(ax_e, {c:  _profile_corr[c] for c in _div_order_f2}, _div_order_f2, 'ProfileCorrDistFCV')
_add_sig_bars(ax_f, {c:  out_net_te_data_df[c] for c in _div_order_f2}, _div_order_f2, r'$\mathrm{TE}_{\mathrm{net}}$')
ax_g.set_ylim(-3.0, ax_g.get_ylim()[1] * 1.15)
_add_sig_bars(ax_g, {c:  out_neigh_net_te_data_df[c] for c in _div_order_f2}, _div_order_f2, 'Neighbor ' + r'$\mathrm{TE}_{\mathrm{net}}$')

# Panel A: dendrogram, division bar, and feature heatmap
dendrogram(Z, ax=ax_dendro, no_labels=True,
           color_threshold=0, above_threshold_color='#333333')
ax_dendro.set_xlim(0, n_regions * 10)
# 상단 여백 제거: 실제 최대 높이 + 5%만 표시
dend_max_h = max(max(d) for d in dend_info['dcoord'])
ax_dendro.set_ylim(0, dend_max_h * 1.05)
ax_dendro.axis('off')


# ── 2. Division bar (full height) ──
for i, div in enumerate(divs_c):
    ax_divbar_a.add_patch(plt.Rectangle((i - 0.5, 0), 1, 1,
                        color=div_color_map.get(div, 'gray'), linewidth=0))

ax_divbar_a.set_xlim(-0.5, n_regions - 0.5)
ax_divbar_a.set_ylim(0, 1)
ax_divbar_a.set_xticks([])
ax_divbar_a.set_yticks([])
for spine in ax_divbar_a.spines.values():
    spine.set_visible(False)

# ── 3. Heatmap ──
im = ax_heat.imshow(out_data_c, aspect='auto', cmap='RdBu_r',
                    vmax=2, vmin=-2, interpolation='nearest')
ax_heat.set_xlim(-0.5, n_regions - 0.5)

ax_heat.set_yticks(range(len(y_labels)))
ax_heat.set_yticklabels(y_labels, fontsize=fs.TICK_FS_2COL, fontweight='bold')
ax_heat.tick_params(axis='y', length=0, labelright=True, labelleft=False)

ax_heat.set_xticks(range(n_regions))
ax_heat.set_xticklabels(regions_c, rotation=90, fontsize=fs.TICK_FS_2COL - 2,
                         ha='center', fontweight='bold')
ax_heat.tick_params(axis='x', length=2, width=0.5)

# x-tick 레이블 색: division 색
for tick_label, div in zip(ax_heat.get_xticklabels(), divs_c):
    tick_label.set_color(div_color_map[div])

# 행 구분선
for y in np.arange(-0.5, len(y_labels), 1):
    ax_heat.axhline(y=y, color='white', linewidth=0.4)

_highlight_regions = {"P", "rP", "SP", "rSP","OB","rOB"}
_highlight_idx = [i for i, name in enumerate(regions_c) if str(name) in _highlight_regions]
if _highlight_idx:
    _x0 = min(_highlight_idx) - 0.5
    _width = max(_highlight_idx) - min(_highlight_idx) + 1
    _highlight_edge = "#0AEB60FF"
    ax_divbar_a.add_patch(
        mpatches.Rectangle(
            (_x0, 0),
            _width,
            1,
            fill=False,
            edgecolor=_highlight_edge,
            linewidth=1.8,
            zorder=20,
            clip_on=False,
        )
    )
    ax_heat.add_patch(
        mpatches.Rectangle(
            (_x0, -0.5),
            _width,
            len(y_labels),
            fill=False,
            edgecolor=_highlight_edge,
            linewidth=1.8,
            zorder=25,
            clip_on=False,
        )
    )

# 외곽 테두리
for spine in ax_heat.spines.values():
    spine.set_linewidth(1.4)
    spine.set_color('black')
ax_heat.add_patch(
    mpatches.Rectangle(
        (-0.5, -0.5), n_regions, len(y_labels),
        fill=False, edgecolor='black', linewidth=1.4,
        zorder=10, clip_on=False
    )
)

# ── 4. Colorbar ──
ax_cbar.set_position([0.48, 0.585, 0.020, 0.205])
                     
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.ax.set_title('Z-score', fontsize=fs.AXIS_LABEL_FS_2COL, fontweight='bold')
cbar.ax.tick_params(labelsize=fs.TICK_FS_2COL)
cbar.set_ticks([-2, -1, 0, 1, 2])


# ── 5. 통합 범례 (heatmap 하단 바깥) ──
div_handles = [mpatches.Patch(color=div_color_map[d], label=d)
               for d in _div_order_f2]
leg = ax_heat.legend(handles=div_handles,
                     loc='upper left',
                     bbox_to_anchor=(0.35, 1.95),
                     ncol=4, fontsize=fs.TICK_FS_2COL,
                     frameon=False, framealpha=0.95,
                     title_fontsize=fs.AXIS_LABEL_FS_2COL)
leg.get_title().set_fontweight('bold')
for text in leg.get_texts():
    text.set_fontweight('bold')

# Panel B: directed FC+TE network across zebrafish root areas.
_white = np.array([1.0, 1.0, 1.0])
net_nodes = fig1_zf_panel_b[["node", "anatomy_group", "EdgeStdFCV", "ObservedNetTE"]].copy()
net_nodes["node"] = net_nodes["node"].astype(str).map(_panel_b_side_label)
net_nodes = (
    net_nodes.drop_duplicates("node")
    .sort_values(["anatomy_group", "node"])
    .reset_index(drop=True)
)
if "rOB" not in set(net_nodes["node"].astype(str)):
    rob_rec = fig1_zf_rec.loc[fig1_zf_rec["node"].astype(str).eq("rOB")]
    if not rob_rec.empty:
        rob_row = {
            "node": "rOB",
            "anatomy_group": "Tel",
            "EdgeStdFCV": float(rob_rec["EdgeStdFCV"].mean()),
            "ObservedNetTE": float(net_nodes.loc[net_nodes["anatomy_group"].eq("Tel"), "ObservedNetTE"].mean()),
        }
        net_nodes = pd.concat([net_nodes, pd.DataFrame([rob_row])], ignore_index=True)
        net_nodes = net_nodes.sort_values(["anatomy_group", "node"]).reset_index(drop=True)
node_names = net_nodes["node"].to_numpy()
node_groups = net_nodes["anatomy_group"].to_numpy()
node_to_pos = {name: idx for idx, name in enumerate(node_names)}
n_nodes_b = len(node_names)
PANEL_B_NODE_ROWS = net_nodes.copy()

te_sum_b = np.zeros((n_nodes_b, n_nodes_b), dtype=float)
te_cnt_b = np.zeros((n_nodes_b, n_nodes_b), dtype=int)
fc_sig_cnt_b = np.zeros((n_nodes_b, n_nodes_b), dtype=int)

for subject_id in SUBJECT_IDS:
    subject_data = _load_subject_causality(subject_id)
    rn = subject_data["region_num"]
    comm_nodes = np.array([
        _panel_b_side_label(region[int(reg_idx)]) if int(reg_idx) < len(region) else ""
        for reg_idx in rn
    ])
    unique_nodes = [name for name in np.unique(comm_nodes) if name in node_to_pos]
    nte = subject_data["net_te_matrix"]
    fc_sig = subject_data["fc_neighbor_mask_fdr"]
    for src_name in unique_nodes:
        src_comm = np.where(comm_nodes == src_name)[0]
        src_pos = node_to_pos[src_name]
        for tgt_name in unique_nodes:
            if src_name == tgt_name:
                continue
            tgt_comm = np.where(comm_nodes == tgt_name)[0]
            tgt_pos = node_to_pos[tgt_name]
            sub_nte = nte[np.ix_(src_comm, tgt_comm)]
            sub_nte = sub_nte[np.isfinite(sub_nte)]
            if len(sub_nte):
                te_sum_b[src_pos, tgt_pos] += float(np.mean(sub_nte))
                te_cnt_b[src_pos, tgt_pos] += 1
            if np.any(fc_sig[np.ix_(src_comm, tgt_comm)]):
                fc_sig_cnt_b[src_pos, tgt_pos] += 1

te_mean_b = np.divide(te_sum_b, np.maximum(te_cnt_b, 1), where=te_cnt_b > 0)
te_mean_b[te_cnt_b == 0] = np.nan
fc_edge_mask_b = fc_sig_cnt_b >= int(np.ceil(len(list(SUBJECT_IDS)) / 2))
finite_fc_te = te_mean_b[np.isfinite(te_mean_b) & fc_edge_mask_b]
positive_fc_te = finite_fc_te[finite_fc_te > 0]
_te_thresh_b = float(np.percentile(positive_fc_te, 85)) if len(positive_fc_te) else np.inf

hier_edges = []
for src_idx in range(n_nodes_b):
    for tgt_idx in range(n_nodes_b):
        if src_idx == tgt_idx or not fc_edge_mask_b[src_idx, tgt_idx]:
            continue
        weight = te_mean_b[src_idx, tgt_idx]
        if np.isfinite(weight) and weight >= _te_thresh_b:
            hier_edges.append((src_idx, tgt_idx, float(weight)))
hier_edges.sort(key=lambda item: item[2])

# Visual cap: TS receives many selected edges and otherwise dominates panel B.
# Keep only the strongest incoming TS edges for the illustrative network.
_target_incoming_cap = {"TS": 6}
_target_counts = {target: 0 for target in _target_incoming_cap}
_capped_edges = []
for src, tgt, weight in sorted(hier_edges, key=lambda item: item[2], reverse=True):
    target_name = node_names[tgt]
    if target_name in _target_incoming_cap:
        if _target_counts[target_name] >= _target_incoming_cap[target_name]:
            continue
        _target_counts[target_name] += 1
    _capped_edges.append((src, tgt, weight))
hier_edges = sorted(_capped_edges, key=lambda item: item[2])

PANEL_B_EDGE_ROWS = [
    {
        "source": node_names[src],
        "target": node_names[tgt],
        "source_group": node_groups[src],
        "target_group": node_groups[tgt],
        "mean_net_te": weight,
        "n_subjects_te": int(te_cnt_b[src, tgt]),
        "n_subjects_fc_sig": int(fc_sig_cnt_b[src, tgt]),
    }
    for src, tgt, weight in hier_edges
]

div_y_hier = {"Tel": 3, "Di": 2, "Mes": 1, "Hind": 0}
hier_pos = {}
for div in ZF_GROUP_ORDER:
    idxs = [idx for idx, group in enumerate(node_groups) if group == div]
    idxs = sorted(idxs, key=lambda idx: node_names[idx])
    for j, idx in enumerate(idxs):
        x = (j / max(len(idxs) - 1, 1)) * 2.0 - 1.0
        hier_pos[idx] = (x, float(div_y_hier[div]))

fcv_rank = rankdata(net_nodes["EdgeStdFCV"].to_numpy(float))
fcv_norm = (fcv_rank - 1) / max(len(fcv_rank) - 1, 1)

if hier_edges:
    te_vals = np.array([edge[2] for edge in hier_edges], dtype=float)
    te_min = float(np.nanmin(te_vals))
    te_max = float(np.nanmax(te_vals))
    te_rng = max(te_max - te_min, 1e-9)
    for src, tgt, weight in hier_edges:
        x0, y0 = hier_pos[src]
        x1, y1 = hier_pos[tgt]
        norm_w = (weight - te_min) / te_rng
        same_layer = abs(y1 - y0) < 0.1
        rad = 0.17 if same_layer else 0.045
        edge_color = div_color_map.get(node_groups[src], "#999999")
        ax_hier.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>",
                color=edge_color,
                lw=0.6 + 1.6 * norm_w,
                alpha=0.22 + 0.38 * norm_w,
                mutation_scale=10 + 7 * norm_w,
                shrinkA=8,
                shrinkB=8,
                connectionstyle=f"arc3,rad={rad}",
            ),
            zorder=2,
        )

for idx, name in enumerate(node_names):
    if idx not in hier_pos:
        continue
    x, y = hier_pos[idx]
    group = node_groups[idx]
    base_rgb = np.array(mcolors.to_rgb(div_color_map.get(group, "#999999")))
    intensity = 0.18 + 0.82 * float(fcv_norm[idx])
    color = np.clip(_white * (1 - intensity) + base_rgb * intensity, 0, 1)
    size = 18 + 185 * float(fcv_norm[idx])
    ax_hier.scatter(x, y, s=size, c=[color],
                    edgecolors="black", linewidths=0.25, zorder=4)
    ax_hier.text(
        x, y - 0.18, name,
        fontsize=fs.TICK_FS_2COL * 0.9,
        ha="center", va="top", rotation=90,
        color=div_color_map.get(group, "#333333"),
        fontweight="bold", zorder=5,
        path_effects=[pe.withStroke(linewidth=3.0, foreground="white")],
    )

for div, y in div_y_hier.items():
    ax_hier.text(-1.12, y, div, fontsize=fs.TICK_FS_2COL * 1.2,
                 ha="right", va="center", color=div_color_map[div],
                 fontweight="bold")

ax_hier.set_xlim(-1.18, 1.08)
ax_hier.set_ylim(-1.20, 3.80)
_hier_pos = ax_hier.get_position()
ax_hier.set_position([
    _hier_pos.x0,
    _hier_pos.y0 - 0.070,
    _hier_pos.width,
    _hier_pos.height * 1.35,
])
ax_hier.set_aspect("auto")
ax_hier.axis("off")




# ============================================================
# 4. Panel position adjustment and panel labels
# ============================================================
_cb_pos = ax_cbar.get_position()
ax_cbar.set_position([_cb_pos.x0, _cb_pos.y0, _cb_pos.width * 0.4, _cb_pos.height])

# 패널 C-G 크기 30% 축소 (중심 유지)
for _ax in [ax_c, ax_d, ax_e, ax_f, ax_g]:
    _p  = _ax.get_position()
    _cx = _p.x0 + _p.width  / 2
    _cy = _p.y0 + _p.height / 2
    _nw = _p.width  * 0.8
    _nh = _p.height * 0.45
    _ax.set_position([_cx - _nw/2, _cy - _nh/2, _nw, _nh])

_add_figure9_panel_labels(ax_dendro, ax_hier, ax_c, ax_d, ax_e, ax_f, ax_g)


# ============================================================
# 5. Save figure and statistics
# ============================================================
fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight', transparent=False)
os.makedirs(STATS_DIR, exist_ok=True)
pd.DataFrame(STATS_ROWS).to_csv(STATS_CSV, index=False)
pd.DataFrame(PANEL_B_EDGE_ROWS).to_csv(
    os.path.join(STATS_DIR, "figure9_panel_b_root_area_fc_te_edges.csv"),
    index=False,
)
PANEL_B_NODE_ROWS.to_csv(
    os.path.join(STATS_DIR, "figure9_panel_b_root_area_nodes.csv"),
    index=False,
)
print(f"Saved {STATS_CSV}")
