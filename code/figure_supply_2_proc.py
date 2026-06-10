import numpy as np
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community import louvain_communities
from networkx.algorithms.community.quality import modularity
import itertools
from pathlib import Path
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from figure_style import set_paper_style, add_panel_label_fig,add_panel_label,plot_division_box_with_stats,darw_region_bar, region as style_regions
set_paper_style()

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
OUTPUT_DIR = PROJECT_ROOT / "figures"
DATA_DIR = PROJECT_ROOT / "data"


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


def load_log_outin_degree_by_region():
    region_fcv_path = first_existing_path(
        DATA_DIR / "fig1_prism_D_FCS_FCV_bar.csv",
        "./fig1_prism_D_FCS_FCV_bar.csv",
    )
    degree_fcv_path = first_existing_path(
        DATA_DIR / "fig4_prism_C_degree_FCV.csv",
        "./fig4_prism_C_degree_FCV.csv",
    )

    region_fcv_df = pd.read_csv(region_fcv_path)
    degree_fcv_df = pd.read_csv(degree_fcv_path)

    fcv_to_region = {
        round(float(fcv), 8): region
        for region, fcv in zip(region_fcv_df["Region"], region_fcv_df["FCV"])
    }

    allowed_regions = set(region_fcv_df["Region"].tolist())
    degree_by_region = {}
    for _, row in degree_fcv_df.iterrows():
        fcv_key = round(float(row["FCV"]), 8)
        region = fcv_to_region.get(fcv_key)
        if region is None:
            continue
        degree_by_region[region] = float(row["log10_OutIn_degree"])

    return degree_by_region, allowed_regions


def load_subject_log_outin_degree_data(allowed_regions):
    region_sc_path = first_existing_path(
        DATA_DIR / "figure1_emp_variation" / "region_sc.npy",
        "./region_sc.npy",
    )
    region_sc = np.load(region_sc_path)
    out_degree = np.sum(region_sc, axis=2)
    in_degree = np.sum(region_sc, axis=1)

    log_outin_degree_data = [[] for _ in range(len(style_regions))]
    log_outin_degree_sel_n = []

    for i, region_name in enumerate(style_regions):
        if region_name not in allowed_regions:
            continue
        values = np.log10((out_degree[:, i] + 1.0) / (in_degree[:, i] + 1.0))
        values = values[np.isfinite(values)]
        if len(values) > 0:
            log_outin_degree_data[i] = values.tolist()
            log_outin_degree_sel_n.append(i)

    return log_outin_degree_data, log_outin_degree_sel_n

brain_division_list = [
    0,  # 1  Medial Octavolateral Nucleus
    0,  # 2  Cerebellum
    0,  # 3  Medulla Oblongata strip1
    0,  # 4  Medulla Oblongata strip2
    0,  # 5  Medulla Oblongata strip3
    0,  # 6  Medulla Oblongata strip4
    0,  # 7  Medulla Oblongata strip5
    0,  # 8  Interpeduncular Nucleus
    0,  # 9  Inferior Olive
    1,  # 10 Caudal Hypothalamus
    0,  # 11 Raphe Nucleus
    0,  # 12 Tegmentum
    0,  # 13 Anterior Reticular Formation
    0,  # 14 Intermediate Reticular Formation
    0,  # 15 Posterior Reticular Formation
    4,  # 16 Glossopharyngeal Ganglion
    1,  # 17 Habenula
    1,  # 18 Intermediate Hypothalamus
    1,  # 19 Rostral Hypothalamus
    4,  # 20 Octaval Ganglion
    2,  # 21 Olfactory Bulb
    2,  # 22 Olfactory Epithelium
    2,  # 23 Pallium
    1,  # 24 Pituitary
    1,  # 25 Posterior Tuberculum
    2,  # 26 Preoptic Region
    1,  # 27 Pretectum
    4,  # 28 Retina
    2,  # 29 Subpallium
    3,  # 30 Tectum
    1,  # 31 Thalamus
    3,  # 32 Torus Longitudinalis
    3,  # 33 Torus Semicircularis
    4,  # 34 Trigeminal Ganglion
    0,  # 35 Vagal Region
    0   # 36 Vagus Motor Neurons
]
brain_division_list=np.array(brain_division_list)
brain_division_list=np.concatenate((brain_division_list,brain_division_list));
'''
def fun_plot_box_scatter_statcomparison(ax,x_data,  ylabel_str):
    
    import itertools
    from scipy.stats import kruskal, mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    # division 선택 (예: 0~3, Peripheral=4 제외)
    valid_divisions = [2, 1, 3,0]

    division_names = {
        0: 'Hind',
        1: 'Dien',
        2: 'Telen',
        3: 'Mesen'
    }

    # region-level mean across files

    t_region_data = x_data        # shape [n_regions]

    new_region=[];
    new_region_data=[];

    #for i in range(72):
        #if(len(t_region_data[i])>2):
            #new_region.append(i);
            #new_region_data.append(t_region_data[i]);

    new_brain_division_list=brain_division_list[new_region];
    #new_region_data=np.array(new_region_data);
    #division_data = {}

    division_data=x_data;
    
    #for div in valid_divisions:
       # idx = np.where(new_brain_division_list == div)[0]
        #division_data[div] = new_region_data[idx]
        
     

    kw_stat, kw_p = kruskal(*[division_data[d] for d in valid_divisions])
    print(f"Kruskal–Wallis {ylabel_str}: H={kw_stat:.3f}, p={kw_p:.4g}")



    pairs = list(itertools.combinations(valid_divisions, 2))

    pvals = []
    pair_labels = []

    for d1, d2 in pairs:
        u, p = mannwhitneyu(
            division_data[d1],
            division_data[d2],
            alternative='two-sided'
        )
        pvals.append(p)
        pair_labels.append((d1, d2))

    # 다중 비교 보정 (Holm)
    reject, pvals_corr, _, _ = multipletests(pvals, method='holm')

    significant_pairs = [
        (pair_labels[i][0], pair_labels[i][1], pvals_corr[i])
        for i in range(len(pairs)) if reject[i]
    ]

    print("Significant pairs (Holm corrected):")
    for d1, d2, p in significant_pairs:
        print(f"{division_names[d1]} vs {division_names[d2]} : p={p:.4g}")
        
    import seaborn as sns

    #plt.figure(figsize=(7,5))
    #plt.rcParams['font.size'] = 14
    plot_data = []
    plot_labels = []

    for d in valid_divisions:
        plot_data.extend(division_data[d])
        plot_labels.extend([division_names[d]] * len(division_data[d]))

    sns.boxplot(x=plot_labels, ax=ax,y=plot_data, linewidth=1.2,width=0.5, showfliers=False,palette=[division_colors[d] for d in valid_divisions])
    sns.stripplot(x=plot_labels, ax=ax,y=plot_data, alpha=0.7,color='black', size=5, jitter=True)

    ax.set_ylabel(ylabel_str)
    #plt.title("Empirical FC mean across brain divisions")

    # significance bar
    y_max = max(plot_data)
    h = 0.05    * (y_max - min(plot_data))
    level = 0
 
    if p < 0.001:
            p_text = "***"
    elif p < 0.01:
            p_text = "**"
    else:
            p_text = "*"
        
        # 혹은 요청하신 대로 (p=0.000) 형식 유지
        #p_label = f"{p_text}\n(p={p:.3f})" if p >= 0.001 else f"{p_text}\n(p<0.001)"
    p_label = f"{p_text}" if p >= 0.001 else f"{p_text}"
        
    for d1, d2, p in significant_pairs:
        x1 = valid_divisions.index(d1)
        x2 = valid_divisions.index(d2)
        y = y_max + h * level

        ax.plot([x1, x1, x2, x2], [y, y+h*0.2, y+h*0.2, y], lw=1.5, c='k')
        ax.text((x1+x2)/2, y+h*0.25, f"* (p={p:.3f})",
                ha='center', va='bottom', fontsize=10)
        level += 1

    #plt.tight_layout()
    #plt.savefig(filename_str,transparent=True, dpi=300)
    #plt.close();
'''
def remove_nan_from_metric(metric_data):
    """
    metric_data: list of length 72
                 each element = list of values (subjects)
    return: same structure, NaN 제거됨
    """
    clean_data = []

    for vals in metric_data:
        arr = np.asarray(vals, dtype=float)
        arr = arr[~np.isnan(arr)]
        clean_data.append(arr.tolist())

    return clean_data
'''
def add_panel_label(ax, label):
    ax.text(-0.12, 1.08, label, transform=ax.transAxes,
    fontsize=22, fontweight='bold', va='top', ha='left')
'''
# ============================================================
# Load metrics data
# ============================================================

network_metrics_path = first_existing_path(
    DATA_DIR / "sc_original_per_area_network_metrics.pkl",
)
with open(network_metrics_path, 'rb') as f:
    data = pickle.load(f)

metrics = {
    'In-degree': data['degree_data'],
    'Clustering': data['clustering_data'],
    'Betweenness': data['BC_mean_data'],
    'GlobalEfficiency': data['Eglob_data'],
    'ModularityQ': data['q_data']
}

metrics_clean = {}
for name, data_ in metrics.items():
    metrics[name] = remove_nan_from_metric(data_)

valid_divisions = [2, 1, 3, 0]


x_cluster_data=[[] for _ in range(4)];
x_between_data=[[] for _ in range(4)];
x_GlobalEfficiency_data=[[] for _ in range(4)];
x_modularityQ_data=[[] for _ in range(4)];

sel_n=[];

for i in range(72):
    d = brain_division_list[i]

    if d <= 3:
        x_cluster_data[d].extend(metrics['Clustering'][i]);
        x_between_data[d].extend(metrics['Betweenness'][i]);
        x_GlobalEfficiency_data[d].extend(metrics['GlobalEfficiency'][i]);
        x_modularityQ_data[d].extend(metrics['ModularityQ'][i]);

sel_n=[];
clus_data=[];
gl_data=[];
mq_data=[];
dac_sel_n=[];

for i in range(72):
    clus_data.append(metrics['Clustering'][i])
    gl_data.append(metrics['GlobalEfficiency'][i])
    mq_data.append(metrics['ModularityQ'][i])
    if(len(metrics['Clustering'][i])>3):
        sel_n.append(i);
    if(len(metrics['Clustering'][i])>2):
        dac_sel_n.append(i);


dac_data_path = first_existing_path(
    DATA_DIR / "total_selected_region_dac_data.npz",
)
load_data=np.load(dac_data_path);
cont_csel_id=load_data['arr_0'];
cont_new_total_dac_out_data=load_data['arr_1'];
cont_new_total_dac_in_data=load_data['arr_2'];
cont_sel_region=load_data['arr_3'];

for i in range(len(cont_new_total_dac_out_data)):
    if(len(cont_new_total_dac_out_data[i])>0):
        cont_new_total_dac_out_data[i]*=1;
        cont_new_total_dac_in_data[i]*=1;
        
         

        

# ============================================================
# Region names
# ============================================================
lregion = ['MON','Cb','MOS1','MOS2','MOS3','MOS4','MOS5','IPN','IO','Hc','Ra','T',
           'aRF','imRF','pRF','GG','Hb','Hi','HR','OG','OB','OE','P','Pi','PT',
           'PO','PrT','R','SP','TeO','Th','TL','TS','TG','VR','NX']

rregion = ['rMON','rCb','rMOS1','rMOS2','rMOS3','rMOS4','rMOS5','rIPN','rIO','rHc',
           'rRa','rT','raRF','rimRF','rpRF','rGG','rHb','rHi','rHR','rOG','rOB',
           'rOE','rP','rPi','rPT','rPO','rPrT','rR','rSP','rTeO','rTh','rTL',
           'rTS','rTG','rVR','rNX']

regions = lregion + rregion
N_REGION = 72
log_outin_degree_by_region, log_outin_allowed_regions = load_log_outin_degree_by_region()
log_outin_degree_data, log_outin_degree_sel_n = load_subject_log_outin_degree_data(
    log_outin_allowed_regions
)

# ============================================================
# Brain divisions
# ============================================================


division_names = {
    0: "Hind",
    1: "Di",
    2: "Tel",
    3: "Mes",
}

division_colors = MAIN_DIVISION_COLORS.copy()


def _region_to_division_code(region_name):
    try:
        idx = regions.index(region_name)
    except ValueError:
        return None
    d = int(brain_division_list[idx])
    return d if d <= 3 else None


def _draw_region_bar_main(ax, region_lists, selected_regions):
    means = []
    sems = []
    ordered_names = []
    ordered_divisions = []

    for idx in selected_regions:
        arr = np.asarray(region_lists[idx], dtype=float)
        arr = arr[np.isfinite(arr)]
        means.append(float(np.mean(arr)) if arr.size else np.nan)
        sems.append(float(np.std(arr, ddof=1) / np.sqrt(arr.size)) if arr.size > 1 else 0.0)
        ordered_names.append(regions[idx])
        ordered_divisions.append(int(brain_division_list[idx]))

    x = 0
    xticks = []
    xticklabels = []
    xtick_divisions = []
    last_division = None
    for d in MAIN_DIVISION_ORDER:
        for i, div in enumerate(ordered_divisions):
            if div != d:
                continue
            if last_division is not None and div != last_division:
                ax.axvline(x - 0.5, color="#B5B5B5", lw=0.8, zorder=0)
            ax.bar(
                x,
                means[i],
                yerr=sems[i],
                color=MAIN_DIVISION_COLORS[d],
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
            xtick_divisions.append(div)
            last_division = div
            x += 1

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=60, ha="right", fontsize=MAIN_TICK_FS)
    for lbl, d in zip(ax.get_xticklabels(), xtick_divisions):
        lbl.set_color(MAIN_DIVISION_COLORS.get(d, "#333333"))
        if d == 2:
            lbl.set_fontweight("bold")
            lbl.set_fontstyle("italic")
    ax.tick_params(axis="both", which="both", direction="out", bottom=True, left=True, length=3, width=0.8, labelsize=MAIN_TICK_FS)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(MAIN_AXIS_LW)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _display_region_name(region_name):
    region_name = str(region_name)
    if len(region_name) > 1 and region_name[0] == "l" and region_name[1].isupper():
        return region_name[1:]
    return region_name


def _zebrafish_division_label(node):
    base = _display_region_name(node).upper()
    if base in {"P", "SP", "OB", "OG", "OE", "PO"}:
        return "Tel"
    if base in {"HB", "HI", "HR", "TH", "PT", "PRT"}:
        return "Di"
    if base in {"TEO", "TL", "TS"}:
        return "Mes"
    return "Hind"


def _load_figure9_region_order():
    fig9_summary_path = first_existing_path(
        DATA_DIR / "final_summary_tables" / "figure1_dynamic_fc_fingerprint_overview_values.csv",
    )
    fig9_summary = pd.read_csv(fig9_summary_path)
    fig9_summary = fig9_summary.loc[fig9_summary["species"].eq("Zebrafish")].copy()
    fig9_summary["display_node"] = fig9_summary["node"].astype(str).map(_display_region_name)
    fig9_summary["division"] = fig9_summary["display_node"].map(_zebrafish_division_label)
    fig9_summary["_division_order"] = fig9_summary["division"].map({"Tel": 0, "Di": 1, "Mes": 2, "Hind": 3})
    fig9_summary = (
        fig9_summary.dropna(subset=["_division_order"])
        .sort_values(["_division_order", "display_node"])
        .drop_duplicates("display_node")
        .reset_index(drop=True)
    )
    selected = []
    for node in fig9_summary["display_node"].astype(str):
        if node in regions:
            selected.append(regions.index(node))
    return selected


def _metric_lists_from_table(df, value_col, selected_regions):
    region_lists = [[] for _ in range(N_REGION)]
    for region_idx in selected_regions:
        region_name = regions[region_idx]
        values = df.loc[df["node"].astype(str).eq(region_name), value_col].to_numpy(float)
        values = values[np.isfinite(values)]
        if values.size:
            region_lists[region_idx] = values.tolist()
    return region_lists, selected_regions


# ============================================================
# File paths for network diagrams
# ============================================================
save_file_patha_synapse_sc_aux_data = []


# ============================================================
# Create combined figure
# ============================================================

sc_values_path = first_existing_path(
    DATA_DIR / "final_summary_tables" / "sc_four_measures_vs_fcv_all_species_values.csv",
)
oo_values_path = first_existing_path(
    DATA_DIR / "final_summary_tables" / "oo_fraction_recomputed_values_by_species.csv",
)
sc_values = pd.read_csv(sc_values_path)
sc_values = sc_values.loc[sc_values["species"].eq("Zebrafish")].replace([np.inf, -np.inf], np.nan)
oo_values = pd.read_csv(oo_values_path)
oo_values = oo_values.loc[oo_values["species"].eq("Zebrafish")].replace([np.inf, -np.inf], np.nan)
selected_region_order = _load_figure9_region_order()

feature_panels = [
    ("A", r"$\mathrm{DCA}_{\mathrm{post}}$", *_metric_lists_from_table(sc_values, "PostDCA", selected_region_order)),
    ("B", r"$\mathrm{DCA}_{\mathrm{pre}}$", *_metric_lists_from_table(sc_values, "PreDCA", selected_region_order)),
    ("C", "Modularity Q", *_metric_lists_from_table(sc_values, "Modularity", selected_region_order)),
    ("D", r"$\log_{10}$(Out/In)", *_metric_lists_from_table(sc_values, "LogOutIn", selected_region_order)),
    ("E", "OO fraction", *_metric_lists_from_table(oo_values, "OO_fraction", selected_region_order)),
]

fig = plt.figure(figsize=(16, 15))
axes = [
    plt.subplot2grid((5, 6), (row, 0), colspan=6)
    for row in range(5)
]

for ax, (panel_label, ylabel, region_lists, selected_regions) in zip(axes, feature_panels):
    _draw_region_bar_main(ax, region_lists, selected_regions)
    ax.set_ylabel(ylabel, fontsize=MAIN_AXIS_FS)
    ax.set_xlim(-1.1, len(selected_regions))
    ax.yaxis.set_label_coords(-0.055, 0.5)
    add_panel_label_fig(fig, ax, panel_label, dx=-0.075, dy=0.008, fontsize=MAIN_PANEL_FS)

axes[0].set_yticks([-0.04, -0.02, 0.00, 0.02])
axes[1].set_yticks([-0.04, -0.02, 0.00, 0.02])

for ax in axes:
    pos = ax.get_position()
    ax.set_position([pos.x0 + 0.02, pos.y0, pos.width * 0.95, pos.height * 0.90])




OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT_DIR / 'figure_supply_2_proc.png', dpi=600, bbox_inches='tight')
