# Figure 9 / Main Fig. 1

## Purpose

This figure summarizes regional functional dynamics and directed interaction
features across zebrafish anatomical divisions. Panel A shows five regional FC
features, Panel B shows the region-level directed TE network, and Panels C-G
show division-level summaries.

## Input Data

- `data/final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`:
  compact regional FC/TE feature summary used for Panel A.
- `data/source_inputs/ncomms_tables/highpass_ce_zf_plot_measures_recording_node.csv`:
  recording-level FC feature table used to restore `rOB` where available.
- `data/region_community_io/subject_*/subject_*_causality.npz`: subject-level
  transfer-entropy matrices and region order.
- `data/region_community_io/subject_*/subject_*_net_te_drive_fc_neighbors.npz`:
  precomputed neighbor-drive summaries for each subject.

## Calculation

The script assigns each region to `Tel`, `Di`, `Mes`, or `Hind`, aligns the
regional features, z-scores each feature across regions, and clusters Panel A
using Ward linkage on the five measured features. The display order is derived
from the dendrogram leaves, with a three-cluster display adjustment used for the
Panel A heatmap order. `rOB` is restored from the recording-level FC table; TE
features missing for `rOB` are filled by the corresponding telencephalic
division mean so that the node can be displayed in the main figure.

Net TE is computed from each subject causality file as the mean outgoing net TE
for each region. Neighbor Net TE is loaded from the precomputed
`fc_neighbor_mean_drive` arrays and aligned in the same way.

Panel A displays the aligned regional matrix:

- FCV
- FCS
- FC partner reconfiguration
- Net TE
- Neighbor Net TE

Panels C-G show division-level boxplots for FCV, FCS, FC partner
reconfiguration, Net TE, and Neighbor Net TE. Pairwise division comparisons use two-sided
Mann-Whitney U tests followed by Holm correction.  Panel B uses the same
regional feature set to show the hierarchical organization/ordering used in
the figure.

## Output

- `figures/figure9_final.png`
- `outputs/figure9_final.png`
- `data/final_summary_tables/figure9_stats.csv`: pairwise division Mann-Whitney U tests
  with uncorrected and Holm-corrected p values for panels C-G.
- `data/final_summary_tables/figure9_panel_b_root_area_nodes.csv`: Panel B
  node table, including the displayed `rOB` node.

Run from the root of `final_figure_pack/`:

```bash
python3.10 code/figure9_clean.py
```
