# Figure 12 / Main Fig. 2

## Purpose

This figure summarizes zebrafish directed structural-connectivity organization.
It combines a region-level SC feature heatmap, a directed region-level SC
network, representative pair diagrams, and division-level structural summaries.

## Input Data

- `data/final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`:
  zebrafish DCA, modularity, and log out/in values.
- `data/final_summary_tables/oo_fraction_recomputed_values_by_species.csv`:
  zebrafish output-output motif fraction values.
- `data/zebrafish_heatmap_matched_region_sc_network_data.npz`: cached
  region-level directed SC network used for Panel B.
- `data/total_selected_region_dac_data.npz`: region-level DCA axis used for
  Panel B node placement.
- `data/network_*_network_diagarm.png`: representative pair diagrams used for
  Panel C.

## Calculation

Panel A builds a structural feature vector for each region containing
post-DCA, pre-DCA, modularity Q, log out/in ratio, and OO fraction. Values are
averaged within region, each feature row is z-scored across regions, and
columns are ordered by Ward linkage. `rOB` is retained where SC values are
available; its missing OO fraction is filled only for this main figure by the
telencephalic division mean so that the node can be displayed with complete
Panel A features.

Panel B shows the cached directed SC network arranged by anatomical division
and post-DCA axis. Panel C shows representative local network diagrams.

Panels D-H summarize division-level structural values:

- D: post-DCA
- E: pre-DCA
- F: modularity Q
- G: log out/in ratio
- H: OO fraction

Division labels are `Tel`, `Di`, `Mes`, and `Hind`.  Pairwise comparisons use
two-sided Mann-Whitney U tests with Holm correction.  These panels are
descriptive summaries of the precomputed structural metrics, not recomputation
of the original graph construction.

## Output

- `figures/figure12_final.png`
- `outputs/figure12_final.png`
- `data/final_summary_tables/figure12_stats.csv`: pairwise division Mann-Whitney U tests
  with uncorrected and Holm-corrected p values for panels D-H.

Run from the root of `final_figure_pack/`:

```bash
python3.10 code/figure12_clean.py
```
