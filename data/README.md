# Data Guide

This folder contains processed, figure-ready inputs for the manuscript figures.
It is not a raw-data archive. Full raw calcium imaging, anatomical
reconstruction, and external connectome datasets should be obtained from the
original sources cited in the manuscript.

## Species-Level Measure Tables

The main cross-species FC and SC summaries are stored as CSV files in
`final_summary_tables/`.

| Species | FC measure data | SC / DCA measure data | Notes |
|---|---|---|---|
| Zebrafish | `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`; recording-level values in `final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`; TE inputs in `final_summary_tables/observed_nette_no_p_recording_level.csv` and `region_community_io/subject_*/` | `final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`; `final_summary_tables/oo_fraction_recomputed_values_by_species.csv`; network cache in `zebrafish_heatmap_matched_region_sc_network_data.npz`; DCA arrays in `total_selected_region_dac_data.npz` | FC measures include FCV, FCS, FC partner reconfiguration, observed NetTE, and neighbor NetTE. SC measures include DCA_post, DCA_pre, modularity, log out/in, and OO fraction. |
| C. elegans | `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`; recording-level values in `final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`; heat/spontaneous FCV files in `source_inputs/external_processed/fcv_postdca_raw_recompute/out_data/celegans/geometric_fcv/` | `final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`; `final_summary_tables/oo_fraction_recomputed_values_by_species.csv`; functional-group DCA tables in `source_inputs/ncomms_tables/` | Values are neuron-level after matching WormWideWeb calcium recordings to connectome neuron labels. |
| Drosophila | `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`; recording-level values in `final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv` | `final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`; `final_summary_tables/oo_fraction_recomputed_values_by_species.csv`; functional-group DCA tables in `source_inputs/ncomms_tables/` | Values are matched side-aware Ito-region summaries derived from Branson calcium data and FlyWire structural annotations. |

## Core Cross-Species Tables

| File | Contents |
|---|---|
| `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv` | Node or region summaries for FCV, FCS, FC partner reconfiguration, observed NetTE, and neighbor NetTE across species. |
| `final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv` | Recording-by-node FC measures used for point distributions and group summaries. |
| `final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv` | Matched FCV and SC measures: DCA_post, DCA_pre, modularity, and log out/in. |
| `final_summary_tables/oo_fraction_recomputed_values_by_species.csv` | Output-output motif fraction and related motif counts. |
| `source_inputs/ncomms_tables/figure7_recording_functional_group_points.csv` | FCV values summarized by species, recording, and functional group. |
| `source_inputs/ncomms_tables/figure7_postdca_functional_group_node_values.csv` | DCA_post values mapped to comparative functional groups. |
| `source_inputs/ncomms_tables/figure7_predca_functional_group_node_values.csv` | DCA_pre values mapped to comparative functional groups. |
| `final_summary_tables/figure7_by_species_with_figure5_abc_group_signature.csv` | Functional-group mean signatures used in manuscript Fig. 6. |
| `final_summary_tables/figure7_by_species_with_figure5_abc_group_permutation_tests.csv` | Functional-group enrichment statistics used in manuscript Fig. 6. |

For a figure-by-figure index of statistics and plotted-value tables, see
`statistics_by_figure/README.md`. The machine-readable version is
`statistics_by_figure/statistics_index.csv`.

## Figure-Specific Inputs

| Figure output | Main code | Primary data inputs |
|---|---|---|
| `figure9_final.png` | `code/figure9_clean.py` | `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`; `final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`; `final_summary_tables/observed_nette_no_p_recording_level.csv`; `region_community_io/subject_*/subject_*_causality.npz` |
| `figure12_final.png` | `code/figure12_clean.py` | `final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`; `final_summary_tables/oo_fraction_recomputed_values_by_species.csv`; `total_selected_region_dac_data.npz`; `zebrafish_heatmap_matched_region_sc_network_data.npz`; helper inputs used by `code/figure_supply_sc_heatmap.py` |
| `figure_sc_fc_final_overview.png` | `code/figure_sc_fc_final_overview.py` | Zebrafish FC measures from `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`; Zebrafish SC feature matrix loaded through `code/figure12_clean.py`; output statistics in `final_summary_tables/figure_sc_fc_final_overview_*.csv` |
| `figure13_final.png` | `code/figure13_clean.py` | `figure13/layer_asymmetric_epsilon_linear_data.npz`; `figure13/figure1_whole_brain_fc_mean_std_data.npz`; `figure13/tfigure5_compare_xy_scatter_data.pkl`; `figure13/wsim_fc_null_p_sp_in_fc_mean_std_data.pkl`; `figure13/wsim_fc_null_p_sp_out_fc_mean_std_data.pkl`; `figure6_NULL_P_SP/*.npz`; `network_diagrams/*.png` |
| `figure7_14_combined_group_panels_AB.png` | `code/figure7_14_combined_group_panels_AB.py` | Cross-species anatomy panels drawn through `code/figure7_14_combined_group_panels.py` and `code/figure14_15_celegans_drosophila_AFHI_combined.py`; primary measure tables are `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`, `final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`, and `final_summary_tables/oo_fraction_recomputed_values_by_species.csv` |
| `figure7_14_combined_group_panels_CDE.png` | `code/figure7_14_combined_group_panels_CDE.py` | Functional-group FCV/DCA tables in `source_inputs/ncomms_tables/`; cached group statistics in `figure7_by_species_with_figure5_abc_*.csv`; zebrafish subject-region DCA in `source_inputs/external_processed/.../zebrafish_rank1_subject_region_post_dca.csv` |

## Supplementary Figure Inputs

| Supplementary output | Main code | Primary data inputs |
|---|---|---|
| `figure_supply_0.png` | `code/figure_supply_0.py` | `figure1_emp_variation/figure1_whole_brain_fc_mean_std_data.npz`; `figure5_simulation/figure5_compare_xy_scatter_data.pkl` |
| `figure_supply_1.png` | `code/figure_supply_1.py` | `figure9_distance_synapse/figure_supply1_fc_dist_data.npz`; `figure9_distance_synapse/figure_supply1_soma_dist_data.npy` |
| `figure_supply_2_proc.png` | `code/figure_supply_2_proc.py` | Zebrafish SC metric inputs: `sc_original_per_area_network_metrics.pkl`, `total_selected_region_dac_data.npz`, `figure1_emp_variation/region_sc.npy`, `fig4_prism_C_degree_FCV.csv`, and `sc_four_measures_vs_fcv_all_species_values.csv` |
| `figure_supply_5.png` | `code/figure_supply_5.py` | `figure6_NULL_P_SP/figure6_network_properites_in.npz`; `figure6_NULL_P_SP/figure6_network_properites_out.npz` |
| `figure_supply_10_proc.png` | `code/figure_supply_10_proc.py` | Zebrafish FC/TE inputs from `final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`, `final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`, `final_summary_tables/observed_nette_no_p_recording_level.csv`, and `region_community_io/subject_*/` |
| `figure_supply_13.png` | `code/figure_supply_13.py` | `region_community_io/subject_12/subject_12_causality.npz`; `figure_supply_13/figure_supply_13_example.npz` |
| `figure_supply_14.png` | `code/figure_supply_14.py` | OMR heatmap from `final_summary_tables/figure7_zebrafish_all_region_stimulus_values.csv`; A-C scatter inputs from `source_inputs/external_processed/fcv_postdca_raw_recompute/out_data/celegans/geometric_fcv/` and the zebrafish stimulus table via `code/figure7_by_species_with_figure5_abc.py` |
| `figure_supply_15.png` | `code/figure_supply_15.py` | `figure_supply_15_layer_trace_cache.npz`; `figure_supply_15_rsp_rmos5_trace.npz`; `figure13/layer_asymmetric_epsilon_linear_data.npz`; `figure13/tfigure5_compare_xy_scatter_data.pkl` |

## Folder Roles

| Folder or file | Purpose |
|---|---|
| `final_summary_tables/` | Core figure-ready statistics, measure tables, correlations, model outputs, and group summaries. |
| `statistics_by_figure/` | Human-readable and machine-readable index of which statistics tables correspond to each manuscript figure. |
| `source_inputs/` | Processed upstream tables copied from the analysis workspace when final figure scripts require them directly. |
| `region_community_io/` | Zebrafish subject-level community and transfer-entropy arrays. |
| `figure13/`, `figure5_simulation/`, `figure6_NULL_P_SP/` | Cached simulation and null-model outputs for manuscript Fig. 4 and related supplementary figures. |
| `figure1_emp_variation/`, `figure9_distance_synapse/`, `figure_supply_13/` | Supplementary figure input bundles. |
| `network_diagrams/` | PNG network diagrams embedded in simulation figures. |

The files included here are sufficient to regenerate the bundled figures with
`scripts/run_all_figures.sh`. They are processed derivatives, not replacements
for the raw public datasets cited in the manuscript.
