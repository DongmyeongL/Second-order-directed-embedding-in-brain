# Statistics Tables By Figure

This folder indexes the statistics, plotted-value tables, and compact cache files retained for each manuscript and supplementary figure.
Most CSV statistics tables remain in `data/final_summary_tables/`; supplementary figures often use compact `.npz`, `.npy`, or `.pkl` cache files and compute simple on-panel statistics during plotting.

## Machine-Readable Index

- `statistics_index.csv`: one row per figure/statistics-table relationship, with the figure, panel/scope, source path, recorded statistic or value type, and key columns/arrays.

## Figure-Level Guide

### Fig. 1 / figure9_final.png
- **C-G**: `data/final_summary_tables/figure9_stats.csv`
  - Pairwise division comparisons for zebrafish FC measures using Mann-Whitney U tests with Holm-corrected p-values.
  - Key columns/arrays: `panel, group_1, group_2, n_group_1, n_group_2, p_uncorrected, p_holm, reject_holm_0.05`
- **B network nodes**: `data/final_summary_tables/figure9_panel_b_root_area_nodes.csv`
  - Node-level zebrafish FCV and observed NetTE values used for the region-level directed FC network panel.
  - Key columns/arrays: `node, anatomy_group, EdgeStdFCV, ObservedNetTE`
- **B network edges**: `data/final_summary_tables/figure9_panel_b_root_area_fc_te_edges.csv`
  - Root-area directed FC/TE edge summaries used in the network panel.
  - Key columns/arrays: `source, target, source_group, target_group, mean_net_te, n_subjects_te, n_subjects_fc_sig`

### Fig. 2 / figure12_final.png
- **D-H**: `data/final_summary_tables/figure12_stats.csv`
  - Pairwise division comparisons for zebrafish SC measures using Mann-Whitney U tests with Holm-corrected p-values.
  - Key columns/arrays: `panel, group_1, group_2, n_group_1, n_group_2, p_uncorrected, p_holm, reject_holm_0.05`

### Fig. 3 / figure_sc_fc_final_overview.png
- **A**: `data/final_summary_tables/figure_sc_fc_final_overview_pairwise_correlations.csv`
  - All 25 FC-SC feature Pearson correlations with Benjamini-Hochberg FDR-corrected p-values.
  - Key columns/arrays: `fc_feature, sc_feature, pearson_r, p_value, p_fdr_bh, fdr_bh_significant, n_regions`
- **B,E**: `data/final_summary_tables/figure_sc_fc_final_overview_stats.csv`
  - PLS LV1 score correlation and fivefold cross-validated SC-to-FCV prediction performance.
  - Key columns/arrays: `panel, test, metric, r, p_value, n_regions, r2_cv`
- **C,D,F**: `data/final_summary_tables/figure_sc_fc_final_overview_weights.csv`
  - PLS feature weights, leave-one-region-out stability weights, and standardized beta coefficients for FCV prediction.
  - Key columns/arrays: `model, feature_set, feature, weight, left_out_region, left_out_division, standardized_beta`
- **E**: `data/final_summary_tables/figure_sc_fc_final_overview_fcv_predictions.csv`
  - Observed FCV, fivefold cross-validated predicted FCV, and residual by zebrafish region.
  - Key columns/arrays: `region, division, observed_fcv, predicted_fcv_5fold, residual`

### Fig. 4 / figure13_final.png
- **E**: `data/final_summary_tables/figure13_stats.csv`
  - Simulation/null-model comparisons including Kruskal-Wallis and bootstrap mean-difference tests with Holm correction where applicable.
  - Key columns/arrays: `metric, test, comparison, statistic, p_value, mean_base, mean_null_in, mean_null_out, ci_low, ci_high, p_holm`
- **B-C layer model**: `data/final_summary_tables/figure13_layer_energy_potential_summary.csv`
  - Layer-wise potential-landscape summaries across epsilon, run, and layer.
  - Key columns/arrays: `epsilon, run, layer, fc_state_mean, fc_state_std, fc_state_width_p95_p05, fcv_from_fc_state`
- **B potential curves**: `data/final_summary_tables/figure13_layer_energy_potential_curves.csv`
  - Effective-potential curve samples for dynamic FC states across epsilon, run, and layer.
  - Key columns/arrays: `epsilon, run, layer, dynamic_fc_state, potential, fc_state_mean, fc_state_std`

### Fig. 5 / figure7_14_combined_group_panels_AB.png
- **A-F source values**: `data/final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`
  - Cross-species node/region FC summaries used for anatomical FCV distributions.
  - Key columns/arrays: `species, node, EdgeStdFCV, FCS, ProfileCorrDistFCV, ObservedNetTE, NeighborNetTE, class_label`
- **B-C,E-F source values**: `data/final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`
  - Matched cross-species FCV and SC measures used for anatomical DCA distributions.
  - Key columns/arrays: `species, node, EdgeStdFCV, PostDCA, PreDCA, Modularity, LogOutIn`

### Fig. 6 / figure7_14_combined_group_panels_CDE.png
- **A-C**: `data/final_summary_tables/figure7_by_species_with_figure5_abc_group_signature.csv`
  - Functional-group means, group-vs-rest differences, relative z means, SEM, and sample sizes for FCV, DCA_post, and DCA_pre.
  - Key columns/arrays: `species, shared_fine_label, metric, mean, other_group_mean, mean_minus_other_group_mean, relative_z_mean, sem, n`
- **A-C**: `data/final_summary_tables/figure7_by_species_with_figure5_abc_group_permutation_tests.csv`
  - Functional-group enrichment permutation tests with Holm and FDR corrected p-values.
  - Key columns/arrays: `species, metric, shared_fine_label, observed_zdelta, p_abs_permutation, p_holm_within_species_metric, q_fdr_within_species_metric`
- **planned contrasts**: `data/final_summary_tables/figure7_by_species_with_figure5_abc_effects.csv`
  - Planned high-vs-low functional-group contrasts with bootstrap confidence intervals and permutation p-values.
  - Key columns/arrays: `species, metric, contrast, n_high, n_low, mean_high, mean_low, mean_diff, bootstrap_ci_low, bootstrap_ci_high, permutation_two_sided_p`

### Supplementary Fig. S14 / figure_supply_14.png
- **OMR heatmap**: `data/final_summary_tables/figure7_zebrafish_all_region_stimulus_values.csv`
  - Subject-region FCV values under OMR stimulus conditions, including within-subject z-scored FCV used for the heatmap.
  - Key columns/arrays: `Subject, StimulusIndex, Region, RegionID, Division, FCV, SubjectZFCV, normalization`

### Cross-species class mapping
- **metadata**: `data/final_summary_tables/current_plot_functional_classes.csv`
  - Node-to-functional-class metadata used to interpret cross-species grouping.
  - Key columns/arrays: `species, node, class_order, class_label, broad_class_label, fine_class, classification_note`

### Cross-species SC motifs
- **source values**: `data/final_summary_tables/oo_fraction_recomputed_values_by_species.csv`
  - OO fraction, motif counts, degree values, and matched SC/FC metadata across species.
  - Key columns/arrays: `species, node, OO_fraction, OO_count, II_fraction, out_degree, PostDCA, PreDCA, EdgeStdFCV`

### Cross-species recording-level FC
- **source values**: `data/final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`
  - Recording-by-node FC measures used by figure scripts and group summaries.
  - Key columns/arrays: `species, recording_id, node, EdgeStdFCV, ProfileCorrDistFCV, FCS, window_config, n_rois`

### Zebrafish TE source values
- **source values**: `data/final_summary_tables/observed_nette_no_p_recording_level.csv`
  - Recording-level NetTE and neighbor NetTE values by node.
  - Key columns/arrays: `species, recording_id, node, NetTE, NeighborNetTE, NetTE_z, NeighborNetTE_z, n_communities`

### Supplementary Fig. S0 / figure_supply_0.png
- **empirical/simulation FC summaries**: `data/figure1_emp_variation/figure1_whole_brain_fc_mean_std_data.npz`
  - Cached empirical whole-brain FC mean/std arrays and region labels used as source values for the supplementary empirical/simulation comparison.
  - Key columns/arrays: `NPZ arrays including new_region_array and FC summary arrays`
- **simulation traces and FC matrices**: `data/figure5_simulation/figure5_reprsesnt_data_data.pkl`
  - Cached empirical and simulated calcium traces, FC matrices, FCV summaries, and selected region/node indices used for plotting.
  - Key columns/arrays: `pickle keys include emp_ca_data, sim_ca_data, save_fc, eave_fc, estd_fc, sstd_fc, sel_region`

### Supplementary Fig. S1 / figure_supply_1.png
- **distance-FC source values**: `data/figure9_distance_synapse/figure_supply1_fc_dist_data.npz`
  - Pairwise distance and FC values used to fit binned lognormal FC distributions and an exponential distance-decay curve in the plotting script.
  - Key columns/arrays: `NPZ arrays dist_data, fc_data`
- **synapse-distance source values**: `data/figure9_distance_synapse/figure_supply1_soma_dist_data.npy`
  - Synapse-to-soma or endpoint-distance values used for the distance distribution panel.
  - Key columns/arrays: `NumPy array of distances`

### Supplementary Fig. S2 / figure_supply_2_proc.png
- **regional SC metric panels**: `data/sc_original_per_area_network_metrics.pkl`
  - Cached zebrafish regional structural network metrics used for clustering, modularity, efficiency, and related SC metric panels.
  - Key columns/arrays: `pickle dictionary of regional network metric arrays`
- **DCA and SC feature panels**: `data/total_selected_region_dac_data.npz`
  - Cached regional DCA arrays used for DCA_post/DCA_pre panels and related division-level plotting.
  - Key columns/arrays: `NPZ arrays including DCA values by region/subject`
- **SC feature source values**: `data/final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`
  - Zebrafish SC feature values used for supplementary regional SC metric panels; division-level bars compute comparisons during plotting.
  - Key columns/arrays: `species, node, PostDCA, PreDCA, Modularity, LogOutIn, recording_id, Subject`
- **OO motif source values**: `data/final_summary_tables/oo_fraction_recomputed_values_by_species.csv`
  - OO-fraction and related motif-count source values used for the supplementary structural motif panel.
  - Key columns/arrays: `species, node, OO_fraction, OO_count, II_fraction, out_degree, Division`

### Supplementary Fig. S5 / figure_supply_5.png
- **null-model DCA comparisons**: `data/figure6_NULL_P_SP/figure6_network_properites_in.npz`
  - Base and Null-Out values used for on-panel Mann-Whitney U comparison of DCA_post-related null-model distributions.
  - Key columns/arrays: `NPZ arrays x_data, y_data`
- **null-model DCA comparisons**: `data/figure6_NULL_P_SP/figure6_network_properites_out.npz`
  - Base and Null-In values used for on-panel Mann-Whitney U comparison of DCA_pre-related null-model distributions.
  - Key columns/arrays: `NPZ arrays x_data, y_data`

### Supplementary Fig. S10 / figure_supply_10_proc.png
- **regional FC/TE summaries**: `data/final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`
  - Zebrafish regional FCV, FCS, FC partner reconfiguration, NetTE, and neighbor NetTE values used for supplementary FC/TE panels.
  - Key columns/arrays: `species, node, EdgeStdFCV, FCS, ProfileCorrDistFCV, ObservedNetTE, NeighborNetTE`
- **recording-level FC/TE source values**: `data/final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`
  - Recording-by-node FC values used for supplementary regional functional metric panels.
  - Key columns/arrays: `species, recording_id, node, EdgeStdFCV, ProfileCorrDistFCV, FCS`
- **recording-level NetTE source values**: `data/final_summary_tables/observed_nette_no_p_recording_level.csv`
  - Recording-level NetTE and neighbor NetTE values used for supplementary regional directed-information panels.
  - Key columns/arrays: `species, recording_id, node, NetTE, NeighborNetTE, NetTE_z, NeighborNetTE_z`

### Supplementary Fig. S13 / figure_supply_13.png
- **example causality/TE arrays**: `data/region_community_io/subject_12/subject_12_causality.npz`
  - Subject-level community causality, FC-neighbor, and NetTE arrays used for the example TE/causality supplementary panel.
  - Key columns/arrays: `NPZ arrays including net_te_matrix, fc_neighbor_mask_fdr, region_num`
- **example extracted values**: `data/figure_supply_13/figure_supply_13_example.npz`
  - Precomputed example traces/indices used by the supplementary causality example panel.
  - Key columns/arrays: `NPZ arrays for example panel data`

### Supplementary Fig. S15 / figure_supply_15.png
- **layer-model trace cache**: `data/figure_supply_15_layer_trace_cache.npz`
  - Cached layer-wise directed linear model traces and FC/FCV summaries used for supplementary model diagnostics.
  - Key columns/arrays: `NPZ arrays for layer traces and model summaries`
- **zebrafish example trace**: `data/figure_supply_15_rsp_rmos5_trace.npz`
  - Cached zebrafish example region trace used in the supplementary simulation/trace panel.
  - Key columns/arrays: `NPZ arrays for rSP/rMOS5 example trace data`

## Notes

- `figure9_stats.csv` and `figure12_stats.csv` contain pairwise division comparisons and Holm-corrected p-values, not omnibus Kruskal-Wallis rows.
- `figure_sc_fc_final_overview_*` files contain the Fig. 3 pairwise correlations, PLS summaries, cross-validated FCV predictions, and model/PLS weights.
- Functional-group statistics for manuscript Fig. 6 are split into group signatures, permutation tests, and planned contrasts.
- Supplementary figures S0, S1, S2, S5, S10, S13, and S15 primarily rely on compact source-value/cache files; where p-values are shown, they are computed inside the plotting scripts from the listed arrays rather than stored as separate CSV tables.
