# Sensitivity And Control Analysis Provenance

This note records where the manuscript/SI sensitivity and control values are stored in the public figure pack.

## Zebrafish SC Threshold Sensitivity

Purpose: verify that the main zebrafish regional `DCA_post`--FCV association was not driven by one narrow endpoint-distance threshold.

Included tables:

- `data/final_summary_tables/zebrafish_region_sc_threshold_dca_fcv_correlations.csv`
- `data/final_summary_tables/zebrafish_region_sc_threshold_summary.csv`
- `data/final_summary_tables/zebrafish_global_sc_post_dca_region_mean_47_correlations.csv`
- `data/final_summary_tables/zebrafish_all_candidate_predictors_fcv_correlations_47.csv`
- `data/final_summary_tables/zebrafish_all_candidate_predictors_fcv_correlations_47_fcv_aligned.csv`
- `docs/zebrafish_global_sc_post_dca_region_mean_47_summary.md`

Upstream source:

- `fcv_postdca_raw_recompute/out_data/zebrafish/region_sc_post_dca/`
- `fcv_postdca_raw_recompute/zebrafish_dca_intra_area_validation/out/`

Current final-pack values:

- In the subject-averaged 47-region validation table, `IntraPostDCA` was positively associated with `ZFCV` (Pearson `r = 0.5868`, `p = 1.46e-05`) and `RawFCV` (Pearson `r = 0.5928`, `p = 1.13e-05`).
- The `region_sc_threshold_*` tables are retained as older thresholded regional-SC sensitivity outputs, but their `RegionSC_PostDCA_scaled_sqrtN` rows do not correspond to the final SI sentence reporting the stronger `DCA_post`--FCV range. Use the rank1/global validation tables above for the final DCA/FCV provenance.

## SC-Based FCV Prediction Residual Controls

Purpose: test whether residual errors from the SC-based FCV prediction were explained by simple measurement or sampling proxies.

Included script and tables:

- `code/figure_sc_fc_fcv_prediction_robustness.py`
- `data/final_summary_tables/figure_sc_fc_fcv_prediction_robustness_summary.csv`
- `data/final_summary_tables/figure_sc_fc_fcv_residual_proxy_summary.csv`

Manuscript/SI summary:

- Absolute prediction residuals were not significantly associated with FCV SEM, calcium-signal variance SEM, or log10 neuron count in the included table.
- Current table values are `p = 0.389` for FCV SEM, `p = 0.581` for calcium-signal variance SEM, and `p = 0.399` for log10 neuron count when using absolute residuals. If manuscript/SI text still reports older values, update it to match this table or restore the older generating table.

## Sliding-Window Sensitivity Of FCV Estimates

Purpose: test whether FCV rankings were preserved across reasonable changes in sliding-window length.

Included script and tables:

- `code/check_current_base_window_sensitivity_highpass.py`
- `data/final_summary_tables/current_base_window_sensitivity_baseline_stability.csv`
- `data/final_summary_tables/current_base_window_sensitivity_node_summary.csv`
- `data/final_summary_tables/current_base_window_sensitivity_qc.csv`
- `data/final_summary_tables/current_base_window_sensitivity_vs_fcv_correlations.csv`

Manuscript/SI summary:

- Zebrafish FCV was stable across short/current/long windows.
- Drosophila FCV was stable across 15/5-, 20/10-, and 30/8-frame window/step settings.
- C. elegans FCV remained positively correlated across the current and longer-window summaries.

## Figure-Level Global Statistics

The main caption-level Kruskal-Wallis and effect-size statistics are stored explicitly in:

- `data/final_summary_tables/figure9_global_stats.csv`
- `data/final_summary_tables/figure12_global_stats.csv`

Pairwise Holm-corrected Mann-Whitney U results remain in:

- `data/final_summary_tables/figure9_stats.csv`
- `data/final_summary_tables/figure12_stats.csv`
