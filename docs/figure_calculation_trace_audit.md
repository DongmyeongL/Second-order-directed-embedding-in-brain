# Figure calculation trace audit

Audit date: 2026-06-13  
Scope: figure files in `final_figure_pack/outputs/`, compared with the authoritative manuscript files in `../real_manuscript/Research_report.tex` and `../real_manuscript/SI_Appendix.tex`.

This audit tracks where each submitted figure is generated, what compact data it reads, what calculations are performed inside the plotting code, and which outputs/statistics are written.

## High-level result

The package contains the PNGs needed by the current manuscript, but the automatic final-pack workflow is partly stale:

- Current manuscript uses `figure7_14_combined_group_panels_AB.png` and `figure7_14_combined_group_panels_CDE.png`.
- `scripts/run_all_figures.sh`, `scripts/sync_outputs.sh`, `README.md`, `outputs/README.md`, and `MANIFEST.md` still treat `figure14_15_celegans_drosophila_AFHI_combined.png` and `figure7_by_species_with_figure5_abc.png` as the final main Fig. 5/6 outputs.
- Therefore, the files are present, but the documented regeneration path does not yet match the current TeX.

## Main figures

### Fig. 1: `figure9_final.png`

TeX reference:

- `Research_report.tex:91`, label `fig:fc_analysis`

Primary code:

- `code/figure9_clean.py`

Main inputs:

- `data/final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`
- `data/source_inputs/ncomms_tables/highpass_ce_zf_plot_measures_recording_node.csv`
- `data/final_summary_tables/observed_nette_no_p_recording_level.csv`
- `data/region_community_io/subject_*/subject_*_causality.npz`

Main calculations performed in code:

- Loads regional FC measures and TE-derived directed-drive summaries.
- Uses within-recording normalized values for division-level panels.
- Computes Kruskal--Wallis and Holm-corrected Mann--Whitney U rows for panels C--G.
- Builds the filtered regional FC-supported TE network for panel B.

Outputs written:

- `figures/figure9_final.png`
- `data/final_summary_tables/figure9_stats.csv`
- `data/final_summary_tables/figure9_panel_b_root_area_fc_te_edges.csv`

Additional caption table:

- `data/final_summary_tables/figure9_global_stats.csv`

Audit status:

- Figure, code, and global statistics are aligned with the current manuscript.
- The main input FC/TE values are compact derived tables, not raw calcium traces. Raw TE computation provenance should remain in Methods/SI.

### Fig. 2: `figure12_final.png`

TeX reference:

- `Research_report.tex:127`, label `fig:sc_analysis`

Primary code:

- `code/figure12_clean.py`

Main inputs:

- `data/final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`
- `data/final_summary_tables/oo_fraction_recomputed_values_by_species.csv`
- `data/total_selected_region_dac_data.npz`
- `data/zebrafish_heatmap_matched_region_sc_network_data.npz`

Main calculations performed in code:

- Loads regional SC feature values for zebrafish.
- Builds panel A feature matrix.
- Loads/caches directed region-level SC network for panel B.
- Draws representative neuron-level network diagrams for panel C.
- Computes Kruskal--Wallis and Holm-corrected Mann--Whitney U rows for panels D--H.

Outputs written:

- `figures/figure12_final.png`
- `data/final_summary_tables/figure12_stats.csv`

Additional caption table:

- `data/final_summary_tables/figure12_global_stats.csv`

Audit status:

- Figure and statistics are aligned with the current manuscript.
- Methods should explicitly state that regional DCA implementation averages DCA over inter-regional edge endpoints, not over only the unique set of source/target neurons.

### Fig. 3: `figure_sc_fc_final_overview.png`

TeX reference:

- `Research_report.tex:185`, label `fig:dca_vs_fcv`

Primary code:

- `code/figure_sc_fc_final_overview.py`

Main inputs:

- `data/final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`
- Indirectly calls `figure12_clean._load_panel_a_feature_matrix()` to obtain the SC feature matrix from the Fig. 2 source tables.

Main calculations performed in code:

- Matches FC and SC regional features.
- Computes all 25 FC--SC Pearson correlations.
- Applies Benjamini--Hochberg FDR correction.
- Fits PLSCanonical LV1 between five FC features and five SC features.
- Computes LV1 score correlation.
- Computes leave-one-region-out PLS weight stability.
- Fits structural-feature linear prediction of FCV with shuffled fivefold cross-validation.
- Computes observed--predicted FCV Pearson r, cross-validated R2, standardized beta coefficients, and leave-one-region-out beta stability.

Outputs written:

- `figures/figure_sc_fc_final_overview.png`
- `data/final_summary_tables/figure_sc_fc_final_overview_stats.csv`
- `data/final_summary_tables/figure_sc_fc_final_overview_pairwise_correlations.csv`
- `data/final_summary_tables/figure_sc_fc_final_overview_weights.csv`
- `data/final_summary_tables/figure_sc_fc_final_overview_fcv_predictions.csv`

Key computed values found in output tables:

- PLS LV1 FC-score vs SC-score: `r=0.685201`, `p=2.08e-07`, `n=45`.
- FCV vs DCA_post: `r=0.576987`, BH-FDR `q=0.000799`.
- FCV vs OO fraction: `r=0.559753`, BH-FDR `q=0.000799`.
- FCV prediction: observed vs predicted `r=0.440760`, `p=0.002443`, `R2_CV=0.134318`.

Audit status:

- Code directly recomputes the correlations, PLS, and cross-validation used in the main text.
- The prediction result is modest but matches the manuscript wording.

### Fig. 4: `figure13_final.png`

TeX reference:

- `Research_report.tex:225`, label `fig:layer-linear-model`

Primary code:

- `code/figure13_clean.py`

Related model scripts:

- `code/analyze_layer_asymmetric_epsilon_linear_model.py`
- `code/plot_layer_asymmetric_epsilon_linear_model.py`

Main inputs:

- `data/figure13/tfigure5_compare_xy_scatter_data.pkl`
- `data/figure13/wsim_fc_null_p_sp_out_fc_mean_std_data.pkl`
- `data/figure13/wsim_fc_null_p_sp_in_fc_mean_std_data.pkl`
- `data/figure13/figure1_whole_brain_fc_mean_std_data.npz`
- `data/final_summary_tables/figure13_layer_energy_potential_curves.csv`
- `data/final_summary_tables/figure13_layer_energy_potential_summary.csv`

Main calculations performed in code:

- Uses cached model and null-perturbation outputs for Base, NULL-In, and NULL-Out.
- Plots layer-wise potential landscapes from cached layer model summaries.
- Computes bootstrap mean differences for FCV/FCS comparisons.
- Writes Kruskal--Wallis and bootstrap contrast rows.

Outputs written:

- `figures/figure13_final.png`
- `data/final_summary_tables/figure13_stats.csv`

Key computed values found in output tables:

- FCV NULL-Out minus Base: `-0.558017`, 95% CI `[-0.746989, -0.364615]`, Holm `p=0`.
- FCV NULL-In minus Base: `-0.315620`, 95% CI `[-0.509803, -0.122503]`, Holm `p=0.0044`.
- FCV NULL-Out minus NULL-In: `-0.242396`, 95% CI `[-0.452665, -0.027238]`, Holm `p=0.0286`.

Audit status:

- Main numerical results match the manuscript.
- The code plots from cached simulation outputs; the full simulation-generation path is represented by separate scripts/data and should be documented as cached model outputs, not raw recomputation from scratch during plotting.

### Fig. 5: `figure7_14_combined_group_panels_AB.png`

TeX reference:

- `Research_report.tex:266`, label `fig:celegans_drosophila_anatomical`

Primary code:

- `code/figure7_14_combined_group_panels_AB.py`

Indirect calculation/plotting modules:

- `code/figure7_14_combined_group_panels.py`
- `code/figure14_15_celegans_drosophila_AFHI_combined.py`

Main inputs:

- Cross-species FCV and DCA summary tables loaded by `figure14_15_celegans_drosophila_AFHI_combined.py`.
- C. elegans and Drosophila processed data under:
  - `data/figure14_celegans/`
  - `data/figure15_drosophila/`
  - `data/final_summary_tables/figure7_species_recording_group_holm_tests*.csv`

Main calculations performed in code:

- AB script is a wrapper: it draws selected anatomical FCV/DCA group panels using functions from `figure14_15_celegans_drosophila_AFHI_combined.py`.
- FCV group bars use Mann--Whitney U tests with Holm correction.
- DCA panels use group-label permutation tests with Holm correction.

Outputs written:

- `figures/figure7_14_combined_group_panels_AB.png`
- `outputs/figure7_14_combined_group_panels_AB.png`

Audit status:

- The file is present and matches the current manuscript.
- The script is not included in `scripts/run_all_figures.sh`.
- The file is not copied by `scripts/sync_outputs.sh`, although the script writes directly to `outputs/`.
- README/MANIFEST still document the old combined cross-species figure instead of this TeX-used split figure.

### Fig. 6: `figure7_14_combined_group_panels_CDE.png`

TeX reference:

- `Research_report.tex:287`, label `fig:functional_group_synthesis`

Primary code:

- `code/figure7_14_combined_group_panels_CDE.py`

Indirect calculation/plotting modules:

- `code/figure7_14_combined_group_panels.py`
- `code/figure7_by_species_with_figure5_abc.py`
- `code/figure7_structure_function_synthesis_effects.py`

Main inputs:

- `data/final_summary_tables/figure7_by_species_with_figure5_abc_effects.csv`
- `data/final_summary_tables/figure7_by_species_with_figure5_abc_group_signature.csv`
- `data/final_summary_tables/figure7_by_species_with_figure5_abc_group_permutation_tests.csv`
- `data/source_inputs/ncomms_tables/figure7_recording_functional_group_points.csv`
- `data/source_inputs/ncomms_tables/figure7_postdca_functional_group_node_values.csv`
- `data/source_inputs/ncomms_tables/figure7_predca_functional_group_node_values.csv`

Main calculations performed in code:

- CDE script is a wrapper around the functional-group drawing functions.
- `figure7_by_species_with_figure5_abc.py` can recompute statistics with `--recompute-stats`, but default execution reads cached CSVs.
- Functional-group values are centered/scaled within species and metric.
- Group enrichment uses group-vs-rest standardized mean difference and label permutation preserving group sizes, with Holm correction.

Outputs written:

- `figures/figure7_14_combined_group_panels_CDE.png`
- `outputs/figure7_14_combined_group_panels_CDE.png`

Key computed values found in output tables:

- C. elegans olf./chemo FCV enrichment: observed z-delta `0.860393`, Holm `p=0.004995`.
- Drosophila olf./chemo FCV enrichment: observed z-delta `0.801396`, Holm `p=0.004995`.
- Tables contain the values quoted in the current comparative Results.

Audit status:

- The file is present and matches the current manuscript.
- The script is not included in `scripts/run_all_figures.sh`.
- The file is not copied by `scripts/sync_outputs.sh`, although the script writes directly to `outputs/`.
- README/MANIFEST still document the old functional-group figure instead of this TeX-used split figure.

## Supplementary figures

### Fig. S1: `figure_supply_1.png`

TeX reference:

- `SI_Appendix.tex:506`, label `fig:SI_synapse`

Primary code:

- `code/figure_supply_1.py`

Main inputs:

- `data/figure9_distance_synapse/figure_supply1_fc_dist_data.npz`
- `data/figure9_distance_synapse/figure_supply1_soma_dist_data.npy`

Audit status:

- Figure, code, and SI reference are aligned.

### Fig. S2: `figure_supply_2_proc.png`

TeX reference:

- `SI_Appendix.tex:517`, label `fig:SI_struct_metrics`

Primary code:

- `code/figure_supply_2_proc.py`

Main inputs:

- Regional FCV/degree and SC metric files.
- `data/sc_original_per_area_network_metrics.pkl`
- `data/total_selected_region_dac_data.npz`
- `data/final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`
- `data/final_summary_tables/oo_fraction_recomputed_values_by_species.csv`

Audit status:

- Figure, code, and SI reference are aligned. `rOB` is included in the displayed
  region order after `OB`; structural metrics are shown where available, while
  absent values such as the `rOB` OO fraction remain blank rather than being
  imputed.

### Fig. S3: `figure_supply_10_proc.png`

TeX reference:

- `SI_Appendix.tex:532`, label `fig:SI_fc_metrics`

Primary code:

- `code/figure_supply_10_proc.py`

Main inputs:

- `data/final_summary_tables/figure1_dynamic_fc_fingerprint_overview_values.csv`
- `data/region_community_io/subject_*/subject_*_causality.npz`
- `data/region_community_io/subject_*/subject_*_net_te_drive_fc_neighbors.npz`

Audit status:

- Figure, code, and SI reference are aligned. `rOB` is included in the displayed
  region order after `OB`; FCV/FCS/FC-partner values are shown where available,
  while missing TE and neighbor-TE values remain blank.

### Fig. S4: `figure_supply_13.png`

TeX reference:

- `SI_Appendix.tex:548`, label `fig:SI_TE`

Primary code:

- `code/figure_supply_13.py`

Main inputs:

- `data/region_community_io/subject_13/subject_13_causality.npz`
- `data/figure_supply_13/figure_supply_13_example.npz`

Audit status:

- Figure, code, and SI reference are aligned.

### Fig. S5: `figure_supply_14.png`

TeX reference:

- `SI_Appendix.tex:565`, label `fig:SI_stimulus_fcv`

Primary code:

- `code/figure_supply_14.py`

Main inputs:

- `data/final_summary_tables/figure7_zebrafish_all_region_stimulus_values.csv`

Main calculations performed in code:

- Loads precomputed C. elegans/zebrafish spontaneous-vs-stimulus summary values from the Fig. 7 tables.
- Draws panels A--C for spontaneous/stimulus FCV relationships.
- Computes region x stimulus mean heatmap values for zebrafish OMR stimuli.
- Writes the plotted heatmap value table.

Outputs written:

- `output/png/figure_supply_14.png`
- `figures/figure_supply_14.png`
- `outputs/supplementary/figure_supply_14.png`
- `output/stats/figure_supply_14_stimulus_fcv_heatmap_values.csv`

Audit status:

- Figure, code, and SI reference are aligned at the file level.
- Methods text should be checked against the actual stimulus-FCV window implementation. Earlier audit found the SI method says 5 s / 2.5 s, while the code uses 20-frame windows and 5-frame steps for stimulus FCV.

### Fig. S6: `figure_supply_0.png`

TeX reference:

- `SI_Appendix.tex:584`, label `fig:SI_sim`

Primary code:

- `code/figure_supply_0.py`

Main inputs:

- `data/figure1_emp_variation/figure1_whole_brain_fc_mean_std_data.npz`
- cached simulation pickle/NPZ files under `data/figure5_simulation/`

Audit status:

- Figure, code, and SI reference are aligned.

### Fig. S7: `figure_supply_5.png`

TeX reference:

- `SI_Appendix.tex:604`, label `fig:SI_null_dca`

Primary code:

- `code/figure_supply_5.py`

Main inputs:

- `data/figure6_NULL_P_SP/figure6_network_properites_in.npz`
- `data/figure6_NULL_P_SP/figure6_network_properites_out.npz`
- network diagram PNGs under `data/network_diagrams/`

Audit status:

- Figure, code, and SI reference are aligned.

### Fig. S8: `figure_supply_15.png`

TeX reference:

- `SI_Appendix.tex:619`, label `fig:SI_traces`

Primary code:

- `code/figure_supply_15.py`

Main inputs:

- `data/figure_supply_15_rsp_rmos5_trace.npz`
- `data/figure_supply_15_layer_trace_cache.npz`
- selected cached/null simulation files.

Audit status:

- Figure, code, and SI reference are aligned.

## Calculation trace concerns

### 1. Split cross-species figures are not in the automated run list

The current manuscript uses AB/CDE split figures. Their scripts exist and save correct PNGs, but `run_all_figures.sh` does not call them. A fresh user running the documented all-figure command would regenerate old cross-species figures but not the two TeX-used split figures.

Recommended change:

```bash
python3.10 code/figure7_14_combined_group_panels_AB.py
python3.10 code/figure7_14_combined_group_panels_CDE.py
```

should be added to `scripts/run_all_figures.sh`.

### 2. `sync_outputs.sh` does not sync AB/CDE

Although AB/CDE scripts save directly to `outputs/`, the sync script should still copy them from `figures/` to make the package behavior consistent.

Recommended change:

```bash
cp figures/figure7_14_combined_group_panels_AB.png outputs/
cp figures/figure7_14_combined_group_panels_CDE.png outputs/
```

### 3. README/MANIFEST point to legacy main Fig. 5/6

Current package documentation says:

- Fig. 5 = `figure14_15_celegans_drosophila_AFHI_combined.png`
- Fig. 6 = `figure7_by_species_with_figure5_abc.png`

Current manuscript says:

- Fig. 5 = `figure7_14_combined_group_panels_AB.png`
- Fig. 6 = `figure7_14_combined_group_panels_CDE.png`

This documentation mismatch should be fixed before upload.

### 4. Some figures plot cached/derived values rather than recomputing raw analysis

This is acceptable for a figure reproduction pack, but it should be stated clearly:

- Fig. 1 uses compact FC/TE summary tables and subject TE NPZ outputs.
- Fig. 2 uses compact SC/DCA/OO summary tables and network caches.
- Fig. 4 uses cached simulation/null outputs.
- Fig. 5/6 use processed cross-species FCV/DCA tables and cached permutation summaries unless recomputation is requested.

The pack therefore reproduces final figures from processed analysis outputs, not from all raw calcium/synapse datasets.

### 5. Known manuscript-method mismatches still need text correction

These are not caused by the PNG files themselves, but they affect whether a reader can reproduce the stated method:

- Stimulus FCV window wording in SI does not match `figure_supply_14.py`.
- Residual-control p-values in SI do not match `figure_sc_fc_fcv_residual_proxy_summary.csv`.
- Threshold-sensitivity statement needs exact provenance wording because the final-pack threshold table does not directly reproduce the full `r=0.580--0.650` range.
- DCA methods should specify edge-endpoint averaging.

## Recommended immediate edits

1. Patch `run_all_figures.sh`, `sync_outputs.sh`, `README.md`, `outputs/README.md`, and `MANIFEST.md` to match the current real manuscript figures.
2. Keep old `figure14_15...png`, `figure7_by_species...png`, and `figure7_14_combined_group_panels.png` either in a `legacy/` folder or remove them from `outputs/` to avoid confusion.
3. Fix the SI method/statistic mismatches listed above.
4. Re-run the final-pack workflow and verify that the TeX-used PNGs are regenerated from code.
