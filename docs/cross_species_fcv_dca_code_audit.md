# Cross-Species FCV/DCA Code Audit

Date checked: 2026-06-16

Scope: code-level provenance check for the C. elegans and Drosophila FCV, DCA, SC/FC construction, and SC-FC matching used in the final cross-species figures.

## Research_report.tex Methods Consistency Check

Authoritative manuscript checked:

- `ncomms_figure_build/real_manuscript/Research_report.tex`

Overall conclusion: the main FCV/DCA/SC-FC matching logic in the Methods is consistent with the implemented analysis, but several wording-level fixes are recommended so the text precisely matches the code. These are wording/provenance issues rather than evidence that the plotted FCV/DCA values were computed incorrectly.

### Items that match the code

| Methods item | Code check | Status |
|---|---|---|
| Zebrafish FCV windowing | `ncomms_figures_rebuild/scripts/compute_highpass_ce_zf_plot_measures.py` uses 20-frame windows and 5-frame steps for zebrafish. | OK |
| Zebrafish FCV definition | Code computes sliding-window Pearson FC, edge-wise temporal SD, then node/community mean FCV and region-level summaries. | OK |
| Zebrafish within-recording normalization | Figure code and source tables use within-subject/recording z-scoring before group summaries. | OK |
| Zebrafish DCA direction | `fcv_postdca_raw_recompute/code/zebrafish/compute_rank1_post_dca.py` interprets raw edges as `(pre, post)` and computes `Rank1PostDCA` from downstream target-cell DCA over outgoing inter-region edges. | OK |
| C. elegans SC/FC matching | Neuron-level chemical SC is matched to WormWideWeb activity by cleaned neuron names. FCV is z-scored within recording and averaged across recordings. | OK |
| C. elegans DCA | Shared DCA helper computes neuron-level DCA, outgoing target-weighted PostDCA, and incoming source-weighted PreDCA. | OK |
| Drosophila FCV | Final Drosophila outputs use side-aware Ito-region traces, detrending/high-pass filtering, 15-frame windows, 5-frame steps, within-recording z-scoring, and region matching. | OK |
| Drosophila matching scale | Structure and function are matched by side-aware Ito anatomical region, not by one-to-one cell identity. | OK |
| Cross-species statistics | FCV group panels use Kruskal-Wallis plus Holm-corrected Mann-Whitney U tests; DCA panels use group-label permutation tests with Holm correction. | OK |

### Methods wording that should be revised

1. **Zebrafish high-pass frequency**

   - Manuscript currently says zebrafish calcium signals were high-pass filtered at `0.01 Hz`.
   - The final FCV/FCS/profile calculation script uses `HIGHPASS_HZ = 0.03`.
   - Recommended fix: either change the Methods to `0.03 Hz` for the FCV/FCS analysis, or explicitly state that the submitted FCV/FCS metrics were computed after an additional `0.03 Hz` high-pass filtering step.

2. **FCS definition**

   - Manuscript currently defines FCS as the mean absolute FC.
   - The high-pass FC measure code computes mean signed FC across windows and partners:
     `fc_strength_series()` uses `np.nanmean(corr, axis=1)` without an absolute value.
   - Recommended fix: replace "mean absolute FC" with "mean signed FC" or "mean FC".

3. **Regional DCA equations**

   - Manuscript equations can be read as averages over unique upstream/downstream neurons.
   - Zebrafish code averages over inter-regional edge endpoints. Thus a downstream target cell contributes once per outgoing edge instance, and C. elegans/Drosophila use synapse/edge-weighted analogues.
   - Recommended fix: define `DCA_post` and `DCA_pre` over edge sets:

     ```text
     E_out(k) = {(u,v): r(u)=k, r(v) != k}
     E_in(k)  = {(u,v): r(u) != k, r(v)=k}

     DCA_post,k = mean_{(u,v) in E_out(k)} DCA_v
     DCA_pre,k  = mean_{(u,v) in E_in(k)} DCA_u
     ```

4. **Stimulus-associated zebrafish FCV window**

   - Manuscript currently says stimulus FCV used 5 s windows and 2.5 s steps.
   - The stimulus source code uses `STIM_WINDOW_SIZE = 20` and `STIM_OVERLAP = 15`, i.e. a 20-frame window with a 5-frame step.
   - Recommended fix: write "20-frame windows with 5-frame steps (approximately 10 s and 2.5 s at the zebrafish imaging rate)."

5. **Drosophila structural-connectivity wording**

   - Manuscript describes the Drosophila SC as a cell-to-cell synapse-count matrix with edges below five synapses excluded.
   - This is accurate for the local cell-level DCA/PostDCA calculation.
   - The region-level Ito SC matrix used for network visualization is built separately by a cell pre-region/post-region presence rule:
     `SC(i,j)+=1` when a cell has presynaptic sites in Ito region `i` and postsynaptic sites in Ito region `j`.
   - Recommended fix: explicitly separate these two products:
     1. region-level Ito pre/post presence SC for anatomical network visualization and region aggregation;
     2. thresholded cell-to-cell synapse-count graph for local DCA/PostDCA.

6. **C. elegans window variants**

   - Methods correctly describe the cross-species FCV table as 20-frame/8-frame.
   - Some standalone C. elegans matrix figures use 60-frame/15-frame windows.
   - Recommended fix: no main Methods change is needed if the 20/8 statement refers to the cross-species FCV analysis, but SI/figure captions should keep the 60/15 sensitivity or standalone matrix convention distinct.

## Data and Code Location Map

The analysis workspace contains three different layers of material:

1. raw or near-raw source data,
2. original calculation code and intermediate outputs,
3. the final figure-reproduction package.

These should not be described interchangeably. The public `final_figure_pack` is the third layer: it is designed to redraw the submitted figures from processed/cached analysis products.

### Final Figure-Reproduction Package

| Purpose | Location |
|---|---|
| Final figure code, processed data, outputs, docs | `ncomms_figure_build/final_figure_pack/` |
| Final processed tables used by figure scripts | `ncomms_figure_build/final_figure_pack/data/final_summary_tables/` |
| Figure-ready C. elegans processed data | `ncomms_figure_build/final_figure_pack/data/figure14_celegans/` |
| Figure-ready Drosophila processed data | `ncomms_figure_build/final_figure_pack/data/figure15_drosophila/` |
| Figure-ready zebrafish processed/cached data | `ncomms_figure_build/final_figure_pack/data/region_community_io/`, `data/source_inputs/`, and `data/final_summary_tables/` |

### Raw or Near-Raw Source Data

| Species | Raw or near-raw data location | Notes |
|---|---|---|
| Zebrafish | `raw/zebrafish/` | Local root for zebrafish raw inputs. |
| Zebrafish | `raw/zebrafish/original_raw_data/` | Original raw-data folder. |
| Zebrafish | `raw/zebrafish/standardized_matrix/` | Standardized matrix inputs used by downstream analyses. |
| C. elegans | `fcv_postdca_raw_recompute/data/celegans/fc/atanas_kim_2023_www_archive.bz2` | WormWideWeb activity archive. |
| C. elegans | `celegans_data/herm_full_edgelist.csv` | Chemical synapse edge list. |
| C. elegans | `figure14_celegans_final/data/herm_full_edgelist.csv` | Analysis-local copy of the chemical synapse edge list. |
| Drosophila | `figure15_drosophila_final/data/raw/flywire_783/` | FlyWire-derived raw/near-raw structural files. |
| Drosophila | `figure15_drosophila_final/data/raw/flywire_783/proofread_connections_783.feather` | Proofread FlyWire connection table. |
| Drosophila | `figure15_drosophila_final/data/raw/flywire_783/per_neuron_neuropil_count_pre_783.feather` | Per-neuron presynaptic neuropil counts. |
| Drosophila | `figure15_drosophila_final/data/raw/flywire_783/per_neuron_neuropil_count_post_783.feather` | Per-neuron postsynaptic neuropil counts. |
| Drosophila | `figure15_drosophila_final/data/raw/turner_mann_clandinin/` | Turner/Mann/Clandinin Branson activity-derived source files. |
| Drosophila | `figure15_drosophila_final/data/raw/turner_mann_clandinin/data_TurnerMannClandinin.tar.gz` | Archived Branson/Turner/Mann/Clandinin calcium data source used by the original Drosophila pipeline. |

### Original Calculation Code

| Analysis step | Main code location |
|---|---|
| Zebrafish baseline FCV, FCS, and FC partner reconfiguration | `ncomms_figures_rebuild/scripts/compute_highpass_ce_zf_plot_measures.py` |
| Zebrafish DCA/PostDCA/PreDCA | `fcv_postdca_raw_recompute/code/zebrafish/compute_rank1_post_dca.py` |
| Shared directed coreness/DCA helpers | `fcv_postdca_raw_recompute/code/common/coreness.py`, `fcv_postdca_raw_recompute/code/common/post_dca.py` |
| C. elegans SC/FC matrix construction | `figure14_celegans_final/code/figure14_celegans_w60_nozero_matrices.py` |
| C. elegans spontaneous FCV/FCS measures | `figure14_celegans_final/code/figure14_celegans_fc_spontaneous_basic_measures.py` |
| C. elegans PostDCA/FCV linking | `figure14_celegans_final/code/figure14_celegans_postdca_fcv_z_nozero.py` |
| Drosophila full Branson FC/FCV | `figure15_drosophila_final/code/drosophila_branson999_full_fc_fcv.py` |
| Drosophila full FC five-measure calculation | `figure15_drosophila_final/code/drosophila_branson999_full_fc_5measures.py` |
| Drosophila FlyWire cell-level Ito DCA/PostDCA | `figure15_drosophila_final/code/drosophila_flywire783_cell_ito48_dca_postdca.py` |
| Drosophila pre/post presence SC by Ito48 | `figure15_drosophila_final/code/drosophila_flywire783_cell_prepost_presence_sc_ito48.py` |

### Processed Manuscript Data

| Purpose | Location |
|---|---|
| Main processed manuscript data | `ncomms_figure_build/data/` |
| Final summary tables used across figures | `ncomms_figure_build/data/final_summary_tables/` |
| Processed source-input copies from upstream pipelines | `ncomms_figure_build/data/source_inputs/` |
| Current manuscript source files and figures | `ncomms_figure_build/manuscript/source_reports/` |
| Current working manuscript files | `ncomms_figure_build/real_manuscript/` |

Bottom line: raw data and original calculation scripts live mostly outside `final_figure_pack`; the final pack contains the processed/cached products needed to regenerate figures and verify plotted values.

## Main Figure Code Paths

- Anatomical group panels: `code/figure7_14_combined_group_panels_AB.py`
- Functional group panels: `code/figure7_14_combined_group_panels_CDE.py`
- Combined source plotting utilities: `code/figure14_15_celegans_drosophila_AFHI_combined.py`
- FCV loading and within-recording normalization: `code/figure1_cross_species_fc_dynamics_measures.py`
- SC/DCA feature loading: `code/figure2_structural_connectivity_predictors_of_fcv.py`
- C. elegans standalone provenance/plot code: `code/figure14_celegans_combined_horizontal.py`, `code/figure14_celegans_full_combined_FINAL.py`
- Drosophila standalone provenance/plot code: `code/figure15_drosophila_full_combined_FINAL.py`

## FC, SC, DCA Calculation Logic and SC-FC Matching

This section traces the main raw/near-raw calculation code, not only the final plotting code.

### Shared DCA Definition

Core DCA helpers are in:

- `fcv_postdca_raw_recompute/code/common/post_dca.py`
- `fcv_postdca_raw_recompute/code/common/coreness.py`

The rank-1 directed DCA implementation computes non-negative output and input coreness vectors and defines:

```text
DCA = c_out - c_in
```

For C. elegans and the shared helper path, `PostDCA` is the outgoing edge-weighted mean DCA of target nodes, and `PreDCA` is the incoming edge-weighted mean DCA of source nodes. This matches the biological interpretation used in the manuscript: PostDCA describes the downstream targets reached by a node/region; PreDCA describes the upstream sources projecting into it.

### Zebrafish

**FC calculation**

- Main original code: `ncomms_figures_rebuild/scripts/compute_highpass_ce_zf_plot_measures.py`
- Raw input folder: `ncomms_figures_rebuild/raw/zebrafish_recording_pkl/`
- Key functions:
  - `zebrafish_louvain_community_root_area_traces()`
  - `measure_table_for_recording()`
  - `fcv_and_profile_corr_distance()`

Calculation trace:

1. Cell traces are loaded from each zebrafish recording pkl.
2. For each root area, within-region neuron FC is computed.
3. Louvain communities are detected within each root area.
4. Community mean traces are retained if community size and activity-variance thresholds are passed.
5. Traces are high-pass filtered at 0.03 Hz and z-scored.
6. Sliding-window Pearson FC is computed with 20-frame windows and 5-frame steps.
7. FCV is the mean, per node/community, of the temporal standard deviation of pairwise FC edges.
8. Community-level FC measures are averaged back to root-area regions.

This is consistent with the current Methods interpretation, provided the Methods state that zebrafish FCV is computed from within-root-area Louvain community mean traces and then summarized to anatomical root areas.

**SC/DCA calculation**

- Main original code: `fcv_postdca_raw_recompute/code/zebrafish/compute_rank1_post_dca.py`
- Raw input folder: `fcv_postdca_raw_recompute/data/zebrafish/sc/original_raw_data/`
- Raw file pattern: `subject_{subject}_data_cellular_synapse_sc_100_data.pkl`

Calculation trace:

1. Raw `cellular_sc_list` edges are interpreted as `(pre, post)`.
2. Each neuron is assigned to an anatomical region using `root_area[neuron_region_id]`.
3. For each region, local directed coreness is computed from the within-region directed subnetwork.
4. Cell-level DCA is `Rank1COut - Rank1CIn`.
5. Inter-region edges are selected by requiring source and target regions to differ.
6. Regional `Rank1PostDCA` is the mean DCA of target cells over outgoing inter-region edges from a source region.
7. Regional `Rank1PreDCA` is the mean DCA of source cells over incoming inter-region edges to a target region.

**SC-FC matching**

- Final matching in `figure2_structural_connectivity_predictors_of_fcv.py`, `expand_zebrafish_subject_points()`.
- FCV rows are matched to SC rows by:

```text
Subject + node
```

where `node` is the anatomical root-area name. FCV is first z-scored within each subject/recording, then merged to subject-region `PostDCA`, `PreDCA`, modularity, log out/in, and OO-fraction tables.

**Audit conclusion**

The zebrafish FCV and DCA calculations are conceptually consistent with the manuscript. The important wording detail is that regional DCA is an inter-region edge-endpoint average, not a simple average over unique downstream/upstream neurons.

### C. elegans

**SC calculation**

- Original SC matrix code: `figure14_celegans_final/code/figure14_celegans_w60_nozero_matrices.py`
- SC edge list: `figure14_celegans_final/data/herm_full_edgelist.csv`
- Helper: `figure14_celegans_validation.load_neuron_graph()`

The SC matrix is built by selecting the final neuron order and filling source-by-target chemical synapse weights from the hermaphrodite chemical synapse edge list.

**Standalone w60 no-zero SC/FC/FCV matrix pipeline**

The file `figure14_celegans_final/code/figure14_celegans_w60_nozero_matrices.py` builds the C. elegans matrix products used by the standalone C. elegans matrix/cluster analyses, not the primary cross-species 20/8 FCV summary.

Key implementation details:

1. Constants are `WINDOW = 60` and `STEP = 15`.
2. The input neuron list is loaded from `figure14_celegans_final/results/figure14_celegans_activity_fcv_neuron_summary_w60.csv`.
3. `load_selected_neurons()` keeps neurons with finite `PostDCA` and `FCV_z`, removes numerically zero `PostDCA` values using `POSTDCA_ZERO_TOL`, and applies the requested `--min-recordings` threshold.
4. The current final matrix set uses suffix `nrec3_spontaneous`, yielding 122 neurons:
   - Sensory: 53
   - Interneuron: 48
   - Motorneuron: 21
5. Neuron order is sorted by cell class and PostDCA.
6. SC is built from the Cook/OpenWorm chemical synapse graph by direct source x target synapse weights.
7. Activity is read from the Atanas/Kim WormWideWeb JSON archive. Labels are cleaned using `clean_label()` and matched to the selected neuron list.
8. The `--phase spontaneous` mode slices traces before the first annotated heat event. `--phase heat` slices after heat onset.
9. Traces are z-scored by `zscore_traces()`; this w60 matrix script does not perform the 0.03-Hz high-pass filtering used in the cross-species high-pass summary.
10. `window_corr_stack()` computes sliding-window Pearson FC matrices, masks self-connections to NaN, and requires at least three windows.
11. Pairwise FC is the mean windowed correlation for a neuron pair.
12. Pairwise FCV is the temporal SD of the same windowed correlation for that neuron pair.
13. Pairwise FC and FCV matrices are averaged over recordings in which both neurons are valid.
14. The code also writes a pair-count matrix recording how many recordings contributed to each FC/FCV edge.

Generated `nrec3_spontaneous` outputs checked:

| Output | Shape/status |
|---|---|
| `figure14_celegans_w60_nozero_matrix_neuron_order_nrec3_spontaneous.csv` | 122 neurons; `n_recordings` range 3-28; zero PostDCA count 0 |
| `figure14_celegans_w60_nozero_SC_matrix_nrec3_spontaneous.csv` | 122 x 122 |
| `figure14_celegans_w60_nozero_FC_matrix_nrec3_spontaneous.csv` | 122 x 122 |
| `figure14_celegans_w60_nozero_FCV_matrix_nrec3_spontaneous.csv` | 122 x 122; diagonal NaN; 13,634 finite off-diagonal entries |
| `figure14_celegans_w60_nozero_FC_pair_n_recordings_nrec3_spontaneous.csv` | 122 x 122 |

Methods implication: this script should be described as a standalone C. elegans matrix/cluster or sensitivity/QC pipeline using 60-frame windows and 15-frame steps. It should not be cited as the direct source of the cross-species 20-frame/8-frame neuron-level FCV values unless the text explicitly notes the different window convention.

**FC/FCV calculation**

There are two related C. elegans FCV products:

1. The standalone C. elegans matrix/heatmap pipeline uses `WINDOW=60`, `STEP=15` in `figure14_celegans_w60_nozero_matrices.py` and related `figure14_celegans_final` scripts.
2. The cross-species final FCV table uses the high-pass pipeline in `ncomms_figures_rebuild/scripts/compute_highpass_ce_zf_plot_measures.py`, with `WINDOWS["C. elegans"] = (20, 8)`.

In both cases, the FCV definition is consistent: z-scored calcium traces are converted to sliding-window Pearson FC matrices, pairwise FC standard deviation is computed across windows, and node FCV is the mean of those edge-wise standard deviations.

**DCA calculation**

- Shared formula: `fcv_postdca_raw_recompute/code/common/post_dca.py`
- Function: `celegans_node_post_dca()`

Calculation trace:

1. The SC matrix is oriented source x target.
2. Diagonal is removed for DCA calculation when `remove_diagonal=True`.
3. Rank-1 directed DCA is computed as `c_out - c_in`.
4. `PostDCA` is the weighted mean of target-node DCA over outgoing synapse-weighted edges.
5. `PreDCA` is the weighted mean of source-node DCA over incoming synapse-weighted edges.
6. OO-fraction is computed from outgoing targets, counting output-like source and target nodes where DCA > 0.

**SC-FC matching**

- Final C. elegans matching is by neuron name.
- The final table keeps neurons present in both:
  - within-recording FCV summaries, and
  - PostDCA/PreDCA DCA summaries.

Observed final matched set:

- 120 C. elegans neurons are used in the final SC-FC table.
- The two SC-subset neurons not present in the final matched values are `ASHR` and `BAGR`.

**Audit conclusion**

No major issue was found in C. elegans FC/SC/DCA matching. The main caveat is to distinguish the 60/15 standalone matrix products from the 20/8 high-pass cross-species FCV table.

### Drosophila

**FC/FCV calculation**

- Original current code: `figure15_drosophila_final/code/drosophila_branson999_full_fc_fcv.py`
- Raw data folder: `figure15_drosophila_final/data/raw/turner_mann_clandinin/data/branson_responses/`
- Final current output setting: `w15_step5_hp030`

Calculation trace:

1. Branson ROI response pickles are loaded.
2. ROI traces are linearly detrended.
3. For the final current setting, traces are high-pass filtered at 0.03 Hz.
4. ROI traces are z-scored.
5. Sliding-window Pearson FC matrices are computed.
6. FC is the mean across windows; FCV is the temporal standard deviation of FC across windows.
7. Pairwise ROI matrices are averaged over recordings in which both ROIs are present.

Important caveat:

- The script default constants are `WINDOW=60`, `STEP=15`, `HIGHPASS_HZ=0.0`, but the final current outputs use command-line settings corresponding to `w15_step5_hp030`. The output config file confirms the final settings. The default constants should therefore not be interpreted as the submitted analysis settings.

**SC calculation**

- Region-level SC presence matrix: `figure15_drosophila_final/code/drosophila_flywire783_cell_prepost_presence_sc_ito48.py`
- Cell-level DCA/PostDCA: `figure15_drosophila_final/code/drosophila_flywire783_cell_ito48_dca_postdca.py`
- Raw FlyWire inputs:
  - `proofread_connections_783.feather`
  - `per_neuron_neuropil_count_pre_783.feather`
  - `per_neuron_neuropil_count_post_783.feather`

SC matrix trace:

1. Each FlyWire proofread cell is mapped to Ito side-aware regions based on pre/post neuropil count tables.
2. The Ito48 SC matrix is built by cell-level pre-region to post-region presence: `SC(i,j)+=1` when a cell has any presynaptic site in region `i` and any postsynaptic site in region `j`.
3. This creates the side-aware directed Ito-region SC matrix used for network plotting and aggregation.

DCA trace:

1. Cells are assigned to dominant Ito side-aware regions, with confidence threshold 0.30.
2. Synaptic connections are thresholded at `syn_count >= 5`.
3. For each Ito side-aware region, cell-level local DCA is computed from internal connections.
4. Cell-level PostDCA is the outgoing synapse-count weighted mean target-cell DCA excluding same-region targets.
5. Region summaries are arithmetic summaries of cell-level values within each Ito region.

**SC-FC matching**

- Matching code: `ncomms_figure_build/final_figure_pack/code/figure15_drosophila_full_combined_FINAL.py`
- Key functions:
  - `matched_node_table()`
  - `aggregate_sc_to_matched_nodes()`
  - `aggregate_999_matrix_to_matched_nodes()`
  - `recording_node_timeseries()`

Matching trace:

1. Branson ROI labels are mapped to `atlas_region` labels.
2. ROI-level FC/FCV matrices are collapsed to side-aware atlas regions by averaging all ROI-block values.
3. FlyWire Ito SC labels are collapsed to matched side-aware atlas labels; MB compartments are collapsed to `MB_L`/`MB_R`, and AOTU labels are collapsed consistently.
4. The aggregated SC matrix is reindexed to the matched FCV atlas-region labels and diagonal is set to zero.
5. Final DCA-linked analyses use the inner intersection of FCV, PostDCA, and PreDCA availability.

Observed final matched sets:

- Drosophila FCV anatomy panel: 67 side-aware FCV nodes.
- Drosophila PostDCA summary: 41 nodes.
- Drosophila PreDCA summary: 37 nodes.
- Final DCA-linked SC-FC table: 37 nodes.
- Missing from final 37-node DCA-matched table among the 41 SC nodes: `NO`, `OTU_L`, `OTU_R`, and `PB`.

**Audit conclusion**

Drosophila FCV, SC, DCA, and matching logic are internally consistent. The main manuscript/README caveat is that FCV-only anatomy panels and DCA-linked SC-FC panels use different node sets.

### Overall Calculation Audit Conclusion

No major calculation error was found in the FC/SC/DCA code paths inspected. The main points to document clearly are:

1. `final_figure_pack` is a processed figure-reproduction package, not a full raw-to-figure recomputation archive.
2. Zebrafish regional DCA is computed over inter-region edge endpoints.
3. C. elegans has both 60/15 standalone matrix products and 20/8 cross-species high-pass FCV products.
4. Drosophila final FCV uses `w15_step5_hp030`, even though the original script defaults differ.
5. Drosophila FCV anatomy panels and DCA-linked SC-FC analyses use different matched node sets.

## FCV Calculation and Normalization

The cross-species FCV panels use `data/final_summary_tables/highpass_ce_zf_plot_measures_recording_node.csv`, then z-score plotted measures within each recording before averaging by node or plotting recording-level distributions. This is implemented in `compute_recording_zscore_node_summary()` and `load_species_recording_table()` in `figure1_cross_species_fc_dynamics_measures.py`.

In this table, the manuscript term `FCV` corresponds to the code/table column `EdgeStdFCV`: the mean edge-wise temporal standard deviation of sliding-window FC.

Observed recording/node counts in the bundled final table:

| Species | Recording-level rows | Recordings | FCV nodes | Window config |
|---|---:|---:|---:|---|
| C. elegans | 1,270 | 18 | 120 | `w20_step8_highpass` |
| Drosophila | 1,326 | 20 | 67 | `w15_step5_highpass_atlas_mean_trace` |
| Zebrafish | 450 | 7 | 66 | `baseline_highpass` |

C. elegans FCV is therefore used at neuron scale. Drosophila FCV is used at side-aware atlas-region scale for recording-level FCV anatomy panels, but only the subset with matched SC/DCA is used in DCA-linked panels.

## C. elegans SC and DCA

Primary bundled inputs:

- `data/figure14_celegans/data/herm_full_edgelist.csv`
- `data/figure14_celegans/matrices/figure14_celegans_w60_nozero_SC_matrix_nrec3_spontaneous.csv`
- `data/figure14_celegans/results/figure14_celegans_sc_cell_measures_full297_subset122.csv`
- `data/figure14_celegans/results/figure14_celegans_fc_spontaneous_5measure_summary.csv`
- `data/figure14_celegans/results/figure14_celegans_fc_spontaneous_basic_recording_level.csv`

QC table:

- `data/figure14_celegans/results/figure14_celegans_sc_cell_measures_full297_subset122_qc.csv`

Key QC values:

- Full SC cells used for measure calculation: 297
- Plotted/matched SC subset: 122
- Full SC nonzero edges: 3,604
- Full SC weight sum: 20,484
- Diagonal removed: yes
- Missing PostDCA/PreDCA/OO fraction in the 122-neuron subset: 0

Code check:

- `DCA` equals `c_out - c_in` in the bundled C. elegans SC table, with numerical max absolute difference `1.18e-16`.
- FC/FCV matching in `figure14_celegans_combined_horizontal.py` keeps neurons with finite FC and FCV matrices, then merges SC `PostDCA` and FCV summaries by neuron name.
- The final SC-FC table contains 120 matched C. elegans neurons. The two SC-subset neurons not present in final matched FCV/DCA values are `ASHR` and `BAGR`.

Conclusion: C. elegans FCV, SC, and DCA matching are internally consistent in the bundled final tables.

## Drosophila SC, FC, and Matching

Primary bundled inputs:

- FlyWire/Ito SC: `data/figure15_drosophila/final_results/SC_FC_FCV_current_standard/SC_flywire783_ito_R_then_L_matrix_FINAL.csv`
- Branson ROI order: `data/figure15_drosophila/final_results/Branson999_full_FC_FCV/Branson999_full_roi_order_FINAL.csv`
- Branson FC/FCV matrices: `Branson999_full_FC_matrix_w15_step5_hp030_FINAL.csv`, `Branson999_full_FCV_matrix_w15_step5_hp030_FINAL.csv`
- Recording-level FCV/FCS table: `Branson999_full_ROI_5measure_recording_level_w15_step5_hp030_FINAL.csv`
- Matched SC/DCA table: `data/figure15_drosophila/results/drosophila_flywire783_matched41_sc_cell_measures.csv`

QC table:

- `data/figure15_drosophila/results/drosophila_flywire783_sc_cell_measure_qc.csv`

Key QC values:

- Confident cells: 119,866
- Inter-region edges after threshold: 735,710
- Synapse threshold: 5
- Ito48 regions: 48
- Matched side-aware regions: 41
- Analysis definition: cell-level within-Ito-region DCA; inter-region weighted Pre/Post-DCA; arithmetic region summaries.

FCV config:

- `data/figure15_drosophila/final_results/Branson999_full_FC_FCV/Branson999_full_FC_FCV_config_w15_step5_hp030_FINAL.csv`
- Recording count: 20
- Atlas ROIs: 999
- Sampling rate: 1.2 Hz
- Window: 15 frames, 12.5 s
- Step: 5 frames, 4.17 s
- Slow fluctuation removal: linear detrend plus Butterworth high-pass
- High-pass: 0.03 Hz

Matching code:

- `matched_node_table()` merges Branson ROI labels with atlas regions and side labels.
- `aggregate_sc_to_matched_nodes()` collapses FlyWire Ito regions to matched side-aware nodes, including MB and AOTU collapsing, then sums row/column weights and removes the diagonal.
- `aggregate_999_matrix_to_matched_nodes()` maps ROI-level FC/FCV matrices to matched atlas regions by averaging ROI blocks.
- `recording_node_timeseries()` averages z-scored ROI traces within each atlas region.

Observed matching:

- Drosophila recording-level FCV table has 67 side-aware FCV nodes.
- Post-DCA summary has 41 nodes.
- Pre-DCA summary has 37 nodes.
- Final SC-FC matched table uses the inner FCV/PostDCA/PreDCA intersection: 37 Drosophila nodes.
- The 41-node FlyWire matched SC table has four regions not present in the final 37-node DCA-matched values: `NO`, `OTU_L`, `OTU_R`, and `PB`.

Conclusion: Drosophila FCV and DCA matching is internally consistent, but FCV-only anatomy panels and DCA-matched SC-FC panels use different node sets. This should be kept explicit in Methods/captions if needed.

## Focused Methods Check: C. elegans and Drosophila

This subsection compares the current `Research_report.tex` Methods wording directly against the calculation code for the two comparative species.

### C. elegans Methods Check

Current Methods statement:

- SC is from adult hermaphrodite chemical synapse data distributed through the OpenWorm CElegansNeuroML mirror of WormWiring/Cook 2019.
- 302 annotated neurons are reduced to 297 neurons with chemical synaptic connections.
- Activity is from WormWideWeb Atanas/Kim 2023.
- Labels are cleaned and matched to connectome neuron names.
- Spontaneous epochs are before heat onset; heat epochs are after heat onset.
- Spontaneous calcium traces are high-pass filtered at 0.03 Hz and z-scored.
- Cross-species FCV uses 20-frame windows and 8-frame steps.
- DCA is computed at neuron level; PostDCA is outgoing synapse-weighted target DCA; PreDCA is incoming synapse-weighted source DCA.

Code agreement:

- `figure14_celegans_validation.py` builds a directed chemical SC matrix from `herm_full_edgelist.csv`, keeps chemical synapses among sensory, interneuron, and motor-neuron classes, and yields 297 neurons after excluding disconnected/non-chemical cases.
- `ncomms_figures_rebuild/scripts/compute_highpass_ce_zf_plot_measures.py::compute_celegans_highpass()` loads `*_raw_traces.pkl`, applies `highpass_filter(..., HIGHPASS_HZ=0.03)`, z-scores rows, and computes measures with `WINDOWS["C. elegans"] = (20, 8)`.
- The final bundled table contains 18 recordings, 1,270 recording-neuron rows, and 120 final FCV nodes for C. elegans.
- `fcv_postdca_raw_recompute/code/common/post_dca.py::celegans_node_post_dca()` computes DCA as `c_out - c_in`, PostDCA as outgoing synapse-weighted target DCA, and PreDCA as incoming synapse-weighted source DCA.

Recommended Methods wording refinements:

1. The sentence "Across the 18 spontaneous recordings used for FCV analysis, 148 cleaned neuron labels were available" is not wrong as an intermediate label-coverage statement, but the final matched FCV/DCA table contains 120 neurons. The Methods should add: "The final matched C. elegans SC-FC analysis retained 120 neurons with both FCV and DCA summaries."
2. Keep the 60-frame/15-frame C. elegans matrix pipeline out of the main cross-species Methods unless it is explicitly described as a standalone matrix/sensitivity analysis. The main cross-species Methods should remain 20-frame/8-frame.

### Drosophila Methods Check

Current Methods statement:

- SC is from FlyWire783.
- Cell-to-cell synapse-count edges below five synapses are excluded.
- Each FlyWire cell is assigned to a dominant side-aware Ito L/R region with at least 30% localization confidence.
- Functional data are from Branson999 Turner/Mann/Clandinin calcium recordings.
- ROI traces are matched to side-aware Ito regions and averaged within region.
- Regional traces are detrended, high-pass filtered at 0.03 Hz, z-scored, and analyzed with 15-frame windows and 5-frame steps.
- FCV-centered summary retains 41 side-aware Ito nodes; final SC-FC comparison uses 37 nodes.

Code agreement:

- `drosophila_branson999_full_fc_fcv.py` confirms the final `w15_step5_hp030` setting: 20 recordings, 999 atlas ROIs, sampling rate 1.2 Hz, 15-frame windows, 5-frame steps, linear detrend plus 0.03-Hz Butterworth high-pass.
- `compute_drosophila_atlas_mean_trace_recording_nodes()` in `compute_highpass_ce_zf_plot_measures.py` averages ROI traces by `atlas_region`, detrends, high-pass filters, z-scores, and computes the same sliding-window measures.
- `drosophila_flywire783_cell_ito48_dca_postdca.py` confirms dominant Ito mapping with confidence threshold 0.30, synapse threshold `syn_count >= 5`, cell-level within-Ito-region DCA, and inter-region synapse-weighted PostDCA.
- The final QC table reports 119,866 confident cells, 735,710 inter-region edges after thresholding, 48 Ito regions, 41 matched regions, and final DCA-linked analyses on 37 nodes.

Recommended Methods wording refinements:

1. Separate the two Drosophila SC products more explicitly:
   - Region-level Ito SC/network visualization uses the cell pre-region/post-region presence matrix (`SC(i,j)+=1` if a cell has presynaptic sites in region `i` and postsynaptic sites in region `j`).
   - Local DCA/PostDCA uses the thresholded cell-to-cell synapse-count graph with `syn_count >= 5`.
2. The Methods already correctly state the 41-node FCV-centered set and 37-node final SC-FC subset. Keep this distinction because FCV anatomy panels and DCA-linked panels do not use identical node sets.
3. If the text says "summed synapse count" for all Drosophila SC analyses, clarify that this applies to the DCA cell-to-cell graph, not the region-level pre/post presence SC matrix.

Suggested replacement for the first Drosophila SC paragraph:

```text
The Drosophila melanogaster structural analysis used FlyWire783 cell-level annotations. Two related SC products were constructed. For anatomical network visualization and region-level aggregation, we built a side-aware Ito-region SC matrix using a cell pre-region/post-region presence rule: for each proofread cell, every Ito region containing at least one presynaptic site was paired with every Ito region containing at least one postsynaptic site, and the corresponding region-level edge count was incremented by one. For local DCA and PostDCA calculations, we used the thresholded cell-to-cell FlyWire synapse graph, in which each edge connected one presynaptic cell to one postsynaptic cell and edge weight was summed synapse count; edges with fewer than five synapses were excluded. Each FlyWire cell was assigned to a dominant side-aware Ito L/R neuropil region from synaptic localization counts, and cells were retained when the dominant region accounted for at least 30% of the cell's total localization count.
```

## Cross-Species SC-FC Feature Table

Final matched table:

- `data/final_summary_tables/sc_four_measures_vs_fcv_all_species_values.csv`

Matched node counts:

| Species | Matched SC-FC nodes |
|---|---:|
| C. elegans | 120 |
| Drosophila | 37 |
| Zebrafish | 66 subject-region rows / 47 unique region-level SC nodes depending on panel |

For C. elegans and Drosophila, `PreDCA` in this table matches the bundled `oo_fraction_recomputed_values_by_species.csv` source exactly. `PostDCA` is taken from `edge_target_dca_distribution_summary_by_unit.csv`.

## Code Fix Applied During Audit

`code/figure2_structural_connectivity_predictors_of_fcv.py` still pointed to older unbundled filenames:

- `prepost_dca_signed_plus_minus_edge_fcv_node_values.csv`
- `oo_ii_fraction_recomputed_values_by_species.csv`
- separate modularity/log-ratio summary files that are not included in the current final pack

The loader was updated to use the bundled `oo_fraction_recomputed_values_by_species.csv` table for `PreDCA`, modularity, log(out/in), and OO fraction. After the fix, `load_values()` successfully rebuilds matched SC-FC values with:

- C. elegans: 120 nodes, no FCV/PostDCA/PreDCA missing values
- Drosophila: 37 nodes, no FCV/PostDCA/PreDCA missing values
- Zebrafish: 47 region-level SC nodes

Drosophila modularity remains missing for two matched nodes, consistent with the existing bundled final values.

## Remaining Notes

- `figure15_drosophila_full_combined_FINAL.py` still contains older `w30/step8` cached paths and a fallback function using `window=30, step=8`. The current final cross-species FCV table and QC files use `w15/step5/hp030`. The old paths appear to be legacy/fallback material, but they should be cleaned or clearly marked if the standalone Drosophila script is presented as a fully current raw-recompute script.
- The final figure pack reproduces figures from processed/cached analysis outputs. It does not fully regenerate C. elegans and Drosophila FCV/DCA from all raw activity/connectome files from scratch.
- For manuscript clarity, state that Drosophila FCV anatomy panels use all available side-aware FCV nodes, whereas DCA-linked Drosophila analyses use the stricter matched 37-node FCV/PostDCA/PreDCA intersection.

## Raw-Data Boundary and Matrix QC

The bundled `data/README.md` correctly states that the package contains processed, figure-ready inputs rather than a complete raw-data archive. This was confirmed directly:

| Raw input type | Bundled in final pack? | Note |
|---|---:|---|
| C. elegans WormWideWeb activity archive | No | `atanas_kim_2023_www_archive.bz2` is not included. |
| Drosophila Branson raw pickle recordings | Partial | Only two `branson_*.pkl` files are bundled, while final FCV uses 20 recordings. |
| Zebrafish `subject_*_data_cellular_synapse_sc_100_data.pkl` raw/cache files | No | Stimulus and SC analyses use processed/cached tables in the final pack. |

Therefore, the final pack should be described as a figure-reproduction package based on processed analysis outputs, not as a complete raw-to-figure recomputation archive.

Processed matrix integrity checks were nevertheless consistent with the manuscript-level calculations:

### C. elegans

- SC edgelist source: `data/figure14_celegans/data/herm_full_edgelist.csv`.
- FC and FCV matrices are symmetric to numerical precision, with self-pairs excluded as non-finite values.
- Recording-level FCV table contains 1,270 rows, 18 recordings, and 120 matched neurons, with no missing `FCV_z` values.
- SC/DCA QC reports diagonal removal for the DCA calculation: original diagonal weight sum 105, final diagonal weight sum 0.
- The displayed/cache SC matrix still contains diagonal entries. This is acceptable if it is treated as a cached SC matrix, but DCA methods should refer to the diagonal-removed DCA calculation rather than implying every stored SC matrix has zero diagonal.

### Drosophila

- FC and FCV matrices for the current final setting (`w15_step5_hp030`) are symmetric to numerical precision.
- Recording-level FCV table contains 13,842 rows from 20 recordings, 719 ROI labels, and 67 atlas regions, with no missing `FCV_z` values.
- FlyWire SC matrix is directed and therefore asymmetric, as expected.
- SC/DCA QC reports 119,866 confident cells, 735,710 inter-region edges after thresholding, a synapse threshold of 5, 48 Ito regions, and 41 matched side-aware regions.
- The final SC-FC matched analysis uses 37 Drosophila regions after requiring FCV, PostDCA, and PreDCA to all be present.

### Zebrafish

- Final-pack zebrafish FC/SC figure scripts primarily read processed tables and cached subject-level TE/SC/FC outputs.
- The original raw SC construction from synapse-like endpoints and neuron assignments is upstream of the final figure pack.
- Existing audit notes identify the relevant original-code trace as `fcv_postdca_raw_recompute/code/zebrafish/compute_rank1_post_dca.py` for DCA and `ncomms_figures_rebuild/scripts/compute_highpass_ce_zf_plot_measures.py` for baseline FCV/FCS/profile reconfiguration.
- At the processed-table level, the final zebrafish values match the reported figure statistics, but the public final pack alone cannot prove a fresh raw endpoint-to-SC reconstruction.

Bottom line: no major internal inconsistency was found in the processed SC/FC/FCV/DCA tables used by the final figures. The main caveat is scope: the final pack validates and reproduces from processed analysis products, while full raw-data regeneration requires upstream raw datasets and scripts outside the bundled figure pack.
