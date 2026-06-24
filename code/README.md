# Code Guide

This folder contains the Python scripts used to regenerate the manuscript figures from the bundled processed inputs.

Run scripts from the repository root, not from inside `code/`:

```bash
cd final_figure_pack
python3.10 code/figure9_clean.py
```

## Main Figure Scripts

| Script | Output |
|---|---|
| `figure9_clean.py` | `figures/figure9_final.png` |
| `figure12_clean.py` | `figures/figure12_final.png` |
| `figure_sc_fc_final_overview.py` | `figures/figure_sc_fc_final_overview.png` |
| `figure13_clean.py` | `figures/figure13_final.png` |
| `figure7_14_combined_group_panels_AB.py` | `figures/figure7_14_combined_group_panels_AB.png` |
| `figure7_14_combined_group_panels_CDE.py` | `figures/figure7_14_combined_group_panels_CDE.png` |

## Supplementary Figure Scripts

Supplementary scripts use two historical working-output conventions. Some write PNGs to `figures/`; others create a temporary `output/png/` folder. Running `scripts/sync_outputs.sh` stages every supplementary PNG into `figures/` and then copies the clean final PNG set to `outputs/supplementary/`. PDF files are not retained in this pack.

| Script | Direct script output | Staged working output | Final synchronized output |
|---|---|---|---|
| `figure_supply_0.py` | `output/png/` | `figures/figure_supply_0.png` | `outputs/supplementary/figure_supply_0.png` |
| `figure_supply_1.py` | `output/png/` | `figures/figure_supply_1.png` | `outputs/supplementary/figure_supply_1.png` |
| `figure_supply_2_proc.py` | `figures/` | `figures/figure_supply_2_proc.png` | `outputs/supplementary/figure_supply_2_proc.png` |
| `figure_supply_5.py` | `output/png/` | `figures/figure_supply_5.png` | `outputs/supplementary/figure_supply_5.png` |
| `figure_supply_10_proc.py` | `figures/` | `figures/figure_supply_10_proc.png` | `outputs/supplementary/figure_supply_10_proc.png` |
| `figure_supply_13.py` | `output/png/` | `figures/figure_supply_13.png` | `outputs/supplementary/figure_supply_13.png` |
| `figure_supply_14.py` | `figures/` and `output/stats/` | `figures/figure_supply_14.png` | `outputs/supplementary/figure_supply_14.png` |
| `figure_supply_15.py` | `output/png/` | `figures/figure_supply_15.png` | `outputs/supplementary/figure_supply_15.png` |

## Shared Helpers

- `figure_style.py`: shared fonts, colors, and panel-label style.
- `figure_supply_sc_heatmap.py`: helper for SC heatmap panels.
- `figure_regionwise_multivariate_coupling.py`: helper for SC-FC overview panels.
- Cross-species helper modules such as `figure14_15_celegans_drosophila_AFHI_combined.py` and `figure7_by_species_with_figure5_abc.py` are retained because the split manuscript Fig. 5/6 scripts import their plotting and statistics utilities.

The scripts assume the bundled folder layout in this package. Avoid moving files between folders unless the corresponding paths in the scripts are updated.
