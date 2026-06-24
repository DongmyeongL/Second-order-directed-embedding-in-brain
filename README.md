# Final Figure Pack

This repository folder contains the code, compact figure-generation data, statistics tables, and current outputs for the manuscript figures. It is intended to be uploaded as a standalone package: a reader should be able to inspect the final figures immediately and rerun the plotting scripts without access to the original analysis workspace.

The pack includes figure-ready processed inputs only. It does not include the full raw calcium-imaging, synapse-reconstruction, FlyWire, or external connectome datasets.

## What To Look At First

- Final main figures: `outputs/`
- Final supplementary figures: `outputs/supplementary/`
- Figure-generation scripts: `code/`
- One-command rerun scripts: `scripts/`
- Processed data used by the scripts: `data/`
- Key statistics and figure tables: `data/final_summary_tables/`
- Package inventory: `MANIFEST.md`

Most users only need `outputs/` to inspect the finished figures. Users who want to reproduce or modify figures should start from the tables below and run the corresponding script from this folder.

## Folder Structure

| Folder | Contents |
|---|---|
| `code/` | Python scripts for all included main and supplementary figures, plus shared style/helper modules. |
| `data/` | Bundled processed inputs required to redraw figures. |
| `data/final_summary_tables/` | Main statistical summaries, correlations, model outputs, and panel-level values. |
| `figures/` | Working figure output folder used by plotting scripts. |
| `outputs/` | Clean final PNG collection copied from `figures/` for easy viewing, manuscript editing, or upload. Use this folder to inspect all final figures. |
| `scripts/` | Convenience shell scripts for rerunning all figures and syncing final outputs. |
| `MANIFEST.md` | Short inventory of included outputs and exclusions. |
| `requirements.txt` | Python packages used by the figure scripts. |

Additional folder-specific guides are available at `code/README.md`, `data/README.md`, and `outputs/README.md`.

## Main Figure Map

| Manuscript figure | Main topic | Final output | Script |
|---|---|---|---|
| Fig. 1 | Zebrafish functional dynamics and directed TE overview | `outputs/figure9_final.png` | `code/figure9_clean.py` |
| Fig. 2 | Zebrafish directed structural-connectivity organization | `outputs/figure12_final.png` | `code/figure12_clean.py` |
| Fig. 3 | SC-FC coupling, PLS, and FCV prediction | `outputs/figure_sc_fc_final_overview.png` | `code/figure_sc_fc_final_overview.py` |
| Fig. 4 | Directed asymmetry simulations and perturbation models | `outputs/figure13_final.png` | `code/figure13_clean.py` |
| Fig. 5 | C. elegans and Drosophila anatomical-class analyses | `outputs/figure7_14_combined_group_panels_AB.png` | `code/figure7_14_combined_group_panels_AB.py` |
| Fig. 6 | Functional-group synthesis across species | `outputs/figure7_14_combined_group_panels_CDE.png` | `code/figure7_14_combined_group_panels_CDE.py` |

Only PNG versions are retained in this pack.

## Supplementary Figure Map

| Supplementary output | Main topic | Script |
|---|---|---|
| `outputs/supplementary/figure_supply_0.png` | FCV/window or supporting overview analysis | `code/figure_supply_0.py` |
| `outputs/supplementary/figure_supply_1.png` | Synapse endpoint-distance and distance--FC model inputs | `code/figure_supply_1.py` |
| `outputs/supplementary/figure_supply_2_proc.png` | Zebrafish regional SC topology and directional bias | `code/figure_supply_2_proc.py` |
| `outputs/supplementary/figure_supply_5.png` | Supporting simulation or null-model analysis | `code/figure_supply_5.py` |
| `outputs/supplementary/figure_supply_10_proc.png` | Zebrafish regional FC dynamics and directed information flow | `code/figure_supply_10_proc.py` |
| `outputs/supplementary/figure_supply_13.png` | Supporting trace/example analysis | `code/figure_supply_13.py` |
| `outputs/supplementary/figure_supply_14.png` | Zebrafish OMR stimulus-associated FCV modulation | `code/figure_supply_14.py` |
| `outputs/supplementary/figure_supply_15.png` | Supporting layer/trace example analysis | `code/figure_supply_15.py` |

## Quick Start

Run scripts from the root of this folder:

```bash
cd final_figure_pack
export MPLCONFIGDIR=/tmp/mpl_config
python3.10 code/figure9_clean.py
```

The `MPLCONFIGDIR` line is optional, but useful on systems where the default Matplotlib cache directory is not writable.

To install dependencies in a fresh Python environment:

```bash
python3.10 -m pip install -r requirements.txt
```

## Reproduce All Included Figures

```bash
cd final_figure_pack
export MPLCONFIGDIR=/tmp/mpl_config

python3.10 code/figure9_clean.py
python3.10 code/figure12_clean.py
python3.10 code/figure_sc_fc_final_overview.py
python3.10 code/figure13_clean.py
python3.10 code/figure7_14_combined_group_panels_AB.py
python3.10 code/figure7_14_combined_group_panels_CDE.py

python3.10 code/figure_supply_0.py
python3.10 code/figure_supply_1.py
python3.10 code/figure_supply_2_proc.py
python3.10 code/figure_supply_5.py
python3.10 code/figure_supply_10_proc.py
python3.10 code/figure_supply_13.py
python3.10 code/figure_supply_14.py
python3.10 code/figure_supply_15.py
```

Regenerated files are written by the individual scripts to their working output locations: mostly `figures/` for main figures, and either `figures/` or a temporary `output/` folder for supplementary figures. Running `scripts/sync_outputs.sh` stages all supplementary PNG outputs into `figures/`, then copies the clean final PNG set into `outputs/`. Any generated PDF files are removed by the sync script.

The same sequence can be run with:

```bash
bash scripts/run_all_figures.sh
```

After manually rerunning selected scripts, synchronize the final upload folder with:

```bash
bash scripts/sync_outputs.sh
```

## Software Environment

The scripts were tested with Python 3.10. They use standard scientific Python packages, including:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `seaborn`
- `networkx`
- `scikit-learn`
- `statsmodels`
- `PyYAML`

Some scripts also use image or plotting utilities such as `Pillow` through Matplotlib workflows.

## Data Scope

This pack is designed for figure reproduction, not full raw-data reanalysis. Included data are processed summaries, figure-specific matrices, cached model outputs, and compact network/trace inputs required by the scripts. Full source datasets should be cited and obtained from the original publications or public repositories described in the manuscript.

## Statistics Tables

Panel-level statistics generated by the final figure scripts are stored in `data/final_summary_tables/`. Examples include `figure9_stats.csv`, `figure12_stats.csv`, `figure13_stats.csv`, `figure_sc_fc_final_overview_stats.csv`, and the functional-group synthesis tables used by manuscript Fig. 6.

Only the files needed to redraw the manuscript main and supplementary figures are retained. Broader audit notes and exploratory outputs were removed from this compact figure pack.

## Notes For Editing

- Shared colors, fonts, and panel-label styling are centralized in `code/figure_style.py`.
- The exploratory 3D functional-group figure is intentionally excluded.
- Current outputs have been regenerated and synchronized into `outputs/`.
- If a script is modified, rerun it from this folder so relative paths resolve correctly.
