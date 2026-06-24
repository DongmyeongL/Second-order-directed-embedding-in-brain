# Manifest

This manifest summarizes the upload-ready contents of `final_figure_pack`.

## Package Purpose

Standalone figure reproduction package for the manuscript main and supplementary figures. It includes processed data, plotting scripts, generated statistics, and final PNG figure outputs.

## Primary Entry Points

- `README.md`: start here.
- `outputs/`: final figures.
- `code/`: plotting scripts.
- `scripts/run_all_figures.sh`: rerun all included figure scripts.
- `scripts/sync_outputs.sh`: copy regenerated outputs into `outputs/`.

## Included Main Figure Outputs

- `outputs/figure9_final.png`
- `outputs/figure12_final.png`
- `outputs/figure_sc_fc_final_overview.png`
- `outputs/figure13_final.png`
- `outputs/figure7_14_combined_group_panels_AB.png`
- `outputs/figure7_14_combined_group_panels_CDE.png`

## Included Supplementary Figure Outputs

- `outputs/supplementary/figure_supply_0.png`
- `outputs/supplementary/figure_supply_1.png`
- `outputs/supplementary/figure_supply_2_proc.png`
- `outputs/supplementary/figure_supply_5.png`
- `outputs/supplementary/figure_supply_10_proc.png`
- `outputs/supplementary/figure_supply_13.png`
- `outputs/supplementary/figure_supply_14.png`
- `outputs/supplementary/figure_supply_15.png`

## Excluded

- Full raw imaging datasets.
- Full raw anatomical reconstruction datasets.
- Exploratory 3D functional-group figure files.
- PDF figure outputs.
- Python bytecode/cache files.

## Data Scope

Only figure-ready processed inputs and statistics needed by the retained main and supplementary figure scripts are included. Broader audit notes, exploratory outputs, full raw datasets, and unused intermediate tables were removed to keep the package compact.
