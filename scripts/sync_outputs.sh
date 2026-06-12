#!/usr/bin/env bash
set -euo pipefail

mkdir -p figures outputs/supplementary

# Stage all supplementary figures into figures/ first. Some supplementary
# scripts write directly to figures/, while others write to output/png; this
# step makes figures/ the complete working PNG figure set.
cp output/png/figure_supply_0.png figures/
cp output/png/figure_supply_1.png figures/
cp output/png/figure_supply_5.png figures/
cp output/png/figure_supply_13.png figures/
cp output/png/figure_supply_15.png figures/

cp figures/figure9_final.png outputs/
cp figures/figure12_final.png outputs/
cp figures/figure_sc_fc_final_overview.png outputs/
cp figures/figure13_final.png outputs/
cp figures/figure14_15_celegans_drosophila_AFHI_combined.png outputs/
cp figures/figure7_by_species_with_figure5_abc.png outputs/

cp figures/figure_supply_2_proc.png outputs/supplementary/
cp figures/figure_supply_10_proc.png outputs/supplementary/
cp figures/figure_supply_14.png outputs/supplementary/

cp figures/figure_supply_0.png outputs/supplementary/
cp figures/figure_supply_1.png outputs/supplementary/
cp figures/figure_supply_5.png outputs/supplementary/
cp figures/figure_supply_13.png outputs/supplementary/
cp figures/figure_supply_15.png outputs/supplementary/

find figures output outputs -type f -name '*.pdf' -delete

echo "Staged all PNG figures in figures/ and synchronized final PNG figures into outputs/."
