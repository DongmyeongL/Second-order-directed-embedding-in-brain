#!/usr/bin/env bash
set -euo pipefail

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mpl_config}"
export PYTHONDONTWRITEBYTECODE=1

mkdir -p figures outputs outputs/supplementary output/png output/stats

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

bash scripts/sync_outputs.sh
