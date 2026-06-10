#!/usr/bin/env bash
set -euo pipefail

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mpl_config}"
export PYTHONDONTWRITEBYTECODE=1

python3.10 code/figure9_clean.py
python3.10 code/figure12_clean.py
python3.10 code/figure_sc_fc_final_overview.py
python3.10 code/figure13_clean.py
python3.10 code/figure14_15_celegans_drosophila_AFHI_combined.py
python3.10 code/figure7_by_species_with_figure5_abc.py

python3.10 code/figure_supply_0.py
python3.10 code/figure_supply_2_proc.py
python3.10 code/figure_supply_5.py
python3.10 code/figure_supply_10_proc.py
python3.10 code/figure_supply_13.py
python3.10 code/figure_supply_14.py
python3.10 code/figure_supply_15.py

bash scripts/sync_outputs.sh

