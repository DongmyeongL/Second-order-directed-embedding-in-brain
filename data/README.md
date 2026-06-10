# Data Guide

This folder contains processed, figure-ready inputs required by the plotting scripts. It is not a raw-data archive.

## Main Subfolders

| Folder or file | Purpose |
|---|---|
| `final_summary_tables/` | Core statistics tables, panel values, model outputs, correlations, and group summaries used across figures. |
| `source_inputs/` | Processed source tables copied from the analysis workspace when a figure needs them directly. |
| `figure13/` | Cached outputs for the layer-wise directed linear model and zebrafish simulation panels. |
| `figure14_celegans/` | Processed C. elegans network, activity, and representation-cache inputs. |
| `figure15_drosophila/` | Processed Drosophila FlyWire/Branson matching and FC/FCV inputs. |
| `region_community_io/` | Zebrafish subject-level community and transfer-entropy inputs for supplementary FC/TE panels. |
| `network_diagrams/` | Network diagram images used in simulation figures. |
| `figure5_simulation/`, `figure6_NULL_P_SP/` | Cached zebrafish simulation/null-model data. |
| `figure1_emp_variation/` | Supporting zebrafish empirical FC/SC variation inputs. |
| top-level `.csv`, `.npz`, `.pkl`, `.png` files | Small legacy-compatible processed inputs still used directly by selected figure scripts. |

## Data Scope

The included files are sufficient for reproducing the bundled figures. Full raw imaging, anatomical reconstruction, and external connectome datasets should be obtained from the original sources cited in the manuscript.

Large files in this folder are still below the standard GitHub single-file limit of 100 MB.

