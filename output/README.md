# Working Supplementary Outputs

This folder is a working output location used by selected supplementary scripts.

- PNG files are written to `output/png/`.
- Supplementary statistics from selected scripts are written to `output/stats/`.

This folder is not the final figure collection. Running `scripts/sync_outputs.sh` stages PNG files into `figures/` and then copies the clean final PNG set to `outputs/supplementary/`. PDF files are not retained in this pack.

After rerunning scripts, update `outputs/` with:

```bash
bash scripts/sync_outputs.sh
```
