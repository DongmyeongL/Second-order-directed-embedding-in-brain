# Working Figure Outputs

This folder is the complete working figure set after synchronization.

Main figure scripts write PNG files here directly. Some supplementary scripts also write PNG files here directly; the remaining supplementary PNG outputs are first generated in `output/png/`, then staged here by `scripts/sync_outputs.sh`.

For the clean final upload/viewing PNG set, use:

- `outputs/` for main figures
- `outputs/supplementary/` for supplementary figures

After rerunning scripts, update `outputs/` with:

```bash
bash scripts/sync_outputs.sh
```
