# Drosophila Branson999 Epoch Check

Checked on 2026-05-27.

## Local file inspection

Raw files:

```text
fcv_postdca_raw_recompute/data/drosophila/fc/branson_responses/branson_*.pkl
```

There are 20 local `branson_*.pkl` recordings. Each pickle loads as a pandas
DataFrame:

- rows: Branson999 ROI indices;
- columns: time frames;
- values: numeric calcium trace values;
- no columns or object fields for `stimulus`, `response`, `trial`, `rest`,
  `spontaneous`, or `epoch`.

Observed shapes:

- most files: approximately 605--705 ROIs by 2000 frames;
- `branson_2017-11-08_1.pkl` and `branson_2017-11-08_2.pkl`: 704 ROIs by 4000 frames;
- no missing values were detected in the inspected files.

## Interpretation

The local folder name `branson_responses` is potentially misleading. The pkl
objects themselves do not contain explicit stimulus or response epoch labels.
The source Mann/Gallen/Clandinin Drosophila intrinsic FC paper describes the
functional recordings as calcium imaging acquired in the absence of sensory
stimuli. Therefore, for the current manuscript these files should be described
as Branson999 resting/intrinsic or spontaneous calcium recordings, not as
stimulus-response recordings.

Consequently, Figures 1--4 can use the full Branson999 traces as the
*Drosophila* baseline FC data. There is no separate resting/spontaneous segment
inside the pkl files to extract.

Primary provenance:

- Mann, Gallen, and Clandinin, "Whole-Brain Calcium Imaging Reveals an Intrinsic Functional Network in Drosophila", Current Biology 2017.
- Local method note: `figure15_drosophila_final/DROSOPHILA_FC_ANALYSIS_METHOD.md`.
