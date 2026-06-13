"""Window sensitivity around the current high-pass plot settings.

Current base:
- C. elegans: high-pass 0.03 Hz, window 20, step 8
- Drosophila: high-pass 0.03 Hz, window 15, step 5
- Zebrafish: high-pass 0.03 Hz, window 20, step 5

Drosophila reuses existing high-pass window-sensitivity outputs. C. elegans
and zebrafish are recomputed with high-pass filtering for short/current/long
windows.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RAW_ROOT = ROOT / "ncomms_figures_rebuild/raw"
TAB_DIR = ROOT / "ncomms_figures_rebuild/tables"

HIGHPASS_HZ = 0.03
TOP_EDGE_FRACTION = 0.20
MAX_LAG = 30
MEASURES = ["EdgeStdFCV", "FCS", "ProfileCorrDistFCV", "FCSHalfLifeWindows", "ParticipationFlex"]

WINDOW_CONFIGS = {
    "short": {
        "C. elegans": (10, 4),
        "Zebrafish": (10, 3),
    },
    "current": {
        "C. elegans": (20, 8),
        "Zebrafish": (20, 5),
    },
    "long": {
        "C. elegans": (40, 16),
        "Zebrafish": (40, 10),
    },
}
DROSOPHILA_CONFIG_MAP = {"short": "half", "current": "baseline", "long": "double"}


def zscore_rows(traces: np.ndarray) -> np.ndarray:
    traces = np.asarray(traces, dtype=float)
    mean = np.nanmean(traces, axis=1, keepdims=True)
    sd = np.nanstd(traces, axis=1, keepdims=True)
    sd[sd <= 1e-12] = np.nan
    return np.nan_to_num((traces - mean) / sd, nan=0.0, posinf=0.0, neginf=0.0)


def highpass_filter(traces: np.ndarray, sampling_rate_hz: float) -> np.ndarray:
    from scipy.signal import butter, filtfilt

    nyq = 0.5 * sampling_rate_hz
    b, a = butter(2, HIGHPASS_HZ / nyq, btype="high")
    return filtfilt(b, a, traces, axis=1)


def iter_corr_windows(traces: np.ndarray, window: int, step: int, diagonal: float = np.nan):
    for start in range(0, traces.shape[1] - window + 1, step):
        corr = np.corrcoef(traces[:, start : start + window])
        np.fill_diagonal(corr, diagonal)
        yield corr.astype(np.float64, copy=False)


def row_corr_dist_excluding_diagonal(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    n_valid = n - 1
    sum_a = a.sum(axis=1)
    sum_b = b.sum(axis=1)
    sum_aa = np.square(a).sum(axis=1)
    sum_bb = np.square(b).sum(axis=1)
    sum_ab = (a * b).sum(axis=1)
    num = sum_ab - (sum_a * sum_b / n_valid)
    den_a = sum_aa - (sum_a * sum_a / n_valid)
    den_b = sum_bb - (sum_b * sum_b / n_valid)
    denom = np.sqrt(np.maximum(den_a, 0.0) * np.maximum(den_b, 0.0))
    corr = np.full(n, np.nan)
    ok = denom > 1e-12
    corr[ok] = np.clip(num[ok] / denom[ok], -1.0, 1.0)
    return 1.0 - corr


def fcv_and_profile_corr_distance(traces: np.ndarray, window: int, step: int) -> tuple[pd.DataFrame, int]:
    n_nodes = traces.shape[0]
    fc_sum = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    fc2_sum = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    n_windows = 0
    for corr in iter_corr_windows(traces, window, step, diagonal=0.0):
        fc_sum += corr
        fc2_sum += corr * corr
        n_windows += 1
    mean_fc = fc_sum / n_windows
    edge_std_mat = np.sqrt(np.maximum(fc2_sum / n_windows - mean_fc * mean_fc, 0.0))
    np.fill_diagonal(edge_std_mat, 0.0)
    edge_std = edge_std_mat.sum(axis=1) / max(n_nodes - 1, 1)

    corr_dist_sum = np.zeros(n_nodes, dtype=np.float64)
    corr_dist_count = np.zeros(n_nodes, dtype=np.int64)
    for corr in iter_corr_windows(traces, window, step, diagonal=0.0):
        cd = row_corr_dist_excluding_diagonal(corr, mean_fc)
        ok = np.isfinite(cd)
        corr_dist_sum[ok] += cd[ok]
        corr_dist_count[ok] += 1
    profile_corr = np.divide(
        corr_dist_sum,
        corr_dist_count,
        out=np.full(n_nodes, np.nan),
        where=corr_dist_count > 0,
    )
    return pd.DataFrame({"EdgeStdFCV": edge_std, "ProfileCorrDistFCV": profile_corr}), int(n_windows)


def fc_strength_series(traces: np.ndarray, window: int, step: int) -> np.ndarray:
    signed = []
    for corr in iter_corr_windows(traces, window, step, diagonal=np.nan):
        signed.append(np.nanmean(corr, axis=1))
    return np.vstack(signed)


def autocorr_at_lags(series: np.ndarray, max_lag: int) -> np.ndarray:
    n_time, n_nodes = series.shape
    max_lag = min(max_lag, n_time - 2)
    out = np.full((max_lag + 1, n_nodes), np.nan)
    out[0, :] = 1.0
    for lag in range(1, max_lag + 1):
        a = series[:-lag, :]
        b = series[lag:, :]
        for node in range(n_nodes):
            x = a[:, node]
            y = b[:, node]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= 4 and np.nanstd(x[ok]) > 1e-12 and np.nanstd(y[ok]) > 1e-12:
                out[lag, node] = np.corrcoef(x[ok], y[ok])[0, 1]
    return out


def exp_half_life(curve: np.ndarray) -> np.ndarray:
    n_lags, n_nodes = curve.shape
    out = np.full(n_nodes, np.nan)
    lags = np.arange(n_lags, dtype=float)
    for node in range(n_nodes):
        y = curve[:, node]
        ok = np.isfinite(y) & (y > 0) & (lags > 0)
        if ok.sum() < 3:
            continue
        x = lags[ok]
        slope = float(np.sum(x * np.log(y[ok])) / np.sum(x * x))
        if np.isfinite(slope) and slope < 0:
            out[node] = (-1.0 / slope) * np.log(2.0)
    return out


def modules_from_mean_fc(mean_fc: np.ndarray) -> np.ndarray:
    weights = np.abs(mean_fc).astype(float, copy=True)
    np.fill_diagonal(weights, 0.0)
    n = weights.shape[0]
    positive = weights[np.triu_indices(n, k=1)]
    positive = positive[positive > 0]
    if positive.size == 0:
        return np.arange(n, dtype=int)
    threshold = float(np.quantile(positive, 1.0 - TOP_EDGE_FRACTION))
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    rows, cols = np.where(np.triu(weights >= threshold, k=1))
    for i, j in zip(rows.tolist(), cols.tolist()):
        if weights[i, j] > 0:
            graph.add_edge(i, j, weight=float(weights[i, j]))
    communities = list(nx.algorithms.community.greedy_modularity_communities(graph, weight="weight"))
    labels = np.full(n, -1, dtype=int)
    for module_idx, community in enumerate(communities):
        for node in community:
            labels[int(node)] = module_idx
    for node in np.where(labels < 0)[0]:
        labels[node] = labels.max() + 1
    return labels


def participation_for_corr(corr: np.ndarray, labels: np.ndarray) -> np.ndarray:
    weights = np.abs(corr).astype(float, copy=True)
    np.fill_diagonal(weights, 0.0)
    total = weights.sum(axis=1)
    out = np.full(weights.shape[0], np.nan)
    valid = total > 1e-12
    frac_sq = np.zeros(weights.shape[0], dtype=float)
    for module in np.unique(labels):
        module_sum = weights[:, labels == module].sum(axis=1)
        frac = np.zeros_like(total)
        frac[valid] = module_sum[valid] / total[valid]
        frac_sq += frac * frac
    out[valid] = 1.0 - frac_sq[valid]
    return out


def participation_flex(traces: np.ndarray, window: int, step: int) -> pd.DataFrame:
    total = None
    n_windows = 0
    for corr in iter_corr_windows(traces, window, step, diagonal=0.0):
        corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
        if total is None:
            total = np.zeros_like(corr)
        total += corr
        n_windows += 1
    labels = modules_from_mean_fc(total / n_windows)
    values = []
    for corr in iter_corr_windows(traces, window, step, diagonal=0.0):
        values.append(participation_for_corr(np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0), labels))
    return pd.DataFrame({"ParticipationFlex": np.nanstd(np.vstack(values), axis=0)})


def measure_table_for_recording(traces: np.ndarray, window: int, step: int) -> tuple[pd.DataFrame, int]:
    geom, n_windows = fcv_and_profile_corr_distance(traces, window, step)
    signed = fc_strength_series(traces, window, step)
    fcs = pd.DataFrame(
        {
            "FCS": np.nanmean(signed, axis=0),
            "FCSHalfLifeWindows": exp_half_life(autocorr_at_lags(signed, MAX_LAG)),
        }
    )
    return pd.concat([geom, fcs, participation_flex(traces, window, step)], axis=1), n_windows


def aggregate_mean(df: pd.DataFrame, keys: list[str], extra: dict | None = None) -> pd.DataFrame:
    aggs = {measure: (measure, "mean") for measure in MEASURES}
    if extra:
        aggs = {**extra, **aggs}
    return df.groupby(keys, as_index=False).agg(**aggs)


def compute_celegans(config: str) -> tuple[pd.DataFrame, list[dict]]:
    window, step = WINDOW_CONFIGS[config]["C. elegans"]
    rows = []
    qcs = []
    for path in sorted((RAW_ROOT / "celegans_recording_pkl").glob("*_raw_traces.pkl")):
        with path.open("rb") as handle:
            raw = pickle.load(handle)
        traces = highpass_filter(np.asarray(raw["traces_spontaneous"], dtype=float), float(raw["sampling_rate_hz"]))
        traces = zscore_rows(traces)
        table, n_windows = measure_table_for_recording(traces, window, step)
        table.insert(0, "node", list(raw["neuron_names"]))
        rows.append(table)
        qcs.append({"species": "C. elegans", "window_config": config, "recording_id": raw["recording_id"], "window": window, "step": step, "n_nodes": traces.shape[0], "n_windows": n_windows})
    summary = aggregate_mean(pd.concat(rows, ignore_index=True), ["node"])
    summary.insert(0, "species", "C. elegans")
    summary["level"] = "neuron"
    return summary, qcs


def compute_zebrafish(config: str) -> tuple[pd.DataFrame, list[dict]]:
    window, step = WINDOW_CONFIGS[config]["Zebrafish"]
    rows = []
    qcs = []
    for path in sorted((RAW_ROOT / "zebrafish_recording_pkl").glob("*_raw_cluster_traces.pkl")):
        with path.open("rb") as handle:
            raw = pickle.load(handle)
        traces = highpass_filter(np.asarray(raw["traces"], dtype=float), float(raw["sampling_rate_hz"]))
        traces = zscore_rows(traces)
        table, n_windows = measure_table_for_recording(traces, window, step)
        table.insert(0, "root_area_name", list(raw["root_area_names"]))
        table.insert(0, "root_area_id", list(raw["root_area_ids"]))
        table.insert(0, "cluster_id", list(raw["cluster_ids"]))
        table.insert(0, "recording_id", raw["recording_id"])
        rows.append(table)
        qcs.append({"species": "Zebrafish", "window_config": config, "recording_id": raw["recording_id"], "window": window, "step": step, "n_nodes": traces.shape[0], "n_windows": n_windows})
    subject_region = aggregate_mean(
        pd.concat(rows, ignore_index=True),
        ["recording_id", "root_area_id", "root_area_name"],
        extra={"n_clusters": ("cluster_id", "size")},
    )
    summary = aggregate_mean(subject_region, ["root_area_id", "root_area_name"])
    summary = summary.rename(columns={"root_area_name": "node"})
    summary.insert(0, "species", "Zebrafish")
    summary["level"] = "root_area"
    return summary, qcs


def load_drosophila_from_existing(config: str) -> pd.DataFrame:
    current = TAB_DIR / "current_base_window_sensitivity_node_summary.csv"
    if current.exists():
        df = pd.read_csv(current)
        df = df[df["species"].eq("Drosophila") & df["window_config"].eq(config)]
        if not df.empty:
            return df[["species", "node", "EdgeStdFCV", "FCS", "ProfileCorrDistFCV", "FCSHalfLifeWindows", "ParticipationFlex", "level"]].copy()

    old = DROSOPHILA_CONFIG_MAP[config]
    core = pd.read_csv(TAB_DIR / "window_duration_sensitivity_core4_node_summary.csv")
    core = core[core["species"].eq("Drosophila") & core["window_config"].eq(old)]
    core = core[["species", "node", "EdgeStdFCV", "FCS", "ProfileCorrDistFCV", "FCSHalfLifeWindows", "level"]]
    part = pd.read_csv(TAB_DIR / "window_duration_sensitivity_node_summary.csv")
    part = part[part["species"].eq("Drosophila") & part["window_config"].eq(old)]
    part = part[["species", "node", "ParticipationFlex"]]
    return core.merge(part, on=["species", "node"], how="left")


def vs_fcv_correlations(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (species, config), sdf in summary.groupby(["species", "window_config"], sort=False):
        for measure in ["FCS", "ProfileCorrDistFCV", "FCSHalfLifeWindows", "ParticipationFlex"]:
            sub = sdf[["EdgeStdFCV", measure]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(sub) >= 3:
                pr, pp = stats.pearsonr(sub[measure], sub["EdgeStdFCV"])
                sr, sp = stats.spearmanr(sub[measure], sub["EdgeStdFCV"])
            else:
                pr = pp = sr = sp = np.nan
            rows.append({"species": species, "window_config": config, "measure": measure, "n": len(sub), "pearson_r_vs_fcv": pr, "pearson_p": pp, "spearman_r_vs_fcv": sr, "spearman_p": sp})
    return pd.DataFrame(rows)


def baseline_stability(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    current = summary[summary["window_config"].eq("current")]
    for species, base in current.groupby("species", sort=False):
        for config in ["short", "long"]:
            other = summary[summary["species"].eq(species) & summary["window_config"].eq(config)]
            merged = base.merge(other, on=["species", "node"], suffixes=("_current", f"_{config}"))
            for measure in MEASURES:
                cols = [f"{measure}_current", f"{measure}_{config}"]
                sub = merged[cols].replace([np.inf, -np.inf], np.nan).dropna()
                if len(sub) >= 3:
                    pr, pp = stats.pearsonr(sub[cols[0]], sub[cols[1]])
                    sr, sp = stats.spearmanr(sub[cols[0]], sub[cols[1]])
                    mad = float(np.nanmean(np.abs(sub[cols[1]] - sub[cols[0]])))
                else:
                    pr = pp = sr = sp = mad = np.nan
                rows.append({"species": species, "compare_window": config, "measure": measure, "n": len(sub), "pearson_r_vs_current": pr, "pearson_p": pp, "spearman_r_vs_current": sr, "spearman_p": sp, "mean_abs_delta_vs_current": mad})
    return pd.DataFrame(rows)


def main() -> None:
    all_summary = []
    all_qc = []
    for config in ["short", "current", "long"]:
        print(f"[config] {config}", flush=True)
        ce, qc = compute_celegans(config)
        all_summary.append(ce.assign(window_config=config, window=WINDOW_CONFIGS[config]["C. elegans"][0], step=WINDOW_CONFIGS[config]["C. elegans"][1], highpass_hz=HIGHPASS_HZ))
        all_qc.extend(qc)
        fly = load_drosophila_from_existing(config)
        fly_window = {"short": (8, 3), "current": (15, 5), "long": (30, 10)}[config]
        all_summary.append(fly.assign(window_config=config, window=fly_window[0], step=fly_window[1], highpass_hz=HIGHPASS_HZ))
        zf, qc = compute_zebrafish(config)
        all_summary.append(zf.assign(window_config=config, window=WINDOW_CONFIGS[config]["Zebrafish"][0], step=WINDOW_CONFIGS[config]["Zebrafish"][1], highpass_hz=HIGHPASS_HZ))
        all_qc.extend(qc)
    summary = pd.concat(all_summary, ignore_index=True)
    vs_fcv = vs_fcv_correlations(summary)
    stability = baseline_stability(summary)
    summary.to_csv(TAB_DIR / "current_base_window_sensitivity_node_summary.csv", index=False)
    vs_fcv.to_csv(TAB_DIR / "current_base_window_sensitivity_vs_fcv_correlations.csv", index=False)
    stability.to_csv(TAB_DIR / "current_base_window_sensitivity_baseline_stability.csv", index=False)
    pd.DataFrame(all_qc).to_csv(TAB_DIR / "current_base_window_sensitivity_qc.csv", index=False)
    print("wrote current_base_window_sensitivity_node_summary.csv")
    print(vs_fcv.to_string(index=False))
    print(stability.to_string(index=False))


if __name__ == "__main__":
    main()
