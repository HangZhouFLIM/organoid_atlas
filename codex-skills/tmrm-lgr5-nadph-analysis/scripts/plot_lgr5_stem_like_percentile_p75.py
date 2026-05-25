#!/usr/bin/env python
"""Stem-like LGR5-GFP/TMRM/NADPH percentile-p75 figures."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


METHOD = "percentile_p75"
ALPHA = 0.05
DROPLET_DENSITY_SCALE = 10000

CONDITION_ORDER = ["control", "rotenone", "oligomycin"]
CONDITION_LABELS = {
    "control": "Control",
    "rotenone": "Rotenone",
    "oligomycin": "Oligomycin",
}
CONDITION_COLORS = {
    "control": "#6F7378",
    "rotenone": "#D8893A",
    "oligomycin": "#3B78B8",
}

GFP_ORDER = ["low", "medium", "high"]
GFP_LABELS = {
    "low": "GFP low\n(<=p25)",
    "medium": "GFP medium\n(p25-p75)",
    "high": "GFP high\n(clean p75)",
}
GFP_COLORS = {
    "low": "#5B7DB8",
    "medium": "#4F9A5D",
    "high": "#D8893A",
}

STEM_ORDER = ["gfp_low", "stem_like"]
STEM_LABELS = {
    "gfp_low": "GFP low\n(<=p25)",
    "stem_like": "Stem-like GFP\n(medium+high)",
}
STEM_COLORS = {
    "gfp_low": "#5B7DB8",
    "stem_like": "#4F9A5D",
}

METRIC_ORDER = ["Int", "TauPhase", "TauModulation", "TauNormalized"]
METRIC_LABELS = {
    "Int": "NADPH intensity",
    "TauPhase": "Tau phase",
    "TauModulation": "Tau modulation",
    "TauNormalized": "Tau normalized",
}
DELTA_LABELS = {
    "Int": "High TMRM - low TMRM\nNADPH intensity",
    "TauPhase": "High TMRM - low TMRM\nTau phase (ns)",
    "TauModulation": "High TMRM - low TMRM\nTau modulation (ns)",
    "TauNormalized": "High TMRM - low TMRM\nTau normalized (ns)",
}
DROPLET_YLABELS = {
    "Int": "Droplet NADPH intensity",
    "TauPhase": "Droplet tau phase (ns)",
    "TauModulation": "Droplet tau modulation (ns)",
    "TauNormalized": "Droplet tau normalized (ns)",
}
NADPH_REGION_YLABELS = {
    "Int": "Real NADPH intensity",
    "TauPhase": "Real NADPH tau phase (ns)",
    "TauModulation": "Real NADPH tau modulation (ns)",
    "TauNormalized": "Real NADPH tau normalized (ns)",
}
TMRM_FEATURE_ORDER = ["tmrm_high_fraction", "Int", "TauPhase", "TauModulation", "TauNormalized"]
TMRM_FEATURE_LABELS = {
    "tmrm_high_fraction": "TMRM-high fraction",
    "Int": "TMRM intensity",
    "TauPhase": "TMRM tau phase",
    "TauModulation": "TMRM tau modulation",
    "TauNormalized": "TMRM tau normalized",
}
TMRM_FEATURE_YLABELS = {
    "tmrm_high_fraction": "TMRM-high fraction",
    "Int": "TMRM intensity",
    "TauPhase": "TMRM tau phase (ns)",
    "TauModulation": "TMRM tau modulation (ns)",
    "TauNormalized": "TMRM tau normalized (ns)",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create stem-like LGR5-GFP/TMRM/NADPH percentile-p75 figures."
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "analysis_TMRM_LGR5_NADPH",
        help="Analysis folder containing final_*_percentile_p75.csv outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output folder. Defaults to ANALYSIS_DIR/figures_lgr5_stem_like_percentile_p75.",
    )
    return parser.parse_args()


def apply_style() -> None:
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


def mean_ci(values: pd.Series | np.ndarray) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return math.nan, math.nan
    mean = float(np.mean(arr))
    if arr.size < 2:
        return mean, math.nan
    ci = float(stats.t.ppf(0.975, arr.size - 1) * stats.sem(arr))
    return mean, ci


def paired_test(a: np.ndarray, b: np.ndarray) -> dict[str, object]:
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    diff = diff[np.isfinite(diff)]
    if diff.size < 2:
        return {"n": int(diff.size), "test": "NA", "statistic": math.nan, "p_value": math.nan}
    if diff.size >= 3 and np.ptp(diff) > 0:
        normal_p = float(stats.shapiro(diff).pvalue)
    else:
        normal_p = math.nan
    if np.isfinite(normal_p) and normal_p > ALPHA:
        result = stats.ttest_1samp(diff, 0.0)
        return {
            "n": int(diff.size),
            "test": "paired t-test",
            "statistic": float(result.statistic),
            "p_value": float(result.pvalue),
            "normality_p": normal_p,
        }
    if np.any(np.abs(diff) > 0):
        result = stats.wilcoxon(diff)
        return {
            "n": int(diff.size),
            "test": "Wilcoxon signed-rank",
            "statistic": float(result.statistic),
            "p_value": float(result.pvalue),
            "normality_p": normal_p,
        }
    return {
        "n": int(diff.size),
        "test": "all paired differences are zero",
        "statistic": math.nan,
        "p_value": 1.0,
        "normality_p": normal_p,
    }


def one_sample_test(values: np.ndarray) -> dict[str, object]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        return {"n": int(arr.size), "test": "NA", "statistic": math.nan, "p_value": math.nan}
    if arr.size >= 3 and np.ptp(arr) > 0:
        normal_p = float(stats.shapiro(arr).pvalue)
    else:
        normal_p = math.nan
    if np.isfinite(normal_p) and normal_p > ALPHA:
        result = stats.ttest_1samp(arr, 0.0)
        return {
            "n": int(arr.size),
            "test": "one-sample t-test vs 0",
            "statistic": float(result.statistic),
            "p_value": float(result.pvalue),
            "normality_p": normal_p,
        }
    if np.any(np.abs(arr) > 0):
        result = stats.wilcoxon(arr)
        return {
            "n": int(arr.size),
            "test": "Wilcoxon signed-rank vs 0",
            "statistic": float(result.statistic),
            "p_value": float(result.pvalue),
            "normality_p": normal_p,
        }
    return {
        "n": int(arr.size),
        "test": "all values are zero",
        "statistic": math.nan,
        "p_value": 1.0,
        "normality_p": normal_p,
    }


def condition_vs_control(values: pd.DataFrame, value_col: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    control = values[values["condition"].astype(str) == "control"][value_col].to_numpy(dtype=float)
    for condition in CONDITION_ORDER:
        if condition == "control":
            continue
        treatment = values[values["condition"].astype(str) == condition][value_col].to_numpy(dtype=float)
        if len(control) < 2 or len(treatment) < 2:
            rows.append(
                {
                    "comparison": f"{condition} vs control",
                    "test": "NA",
                    "statistic": math.nan,
                    "p_value": math.nan,
                    "effect_treatment_minus_control": math.nan,
                }
            )
            continue
        result = stats.mannwhitneyu(treatment, control, alternative="two-sided")
        rows.append(
            {
                "comparison": f"{condition} vs control",
                "test": "Mann-Whitney U",
                "statistic": float(result.statistic),
                "p_value": float(result.pvalue),
                "effect_treatment_minus_control": float(np.nanmean(treatment) - np.nanmean(control)),
            }
        )
    return pd.DataFrame(rows)


def parse_compartment(compartment: str) -> tuple[str, str]:
    if compartment.startswith("high_tmrm_lgr5_"):
        return "high_tmrm", compartment.replace("high_tmrm_lgr5_", "")
    if compartment.startswith("low_tmrm_lgr5_"):
        return "low_tmrm", compartment.replace("low_tmrm_lgr5_", "")
    raise ValueError(f"Unexpected compartment: {compartment}")


def load_tables(analysis_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    region = pd.read_csv(analysis_dir / f"final_nadph_compartment_summary_{METHOD}.csv")
    gfp_qc = pd.read_csv(analysis_dir / "lgr5_gfp_mask_qc.csv")
    tmrm_gfp = pd.read_csv(analysis_dir / "tmrm_gfp_group_summary.csv")
    lipid_samples = pd.read_csv(analysis_dir / "oxidized_lipid_droplet_sample_summary.csv")
    lipid_objects = pd.read_csv(analysis_dir / "oxidized_lipid_droplet_objects.csv")
    return region, gfp_qc, tmrm_gfp, lipid_samples, lipid_objects


def build_three_group_tmrm_allocation(region: pd.DataFrame) -> pd.DataFrame:
    base = region[
        (region["method"] == METHOD)
        & (region["nadph_class"] == "real_nadph")
        & (region["metric"] == "Int")
    ].copy()
    parsed = base["compartment"].map(parse_compartment)
    base["tmrm_region"] = parsed.map(lambda item: item[0])
    base["gfp_group"] = parsed.map(lambda item: item[1])
    grouped = base.pivot_table(
        index=["condition", "sample", "gfp_group"],
        columns="tmrm_region",
        values="compartment_pixels",
        aggfunc="first",
    ).reset_index()
    grouped["high_tmrm"] = grouped["high_tmrm"].fillna(0)
    grouped["low_tmrm"] = grouped["low_tmrm"].fillna(0)
    grouped["total_pixels"] = grouped["high_tmrm"] + grouped["low_tmrm"]
    grouped["tmrm_high_fraction"] = grouped["high_tmrm"] / grouped["total_pixels"]
    grouped["gfp_group"] = pd.Categorical(grouped["gfp_group"], GFP_ORDER, ordered=True)
    grouped["condition"] = pd.Categorical(grouped["condition"], CONDITION_ORDER, ordered=True)
    return grouped.sort_values(["condition", "sample", "gfp_group"]).reset_index(drop=True)


def build_stem_tmrm_allocation(region: pd.DataFrame) -> pd.DataFrame:
    three = build_three_group_tmrm_allocation(region)
    low = three[three["gfp_group"].astype(str) == "low"].copy()
    low["gfp_group"] = "gfp_low"
    stem = (
        three[three["gfp_group"].astype(str).isin(["medium", "high"])]
        .groupby(["condition", "sample"], as_index=False, observed=False)
        .agg(
            high_tmrm=("high_tmrm", "sum"),
            low_tmrm=("low_tmrm", "sum"),
            total_pixels=("total_pixels", "sum"),
        )
    )
    stem["tmrm_high_fraction"] = stem["high_tmrm"] / stem["total_pixels"]
    stem["gfp_group"] = "stem_like"
    grouped = pd.concat([low, stem], ignore_index=True, sort=False)
    grouped["gfp_group"] = pd.Categorical(grouped["gfp_group"], STEM_ORDER, ordered=True)
    grouped["condition"] = pd.Categorical(grouped["condition"], CONDITION_ORDER, ordered=True)
    return grouped.sort_values(["condition", "sample", "gfp_group"]).reset_index(drop=True)


def build_three_group_nadph_delta(region: pd.DataFrame) -> pd.DataFrame:
    base = region[
        (region["method"] == METHOD)
        & (region["nadph_class"] == "real_nadph")
        & (region["metric"].isin(METRIC_ORDER))
    ].copy()
    parsed = base["compartment"].map(parse_compartment)
    base["tmrm_region"] = parsed.map(lambda item: item[0])
    base["gfp_group"] = parsed.map(lambda item: item[1])
    pivot = base.pivot_table(
        index=["condition", "sample", "gfp_group", "metric"],
        columns="tmrm_region",
        values=["mean", "n_pixels"],
        aggfunc="first",
    )
    pivot.columns = [f"{value}_{region_name}" for value, region_name in pivot.columns]
    pivot = pivot.reset_index()
    pivot["delta_high_minus_low"] = pivot["mean_high_tmrm"] - pivot["mean_low_tmrm"]
    pivot["condition"] = pd.Categorical(pivot["condition"], CONDITION_ORDER, ordered=True)
    pivot["gfp_group"] = pd.Categorical(pivot["gfp_group"], GFP_ORDER, ordered=True)
    pivot["metric"] = pd.Categorical(pivot["metric"], METRIC_ORDER, ordered=True)
    return pivot.sort_values(["metric", "condition", "sample", "gfp_group"]).reset_index(drop=True)


def build_stem_nadph_delta(region: pd.DataFrame) -> pd.DataFrame:
    base = region[
        (region["method"] == METHOD)
        & (region["nadph_class"] == "real_nadph")
        & (region["metric"].isin(METRIC_ORDER))
    ].copy()
    parsed = base["compartment"].map(parse_compartment)
    base["tmrm_region"] = parsed.map(lambda item: item[0])
    base["source_gfp_group"] = parsed.map(lambda item: item[1])
    base["gfp_group"] = base["source_gfp_group"].map(
        {"low": "gfp_low", "medium": "stem_like", "high": "stem_like"}
    )

    rows: list[dict[str, object]] = []
    group_cols = ["condition", "sample", "gfp_group", "metric", "tmrm_region"]
    for keys, group in base.groupby(group_cols, observed=False, sort=True):
        condition, sample, gfp_group, metric, tmrm_region = keys
        weights = group["n_pixels"].to_numpy(dtype=float)
        means = group["mean"].to_numpy(dtype=float)
        valid = np.isfinite(weights) & np.isfinite(means) & (weights > 0)
        total_n = float(weights[valid].sum())
        weighted_mean = (
            float(np.average(means[valid], weights=weights[valid]))
            if total_n > 0
            else math.nan
        )
        rows.append(
            {
                "condition": condition,
                "sample": sample,
                "gfp_group": gfp_group,
                "metric": metric,
                "tmrm_region": tmrm_region,
                "mean": weighted_mean,
                "n_pixels": int(total_n),
            }
        )
    weighted = pd.DataFrame(rows)
    pivot = weighted.pivot_table(
        index=["condition", "sample", "gfp_group", "metric"],
        columns="tmrm_region",
        values=["mean", "n_pixels"],
        aggfunc="first",
    )
    pivot.columns = [f"{value}_{region_name}" for value, region_name in pivot.columns]
    pivot = pivot.reset_index()
    pivot["delta_high_minus_low"] = pivot["mean_high_tmrm"] - pivot["mean_low_tmrm"]
    pivot["condition"] = pd.Categorical(pivot["condition"], CONDITION_ORDER, ordered=True)
    pivot["gfp_group"] = pd.Categorical(pivot["gfp_group"], STEM_ORDER, ordered=True)
    pivot["metric"] = pd.Categorical(pivot["metric"], METRIC_ORDER, ordered=True)
    return pivot.sort_values(["metric", "condition", "sample", "gfp_group"]).reset_index(drop=True)


def build_nadph_high_low_long(delta: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for tmrm_region, mean_col, n_col in [
        ("low_tmrm", "mean_low_tmrm", "n_pixels_low_tmrm"),
        ("high_tmrm", "mean_high_tmrm", "n_pixels_high_tmrm"),
    ]:
        frame = delta[
            ["condition", "sample", "gfp_group", "metric", mean_col, n_col]
        ].rename(columns={mean_col: "mean", n_col: "n_pixels"})
        frame["tmrm_region"] = tmrm_region
        rows.append(frame)
    paired = pd.concat(rows, ignore_index=True)
    paired["condition"] = pd.Categorical(paired["condition"], CONDITION_ORDER, ordered=True)
    paired["gfp_group"] = pd.Categorical(paired["gfp_group"], STEM_ORDER, ordered=True)
    paired["metric"] = pd.Categorical(paired["metric"], METRIC_ORDER, ordered=True)
    paired["tmrm_region"] = pd.Categorical(
        paired["tmrm_region"], ["low_tmrm", "high_tmrm"], ordered=True
    )
    return paired.sort_values(
        ["metric", "condition", "sample", "gfp_group", "tmrm_region"]
    ).reset_index(drop=True)


def build_lipid_droplet_density(lipid_samples: pd.DataFrame) -> pd.DataFrame:
    lipid = lipid_samples.copy()
    density_col = "lipid_droplet_density_per_10000_valid_pixels"
    if density_col not in lipid.columns:
        lipid[density_col] = (
            lipid["n_lipid_droplets"] / lipid["valid_common_pixels"] * DROPLET_DENSITY_SCALE
        )
    lipid["condition"] = pd.Categorical(lipid["condition"], CONDITION_ORDER, ordered=True)
    return lipid.sort_values(["condition", "sample"]).reset_index(drop=True)


def build_lipid_droplet_metrics(lipid_objects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[pd.DataFrame] = []
    for metric in METRIC_ORDER:
        value_col = f"mean_{metric}"
        if value_col not in lipid_objects.columns:
            continue
        metric_df = lipid_objects[
            ["condition", "sample", "droplet_id", "area_pixels", value_col]
        ].rename(columns={value_col: "value"})
        metric_df["metric"] = metric
        rows.append(metric_df)
    if rows:
        droplet_metrics = pd.concat(rows, ignore_index=True)
    else:
        droplet_metrics = pd.DataFrame(
            columns=["condition", "sample", "droplet_id", "area_pixels", "value", "metric"]
        )
    droplet_metrics["condition"] = pd.Categorical(
        droplet_metrics["condition"], CONDITION_ORDER, ordered=True
    )
    droplet_metrics["metric"] = pd.Categorical(
        droplet_metrics["metric"], METRIC_ORDER, ordered=True
    )
    sample_means = (
        droplet_metrics.groupby(["condition", "sample", "metric"], observed=True, as_index=False)
        .agg(
            value=("value", "mean"),
            n_droplets=("droplet_id", "count"),
            total_droplet_area_pixels=("area_pixels", "sum"),
            mean_droplet_area_pixels=("area_pixels", "mean"),
        )
    )
    sample_means["condition"] = pd.Categorical(
        sample_means["condition"], CONDITION_ORDER, ordered=True
    )
    sample_means["metric"] = pd.Categorical(sample_means["metric"], METRIC_ORDER, ordered=True)
    return (
        droplet_metrics.sort_values(["metric", "condition", "sample", "droplet_id"]).reset_index(drop=True),
        sample_means.sort_values(["metric", "condition", "sample"]).reset_index(drop=True),
    )


def build_tmrm_features(tmrm_stem: pd.DataFrame, tmrm_gfp: pd.DataFrame) -> pd.DataFrame:
    fraction = tmrm_stem[
        ["condition", "sample", "gfp_group", "total_pixels", "tmrm_high_fraction"]
    ].rename(
        columns={
            "total_pixels": "n_pixels",
            "tmrm_high_fraction": "value",
        }
    )
    fraction["feature"] = "tmrm_high_fraction"
    fraction["metric"] = "tmrm_high_fraction"

    metric_rows = tmrm_gfp[
        tmrm_gfp["gfp_group"].astype(str).isin(["low", "stem_like"])
        & tmrm_gfp["metric"].astype(str).isin(METRIC_ORDER)
    ].copy()
    metric_rows["gfp_group"] = metric_rows["gfp_group"].map(
        {"low": "gfp_low", "stem_like": "stem_like"}
    )
    metric_rows = metric_rows.rename(columns={"mean": "value"})
    metric_rows["feature"] = metric_rows["metric"].astype(str)
    metric_rows = metric_rows[
        ["condition", "sample", "gfp_group", "metric", "feature", "n_pixels", "value"]
    ]

    feature_df = pd.concat(
        [
            fraction[["condition", "sample", "gfp_group", "metric", "feature", "n_pixels", "value"]],
            metric_rows,
        ],
        ignore_index=True,
    )
    feature_df["condition"] = pd.Categorical(
        feature_df["condition"], CONDITION_ORDER, ordered=True
    )
    feature_df["gfp_group"] = pd.Categorical(feature_df["gfp_group"], STEM_ORDER, ordered=True)
    feature_df["feature"] = pd.Categorical(
        feature_df["feature"], TMRM_FEATURE_ORDER, ordered=True
    )
    return feature_df.sort_values(["feature", "condition", "sample", "gfp_group"]).reset_index(drop=True)


def summarize_groups(df: pd.DataFrame, group_cols: list[str], value_col: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for keys, group in df.groupby(group_cols, observed=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values = group[value_col].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        mean, ci = mean_ci(values)
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n": int(finite.size),
                "mean": mean,
                "ci95": ci,
                "median": float(np.nanmedian(finite)) if finite.size else math.nan,
                "sd": float(np.nanstd(finite, ddof=1)) if finite.size > 1 else math.nan,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def tmrm_stem_stats(tmrm: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition, group in tmrm.groupby("condition", observed=False):
        wide = group.pivot(index="sample", columns="gfp_group", values="tmrm_high_fraction")
        paired = wide[["stem_like", "gfp_low"]].dropna()
        result = paired_test(
            paired["stem_like"].to_numpy(dtype=float),
            paired["gfp_low"].to_numpy(dtype=float),
        )
        row = {
            "analysis_level": "within_condition_pairwise",
            "condition": condition,
            "comparison": "stem_like - gfp_low",
        }
        row.update(result)
        rows.append(row)
    return pd.DataFrame(rows)


def tmrm_feature_stats(features: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (condition, feature), group in features.groupby(["condition", "feature"], observed=False):
        wide = group.pivot(index="sample", columns="gfp_group", values="value")
        if not all(group_name in wide for group_name in STEM_ORDER):
            continue
        paired = wide[["stem_like", "gfp_low"]].dropna()
        result = paired_test(
            paired["stem_like"].to_numpy(dtype=float),
            paired["gfp_low"].to_numpy(dtype=float),
        )
        row = {
            "analysis_level": "within_condition_pairwise",
            "condition": condition,
            "feature": feature,
            "comparison": "stem_like - gfp_low",
        }
        row.update(result)
        rows.append(row)
    return pd.DataFrame(rows)


def nadph_delta_stats(delta: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (condition, gfp_group, metric), group in delta.groupby(
        ["condition", "gfp_group", "metric"], observed=False
    ):
        result = one_sample_test(group["delta_high_minus_low"].to_numpy(dtype=float))
        row = {
            "analysis_level": "delta_vs_zero",
            "condition": condition,
            "gfp_group": gfp_group,
            "metric": metric,
            "comparison": "high TMRM - low TMRM vs 0",
        }
        row.update(result)
        rows.append(row)
    return pd.DataFrame(rows)


def nadph_stem_group_stats(delta: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (condition, metric), group in delta.groupby(["condition", "metric"], observed=False):
        wide = group.pivot(index="sample", columns="gfp_group", values="delta_high_minus_low")
        paired = wide[["stem_like", "gfp_low"]].dropna()
        result = paired_test(
            paired["stem_like"].to_numpy(dtype=float),
            paired["gfp_low"].to_numpy(dtype=float),
        )
        row = {
            "analysis_level": "within_condition_group_pairwise",
            "condition": condition,
            "metric": metric,
            "comparison": "stem_like delta - gfp_low delta",
        }
        row.update(result)
        rows.append(row)
    return pd.DataFrame(rows)


def lipid_density_stats(lipid: pd.DataFrame) -> pd.DataFrame:
    value_col = "lipid_droplet_density_per_10000_valid_pixels"
    rows: list[dict[str, object]] = []
    summary = summarize_groups(lipid, ["condition"], value_col)
    for _, row in summary.iterrows():
        rows.append(
            {
                "analysis_level": "condition_summary",
                "condition": row["condition"],
                "comparison": "condition mean",
                "test": "summary",
                "statistic": row["mean"],
                "p_value": math.nan,
                "n_samples": row["n"],
                "ci95": row["ci95"],
            }
        )
    arrays = [
        lipid[lipid["condition"].astype(str) == condition][value_col].to_numpy(dtype=float)
        for condition in CONDITION_ORDER
    ]
    arrays = [arr[np.isfinite(arr)] for arr in arrays if len(arr[np.isfinite(arr)]) > 0]
    if len(arrays) >= 2:
        result = stats.kruskal(*arrays)
        rows.append(
            {
                "analysis_level": "between_condition_omnibus",
                "condition": "ALL",
                "comparison": "condition effect",
                "test": "Kruskal-Wallis",
                "statistic": float(result.statistic),
                "p_value": float(result.pvalue),
                "n_samples": int(sum(len(arr) for arr in arrays)),
                "ci95": math.nan,
            }
        )
    posthoc = condition_vs_control(lipid, value_col)
    for _, row in posthoc.iterrows():
        out = {
            "analysis_level": "posthoc_vs_control",
            "condition": row["comparison"].replace(" vs control", ""),
            "comparison": row["comparison"],
            "test": row["test"],
            "statistic": row["statistic"],
            "p_value": row["p_value"],
            "n_samples": math.nan,
            "ci95": math.nan,
            "effect_treatment_minus_control": row["effect_treatment_minus_control"],
        }
        rows.append(out)
    return pd.DataFrame(rows)


def lipid_droplet_metric_stats(sample_means: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for metric, metric_df in sample_means.groupby("metric", observed=False):
        summary = summarize_groups(metric_df, ["condition"], "value")
        for _, row in summary.iterrows():
            rows.append(
                {
                    "analysis_level": "condition_summary",
                    "condition": row["condition"],
                    "metric": metric,
                    "comparison": "condition mean of sample-level droplet means",
                    "test": "summary",
                    "statistic": row["mean"],
                    "p_value": math.nan,
                    "n_samples": row["n"],
                    "ci95": row["ci95"],
                }
            )
        arrays = [
            metric_df[metric_df["condition"].astype(str) == condition]["value"].to_numpy(dtype=float)
            for condition in CONDITION_ORDER
        ]
        arrays = [arr[np.isfinite(arr)] for arr in arrays if len(arr[np.isfinite(arr)]) > 0]
        if len(arrays) >= 2:
            result = stats.kruskal(*arrays)
            rows.append(
                {
                    "analysis_level": "between_condition_omnibus",
                    "condition": "ALL",
                    "metric": metric,
                    "comparison": "condition effect",
                    "test": "Kruskal-Wallis",
                    "statistic": float(result.statistic),
                    "p_value": float(result.pvalue),
                    "n_samples": int(sum(len(arr) for arr in arrays)),
                    "ci95": math.nan,
                }
            )
        posthoc = condition_vs_control(metric_df, "value")
        for _, row in posthoc.iterrows():
            rows.append(
                {
                    "analysis_level": "posthoc_vs_control",
                    "condition": row["comparison"].replace(" vs control", ""),
                    "metric": metric,
                    "comparison": row["comparison"],
                    "test": row["test"],
                    "statistic": row["statistic"],
                    "p_value": row["p_value"],
                    "n_samples": math.nan,
                    "ci95": math.nan,
                    "effect_treatment_minus_control": row["effect_treatment_minus_control"],
                }
            )
    return pd.DataFrame(rows)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=8,
        va="bottom",
        ha="left",
    )


def draw_mean_ci(ax: plt.Axes, x: float, values: np.ndarray, color: str) -> None:
    mean, ci = mean_ci(values)
    if not np.isfinite(mean):
        return
    ax.scatter([x], [mean], s=25, color=color, edgecolors="black", linewidths=0.45, zorder=5)
    if np.isfinite(ci):
        ax.errorbar(
            [x],
            [mean],
            yerr=[[ci], [ci]],
            color="black",
            linewidth=0.8,
            capsize=2.2,
            zorder=4,
        )


def plot_tmrm_allocation_stem(tmrm: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.45), sharey=True, constrained_layout=True)
    x_positions = {group: idx for idx, group in enumerate(STEM_ORDER)}
    for ax, condition, panel_label in zip(axes, CONDITION_ORDER, ["a", "b", "c"]):
        subset = tmrm[tmrm["condition"].astype(str) == condition]
        wide = subset.pivot(index="sample", columns="gfp_group", values="tmrm_high_fraction")
        for _, row in wide.iterrows():
            if all(group in row.index and np.isfinite(row[group]) for group in STEM_ORDER):
                ax.plot(
                    [x_positions[group] for group in STEM_ORDER],
                    [row[group] for group in STEM_ORDER],
                    color="#B8B8B8",
                    linewidth=0.55,
                    alpha=0.75,
                    zorder=1,
                )
        for group in STEM_ORDER:
            values = subset[subset["gfp_group"].astype(str) == group]["tmrm_high_fraction"].to_numpy(dtype=float)
            x = x_positions[group]
            jitter = np.linspace(-0.035, 0.035, max(len(values), 1))[: len(values)]
            ax.scatter(
                np.full(len(values), x) + jitter,
                values,
                s=17,
                color=STEM_COLORS[group],
                edgecolors="black",
                linewidths=0.35,
                alpha=0.92,
                zorder=3,
            )
            draw_mean_ci(ax, x, values, STEM_COLORS[group])
        ax.set_title(CONDITION_LABELS[condition], fontsize=7.7, fontweight="bold")
        ax.set_xticks([x_positions[group] for group in STEM_ORDER])
        ax.set_xticklabels([STEM_LABELS[group] for group in STEM_ORDER], rotation=0)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("TMRM-high fraction" if ax is axes[0] else "")
        add_panel_label(ax, panel_label)
    save_figure(fig, output_dir / "tmrm_high_fraction_stem_like_percentile_p75")
    plt.close(fig)


def plot_tmrm_features(features: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(7.3, 5.15), constrained_layout=True)
    axes_flat = axes.ravel()
    x_base = {group: idx for idx, group in enumerate(STEM_ORDER)}
    offsets = {"control": -0.22, "rotenone": 0.0, "oligomycin": 0.22}

    for ax, feature, panel_label in zip(axes_flat, TMRM_FEATURE_ORDER, ["a", "b", "c", "d", "e"]):
        feature_df = features[features["feature"].astype(str) == feature]
        for condition in CONDITION_ORDER:
            condition_df = feature_df[feature_df["condition"].astype(str) == condition]
            wide = condition_df.pivot(index="sample", columns="gfp_group", values="value")
            for _, row in wide.iterrows():
                if all(group in row.index and np.isfinite(row[group]) for group in STEM_ORDER):
                    ax.plot(
                        [x_base[group] + offsets[condition] for group in STEM_ORDER],
                        [row[group] for group in STEM_ORDER],
                        color=CONDITION_COLORS[condition],
                        linewidth=0.45,
                        alpha=0.22,
                        zorder=1,
                    )
            for group in STEM_ORDER:
                values = condition_df[condition_df["gfp_group"].astype(str) == group][
                    "value"
                ].to_numpy(dtype=float)
                x = x_base[group] + offsets[condition]
                jitter = np.linspace(-0.035, 0.035, max(len(values), 1))[: len(values)]
                ax.scatter(
                    np.full(len(values), x) + jitter,
                    values,
                    s=16,
                    color=CONDITION_COLORS[condition],
                    edgecolors="black",
                    linewidths=0.3,
                    alpha=0.9,
                    zorder=3,
                    label=CONDITION_LABELS[condition] if group == STEM_ORDER[0] and feature == TMRM_FEATURE_ORDER[0] else None,
                )
                draw_mean_ci(ax, x, values, CONDITION_COLORS[condition])
        ax.set_title(TMRM_FEATURE_LABELS[feature], fontsize=7.7, fontweight="bold")
        ax.set_xticks([x_base[group] for group in STEM_ORDER])
        ax.set_xticklabels([STEM_LABELS[group] for group in STEM_ORDER])
        ax.set_ylabel(TMRM_FEATURE_YLABELS[feature])
        if feature == "tmrm_high_fraction":
            ax.set_ylim(bottom=0)
        add_panel_label(ax, panel_label)

    axes_flat[-1].axis("off")
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    save_figure(fig, output_dir / "tmrm_features_stem_like_percentile_p75")
    plt.close(fig)


def plot_nadph_delta(
    delta: pd.DataFrame,
    output_dir: Path,
    group_order: list[str],
    group_labels: dict[str, str],
    output_name: str,
    figsize: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    axes_flat = axes.ravel()
    x_base = {group: idx for idx, group in enumerate(group_order)}
    offsets = {"control": -0.22, "rotenone": 0.0, "oligomycin": 0.22}
    for ax, metric, panel_label in zip(axes_flat, METRIC_ORDER, ["a", "b", "c", "d"]):
        metric_df = delta[delta["metric"].astype(str) == metric]
        ax.axhline(0, color="#A6A6A6", linewidth=0.7, linestyle="--", zorder=0)
        for condition in CONDITION_ORDER:
            condition_df = metric_df[metric_df["condition"].astype(str) == condition]
            for group in group_order:
                values = condition_df[condition_df["gfp_group"].astype(str) == group][
                    "delta_high_minus_low"
                ].to_numpy(dtype=float)
                x = x_base[group] + offsets[condition]
                jitter = np.linspace(-0.035, 0.035, max(len(values), 1))[: len(values)]
                ax.scatter(
                    np.full(len(values), x) + jitter,
                    values,
                    s=16,
                    color=CONDITION_COLORS[condition],
                    edgecolors="black",
                    linewidths=0.3,
                    alpha=0.9,
                    zorder=3,
                    label=CONDITION_LABELS[condition] if group == group_order[0] else None,
                )
                draw_mean_ci(ax, x, values, CONDITION_COLORS[condition])
        ax.set_title(METRIC_LABELS[metric], fontsize=7.7, fontweight="bold")
        ax.set_xticks([x_base[group] for group in group_order])
        ax.set_xticklabels([group_labels[group] for group in group_order])
        ax.set_ylabel(DELTA_LABELS[metric])
        add_panel_label(ax, panel_label)
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    save_figure(fig, output_dir / output_name)
    plt.close(fig)


def plot_nadph_high_low_paired(paired: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.35), constrained_layout=True)
    axes_flat = axes.ravel()
    x_base = {group: idx for idx, group in enumerate(STEM_ORDER)}
    condition_offsets = {"control": -0.23, "rotenone": 0.0, "oligomycin": 0.23}
    region_offsets = {"low_tmrm": -0.045, "high_tmrm": 0.045}
    region_order = ["low_tmrm", "high_tmrm"]

    for ax, metric, panel_label in zip(axes_flat, METRIC_ORDER, ["a", "b", "c", "d"]):
        metric_df = paired[paired["metric"].astype(str) == metric]
        for condition in CONDITION_ORDER:
            condition_df = metric_df[metric_df["condition"].astype(str) == condition]
            color = CONDITION_COLORS[condition]
            for gfp_group in STEM_ORDER:
                group_df = condition_df[condition_df["gfp_group"].astype(str) == gfp_group]
                center = x_base[gfp_group] + condition_offsets[condition]
                for _, sample_df in group_df.groupby("sample", observed=False):
                    points = []
                    for region_name in region_order:
                        row = sample_df[sample_df["tmrm_region"].astype(str) == region_name]
                        if not row.empty:
                            points.append(
                                (
                                    center + region_offsets[region_name],
                                    float(row.iloc[0]["mean"]),
                                )
                            )
                    if len(points) == 2:
                        ax.plot(
                            [points[0][0], points[1][0]],
                            [points[0][1], points[1][1]],
                            color=color,
                            linewidth=0.55,
                            alpha=0.42,
                            zorder=1,
                        )
                for region_name in region_order:
                    values = group_df.loc[
                        group_df["tmrm_region"].astype(str) == region_name, "mean"
                    ].to_numpy(dtype=float)
                    x = center + region_offsets[region_name]
                    face = "white" if region_name == "low_tmrm" else color
                    ax.scatter(
                        np.full(values.size, x),
                        values,
                        s=16,
                        facecolor=face,
                        edgecolor=color,
                        linewidth=0.7,
                        zorder=3,
                    )
                    mean, ci = mean_ci(values)
                    if np.isfinite(mean):
                        ax.errorbar(
                            x,
                            mean,
                            yerr=ci,
                            fmt="o",
                            color="black",
                            markerfacecolor=face,
                            markeredgecolor="black",
                            markersize=3.5,
                            elinewidth=0.7,
                            capsize=2,
                            zorder=4,
                        )
        ax.set_title(METRIC_LABELS[metric], fontsize=7.7, fontweight="bold")
        ax.set_xticks([x_base[group] for group in STEM_ORDER])
        ax.set_xticklabels([STEM_LABELS[group] for group in STEM_ORDER])
        ax.set_ylabel(NADPH_REGION_YLABELS[metric])
        ax.margins(x=0.08, y=0.14)
        add_panel_label(ax, panel_label)

    condition_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            color=color,
            markerfacecolor=color,
            markeredgecolor="black",
            lw=0,
            markersize=4,
            label=CONDITION_LABELS[condition],
        )
        for condition, color in CONDITION_COLORS.items()
    ]
    region_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            color="#333333",
            markerfacecolor="white",
            markeredgecolor="#333333",
            lw=0,
            markersize=4,
            label="Low-TMRM",
        ),
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            color="#333333",
            markerfacecolor="#333333",
            markeredgecolor="#333333",
            lw=0,
            markersize=4,
            label="TMRM-high",
        ),
    ]
    fig.legend(
        handles=condition_handles + region_handles,
        loc="upper center",
        ncol=5,
        bbox_to_anchor=(0.5, 1.03),
    )
    save_figure(fig, output_dir / "real_nadph_high_low_paired_stem_like_percentile_p75")
    plt.close(fig)


def plot_lipid_droplet_density(lipid: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(3.35, 2.55), constrained_layout=True)
    x_positions = {condition: idx for idx, condition in enumerate(CONDITION_ORDER)}
    value_col = "lipid_droplet_density_per_10000_valid_pixels"
    for condition in CONDITION_ORDER:
        values = lipid[lipid["condition"].astype(str) == condition][value_col].to_numpy(dtype=float)
        x = x_positions[condition]
        jitter = np.linspace(-0.045, 0.045, max(len(values), 1))[: len(values)]
        ax.scatter(
            np.full(len(values), x) + jitter,
            values,
            s=20,
            color=CONDITION_COLORS[condition],
            edgecolors="black",
            linewidths=0.35,
            alpha=0.92,
            zorder=3,
        )
        draw_mean_ci(ax, x, values, CONDITION_COLORS[condition])
    ax.set_xticks([x_positions[condition] for condition in CONDITION_ORDER])
    ax.set_xticklabels([CONDITION_LABELS[condition] for condition in CONDITION_ORDER], rotation=25, ha="right")
    ax.set_ylabel(f"Lipid-like droplets per {DROPLET_DENSITY_SCALE:,} valid pixels")
    ax.set_title("Oxidized-lipid-like droplet density", fontsize=7.7, fontweight="bold")
    save_figure(fig, output_dir / "oxidized_lipid_droplet_density_by_condition_percentile_p75")
    plt.close(fig)


def plot_lipid_droplet_metrics(
    droplet_metrics: pd.DataFrame,
    sample_means: pd.DataFrame,
    output_dir: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(6.7, 5.25), constrained_layout=True)
    axes_flat = axes.ravel()
    x_positions = {condition: idx for idx, condition in enumerate(CONDITION_ORDER)}
    sample_offsets = np.linspace(-0.055, 0.055, 5)

    for ax, metric, panel_label in zip(axes_flat, METRIC_ORDER, ["a", "b", "c", "d"]):
        droplet_df = droplet_metrics[droplet_metrics["metric"].astype(str) == metric]
        sample_df = sample_means[sample_means["metric"].astype(str) == metric]
        for condition in CONDITION_ORDER:
            x = x_positions[condition]
            condition_droplets = droplet_df[droplet_df["condition"].astype(str) == condition]
            if not condition_droplets.empty:
                ranks = condition_droplets.groupby("sample", observed=False).cumcount().to_numpy(dtype=float)
                jitter = ((ranks % 9) - 4) * 0.012
                ax.scatter(
                    np.full(len(condition_droplets), x) + jitter,
                    condition_droplets["value"].to_numpy(dtype=float),
                    s=8,
                    color=CONDITION_COLORS[condition],
                    edgecolors="none",
                    alpha=0.28,
                    zorder=2,
                )
            condition_samples = sample_df[sample_df["condition"].astype(str) == condition]
            values = condition_samples["value"].to_numpy(dtype=float)
            offsets = sample_offsets[: len(values)]
            if len(offsets) < len(values):
                offsets = np.linspace(-0.055, 0.055, len(values))
            ax.scatter(
                np.full(len(values), x) + offsets,
                values,
                s=22,
                color=CONDITION_COLORS[condition],
                edgecolors="black",
                linewidths=0.35,
                alpha=0.95,
                zorder=4,
                label=CONDITION_LABELS[condition] if metric == METRIC_ORDER[0] else None,
            )
            draw_mean_ci(ax, x, values, CONDITION_COLORS[condition])
        ax.set_title(METRIC_LABELS[metric], fontsize=7.7, fontweight="bold")
        ax.set_xticks([x_positions[condition] for condition in CONDITION_ORDER])
        ax.set_xticklabels([CONDITION_LABELS[condition] for condition in CONDITION_ORDER], rotation=25, ha="right")
        ax.set_ylabel(DROPLET_YLABELS[metric])
        add_panel_label(ax, panel_label)

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    save_figure(fig, output_dir / "oxidized_lipid_droplet_metrics_by_condition_percentile_p75")
    plt.close(fig)


def validate_outputs(gfp_qc: pd.DataFrame, tmrm_stem: pd.DataFrame, lipid: pd.DataFrame, lipid_objects: pd.DataFrame) -> None:
    gfp_columns = {"low_pixels", "medium_pixels", "cleaned_high_pixels", "valid_common_pixels"}
    if gfp_columns.issubset(gfp_qc.columns):
        total = gfp_qc["low_pixels"] + gfp_qc["medium_pixels"] + gfp_qc["cleaned_high_pixels"]
        if not np.all(total.to_numpy(dtype=float) == gfp_qc["valid_common_pixels"].to_numpy(dtype=float)):
            raise ValueError("GFP groups do not sum to valid common pixels")
    sample_valid = gfp_qc.set_index(["condition", "sample"])["valid_common_pixels"]
    plot_total = tmrm_stem.groupby(["condition", "sample"], observed=False)["total_pixels"].sum()
    for key, value in plot_total.items():
        if key in sample_valid.index and int(value) != int(sample_valid.loc[key]):
            raise ValueError(f"Stem-like plus low-GFP pixels do not sum to valid pixels for {key}")
    if len(lipid_objects) and lipid_objects["area_pixels"].min() < lipid["lipid_min_object_size"].min():
        raise ValueError("Lipid droplet object table contains objects below the configured size filter")


def main() -> None:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else analysis_dir / "figures_lgr5_stem_like_percentile_p75"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    apply_style()

    region, gfp_qc, tmrm_gfp, lipid_samples, lipid_objects = load_tables(analysis_dir)
    tmrm_stem = build_stem_tmrm_allocation(region)
    tmrm_features = build_tmrm_features(tmrm_stem, tmrm_gfp)
    delta_stem = build_stem_nadph_delta(region)
    delta_three = build_three_group_nadph_delta(region)
    paired_stem = build_nadph_high_low_long(delta_stem)
    lipid_density = build_lipid_droplet_density(lipid_samples)
    lipid_metric_objects, lipid_metric_sample_means = build_lipid_droplet_metrics(lipid_objects)
    validate_outputs(gfp_qc, tmrm_stem, lipid_density, lipid_objects)

    tmrm_stem.to_csv(output_dir / "source_data_tmrm_high_fraction_stem_like_percentile_p75.csv", index=False)
    tmrm_features.to_csv(output_dir / "source_data_tmrm_features_stem_like_percentile_p75.csv", index=False)
    delta_stem.to_csv(output_dir / "source_data_real_nadph_delta_stem_like_percentile_p75.csv", index=False)
    delta_three.to_csv(output_dir / "source_data_real_nadph_delta_three_gfp_groups_percentile_p75.csv", index=False)
    paired_stem.to_csv(
        output_dir / "source_data_real_nadph_high_low_paired_stem_like_percentile_p75.csv",
        index=False,
    )
    lipid_density.to_csv(output_dir / "source_data_oxidized_lipid_droplet_density_by_condition_percentile_p75.csv", index=False)
    lipid_objects.to_csv(output_dir / "source_data_oxidized_lipid_droplet_objects_percentile_p75.csv", index=False)
    lipid_metric_objects.to_csv(
        output_dir / "source_data_oxidized_lipid_droplet_metrics_by_condition_percentile_p75.csv",
        index=False,
    )
    lipid_metric_sample_means.to_csv(
        output_dir / "source_data_oxidized_lipid_droplet_metric_sample_means_percentile_p75.csv",
        index=False,
    )
    gfp_qc.to_csv(output_dir / "source_data_lgr5_gfp_group_qc_percentile_p75.csv", index=False)

    summarize_groups(tmrm_stem, ["condition", "gfp_group"], "tmrm_high_fraction").to_csv(
        output_dir / "summary_tmrm_high_fraction_stem_like_percentile_p75.csv",
        index=False,
    )
    summarize_groups(tmrm_features, ["condition", "gfp_group", "feature"], "value").to_csv(
        output_dir / "summary_tmrm_features_stem_like_percentile_p75.csv",
        index=False,
    )
    summarize_groups(delta_stem, ["condition", "gfp_group", "metric"], "delta_high_minus_low").to_csv(
        output_dir / "summary_real_nadph_delta_stem_like_percentile_p75.csv",
        index=False,
    )
    summarize_groups(paired_stem, ["condition", "gfp_group", "metric", "tmrm_region"], "mean").to_csv(
        output_dir / "summary_real_nadph_high_low_paired_stem_like_percentile_p75.csv",
        index=False,
    )
    summarize_groups(delta_three, ["condition", "gfp_group", "metric"], "delta_high_minus_low").to_csv(
        output_dir / "summary_real_nadph_delta_three_gfp_groups_percentile_p75.csv",
        index=False,
    )
    summarize_groups(
        lipid_density,
        ["condition"],
        "lipid_droplet_density_per_10000_valid_pixels",
    ).to_csv(
        output_dir / "summary_oxidized_lipid_droplet_density_by_condition_percentile_p75.csv",
        index=False,
    )
    summarize_groups(lipid_objects, ["condition"], "area_pixels").to_csv(
        output_dir / "summary_oxidized_lipid_droplet_area_by_condition_percentile_p75.csv",
        index=False,
    )
    summarize_groups(lipid_metric_sample_means, ["condition", "metric"], "value").to_csv(
        output_dir / "summary_oxidized_lipid_droplet_metrics_by_condition_percentile_p75.csv",
        index=False,
    )

    pd.concat(
        [
            tmrm_stem_stats(tmrm_stem),
        ],
        ignore_index=True,
    ).to_csv(output_dir / "stats_tmrm_high_fraction_stem_like_percentile_p75.csv", index=False)
    tmrm_feature_stats(tmrm_features).to_csv(
        output_dir / "stats_tmrm_features_stem_like_percentile_p75.csv",
        index=False,
    )
    pd.concat(
        [
            nadph_delta_stats(delta_stem),
            nadph_stem_group_stats(delta_stem),
        ],
        ignore_index=True,
    ).to_csv(output_dir / "stats_real_nadph_delta_stem_like_percentile_p75.csv", index=False)
    nadph_delta_stats(delta_stem).to_csv(
        output_dir / "stats_real_nadph_high_low_paired_stem_like_percentile_p75.csv",
        index=False,
    )
    nadph_delta_stats(delta_three).to_csv(
        output_dir / "stats_real_nadph_delta_three_gfp_groups_percentile_p75.csv",
        index=False,
    )
    lipid_density_stats(lipid_density).to_csv(
        output_dir / "stats_oxidized_lipid_droplet_density_by_condition_percentile_p75.csv",
        index=False,
    )
    lipid_droplet_metric_stats(lipid_metric_sample_means).to_csv(
        output_dir / "stats_oxidized_lipid_droplet_metrics_by_condition_percentile_p75.csv",
        index=False,
    )

    plot_tmrm_allocation_stem(tmrm_stem, output_dir)
    plot_tmrm_features(tmrm_features, output_dir)
    plot_nadph_delta(
        delta_stem,
        output_dir,
        STEM_ORDER,
        STEM_LABELS,
        "real_nadph_delta_stem_like_percentile_p75",
        (6.7, 5.25),
    )
    plot_nadph_high_low_paired(paired_stem, output_dir)
    plot_nadph_delta(
        delta_three,
        output_dir,
        GFP_ORDER,
        GFP_LABELS,
        "support_real_nadph_delta_three_gfp_groups_percentile_p75",
        (7.3, 5.25),
    )
    plot_lipid_droplet_density(lipid_density, output_dir)
    plot_lipid_droplet_metrics(lipid_metric_objects, lipid_metric_sample_means, output_dir)

    expected_samples = tmrm_stem[["condition", "sample"]].drop_duplicates().shape[0]
    print(f"Wrote stem-like figures for {expected_samples} samples to: {output_dir}")


if __name__ == "__main__":
    main()
