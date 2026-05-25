#!/usr/bin/env python
"""Reusable percentile-p75 TMRM/NAD(P)H plotting workflow."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


METHOD = "percentile_p75"
ALPHA = 0.05
X_SPACING = 0.72

CONDITION_ORDER = ["control", "rotenone", "AA", "oligomycin"]
CONDITION_LABELS = {
    "control": "Control",
    "oligomycin": "Oligomycin",
    "rotenone": "Rotenone",
    "AA": "Antimycin A",
}
CONDITION_COLORS = {
    "control": "#6F7378",
    "oligomycin": "#3B78B8",
    "rotenone": "#D8893A",
    "AA": "#C95A71",
}
METRIC_ORDER = ["Int", "TauPhase", "TauModulation", "TauNormalized"]
METRIC_LABELS = {
    "Int": "NADPH intensity",
    "TauPhase": "Tau phase",
    "TauModulation": "Tau modulation",
    "TauNormalized": "Tau normalized",
}
DELTA_YLABELS = {
    "Int": "Δ NAD(P)H photon/pixel\n(high-TMRM region - low-TMRM region)",
    "TauPhase": "Δ NAD(P)H τ(φ) ns\n(high-TMRM region - low-TMRM region)",
    "TauModulation": "Δ NAD(P)H τ(m) ns\n(high-TMRM region - low-TMRM region)",
    "TauNormalized": "Δ NAD(P)H τ(N) ns\n(high-TMRM region - low-TMRM region)",
}
REGION_YLABELS = {
    "Int": "NAD(P)H photon/pixel",
    "TauPhase": "NAD(P)H τ(φ) ns",
    "TauModulation": "NAD(P)H τ(m) ns",
    "TauNormalized": "NAD(P)H τ(N) ns",
}


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
        "legend.frameon": False,
    }
)


def p_to_star(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "NA"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def p_to_text(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "p = NA"
    if p_value < 0.001:
        return "p < 0.001"
    return f"p = {p_value:.3f}"


def mean_ci(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return math.nan, math.nan
    mean = float(np.mean(values))
    if values.size < 2:
        return mean, math.nan
    ci = float(stats.t.ppf(0.975, values.size - 1) * stats.sem(values))
    return mean, ci


def normality(values: np.ndarray) -> dict[str, object]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3 or np.ptp(values) == 0:
        return {"n": int(values.size), "statistic": math.nan, "p_value": math.nan, "normality_pass": False}
    result = stats.shapiro(values)
    return {
        "n": int(values.size),
        "statistic": float(result.statistic),
        "p_value": float(result.pvalue),
        "normality_pass": bool(result.pvalue > ALPHA),
    }


def holm_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return adjusted.tolist()
    finite_indices = np.where(finite)[0]
    order = finite_indices[np.argsort(p[finite])]
    m = len(order)
    running_max = 0.0
    for rank, idx in enumerate(order):
        running_max = max(running_max, p[idx] * (m - rank))
        adjusted[idx] = min(running_max, 1.0)
    return adjusted.tolist()


def welch_anova(groups: dict[str, np.ndarray]) -> dict[str, float]:
    arrays = [np.asarray(groups[label], dtype=float) for label in groups]
    arrays = [arr[np.isfinite(arr)] for arr in arrays]
    n = np.asarray([arr.size for arr in arrays], dtype=float)
    means = np.asarray([np.mean(arr) for arr in arrays], dtype=float)
    variances = np.asarray([np.var(arr, ddof=1) for arr in arrays], dtype=float)
    if np.any(n < 2) or np.any(variances <= 0):
        return {"statistic": math.nan, "df1": math.nan, "df2": math.nan, "p_value": math.nan}
    weights = n / variances
    weight_sum = float(weights.sum())
    weighted_mean = float(np.sum(weights * means) / weight_sum)
    k = len(arrays)
    df1 = k - 1
    statistic = float(np.sum(weights * (means - weighted_mean) ** 2) / df1)
    w_rel = weights / weight_sum
    tmp = float(np.sum(((1 - w_rel) ** 2) / (n - 1)) / (k**2 - 1))
    if tmp <= 0:
        return {"statistic": statistic, "df1": float(df1), "df2": math.nan, "p_value": math.nan}
    statistic = statistic / (1 + 2 * (k - 2) * tmp)
    df2 = 1 / (3 * tmp)
    p_value = float(stats.f.sf(statistic, df1, df2))
    return {"statistic": statistic, "df1": float(df1), "df2": float(df2), "p_value": p_value}


def dunnett_t3_vs_control(groups: dict[str, np.ndarray], metric: str) -> list[dict[str, object]]:
    control_values = np.asarray(groups["control"], dtype=float)
    control_values = control_values[np.isfinite(control_values)]
    rows: list[dict[str, object]] = []
    k = len(groups)
    for condition in groups:
        if condition == "control":
            continue
        treatment_values = np.asarray(groups[condition], dtype=float)
        treatment_values = treatment_values[np.isfinite(treatment_values)]
        n1 = control_values.size
        n2 = treatment_values.size
        v1 = float(np.var(control_values, ddof=1))
        v2 = float(np.var(treatment_values, ddof=1))
        mean_diff = float(np.mean(treatment_values) - np.mean(control_values))
        se2 = v1 / n1 + v2 / n2
        if n1 < 2 or n2 < 2 or se2 <= 0:
            statistic = df = p_value = math.nan
        else:
            df = float(se2**2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)))
            statistic = float(abs(mean_diff) / math.sqrt(se2))
            p_value = float(stats.studentized_range.sf(math.sqrt(2) * statistic, k, df))
        rows.append(
            {
                "analysis_level": "posthoc_vs_control",
                "metric": metric,
                "condition": condition,
                "comparison": f"{CONDITION_LABELS[condition]} vs Control",
                "test": "Dunnett T3",
                "p_adjustment": "Dunnett T3 familywise",
                "n": int(n2),
                "statistic": statistic,
                "df1": math.nan,
                "df2": df,
                "p_value": p_value,
                "p_adjusted": p_value,
                "effect": mean_diff,
                "normality_pass": True,
            }
        )
    return rows


def dunn_vs_control(groups: dict[str, np.ndarray], metric: str) -> list[dict[str, object]]:
    records = []
    for condition, values in groups.items():
        for value in values:
            if np.isfinite(value):
                records.append({"condition": condition, "value": float(value)})
    rank_df = pd.DataFrame(records)
    rank_df["rank"] = stats.rankdata(rank_df["value"].to_numpy(dtype=float), method="average")
    n_total = len(rank_df)
    tie_counts = rank_df["value"].value_counts().to_numpy(dtype=float)
    tie_correction = 1 - np.sum(tie_counts**3 - tie_counts) / (n_total**3 - n_total)
    if tie_correction <= 0:
        tie_correction = math.nan
    rank_means = rank_df.groupby("condition")["rank"].mean()
    n_by_group = rank_df.groupby("condition")["rank"].size()
    raw_rows: list[dict[str, object]] = []
    raw_p: list[float] = []
    for condition in groups:
        if condition == "control":
            continue
        denom = math.sqrt(
            (n_total * (n_total + 1) / 12)
            * tie_correction
            * (1 / n_by_group["control"] + 1 / n_by_group[condition])
        )
        effect = float(rank_means[condition] - rank_means["control"])
        if denom > 0:
            z_value = effect / denom
            p_value = float(2 * stats.norm.sf(abs(z_value)))
        else:
            z_value = p_value = math.nan
        raw_p.append(p_value)
        raw_rows.append(
            {
                "analysis_level": "posthoc_vs_control",
                "metric": metric,
                "condition": condition,
                "comparison": f"{CONDITION_LABELS[condition]} vs Control",
                "test": "Dunn",
                "p_adjustment": "Holm",
                "n": int(n_by_group[condition]),
                "statistic": z_value,
                "df1": math.nan,
                "df2": math.nan,
                "p_value": p_value,
                "effect": effect,
                "normality_pass": False,
            }
        )
    for row, p_adj in zip(raw_rows, holm_adjust(raw_p)):
        row["p_adjusted"] = p_adj
    return raw_rows


def compute_statistics(delta_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for metric in METRIC_ORDER:
        metric_df = delta_df[delta_df["metric"] == metric]
        groups = {
            condition: metric_df.loc[
                metric_df["condition"] == condition, "mean_diff_high_minus_low"
            ].to_numpy(dtype=float)
            for condition in CONDITION_ORDER
        }
        normality_passes = {}
        for condition, values in groups.items():
            norm = normality(values)
            normality_passes[condition] = bool(norm["normality_pass"])
            rows.append(
                {
                    "analysis_level": "normality",
                    "metric": metric,
                    "condition": condition,
                    "comparison": "",
                    "test": "Shapiro-Wilk",
                    "p_adjustment": "",
                    "n": norm["n"],
                    "statistic": norm["statistic"],
                    "df1": math.nan,
                    "df2": math.nan,
                    "p_value": norm["p_value"],
                    "p_adjusted": math.nan,
                    "effect": math.nan,
                    "normality_pass": norm["normality_pass"],
                }
            )

        for condition, values in groups.items():
            values = values[np.isfinite(values)]
            if normality_passes[condition]:
                result = stats.ttest_1samp(values, popmean=0.0)
                test_name = "one-sample two-tailed t-test"
                statistic = float(result.statistic)
                p_value = float(result.pvalue)
            else:
                result = stats.wilcoxon(values, zero_method="wilcox", alternative="two-sided")
                test_name = "Wilcoxon signed-rank"
                statistic = float(result.statistic)
                p_value = float(result.pvalue)
            rows.append(
                {
                    "analysis_level": "within_delta_zero",
                    "metric": metric,
                    "condition": condition,
                    "comparison": f"{CONDITION_LABELS[condition]} delta vs 0",
                    "test": test_name,
                    "p_adjustment": "",
                    "n": int(values.size),
                    "statistic": statistic,
                    "df1": math.nan,
                    "df2": math.nan,
                    "p_value": p_value,
                    "p_adjusted": math.nan,
                    "effect": float(np.mean(values)),
                    "normality_pass": normality_passes[condition],
                }
            )

        all_normal = all(normality_passes.values())
        if all_normal:
            bf = stats.levene(*[groups[c] for c in CONDITION_ORDER], center="median")
            rows.append(
                {
                    "analysis_level": "omnibus_variance",
                    "metric": metric,
                    "condition": "all",
                    "comparison": "all groups",
                    "test": "Brown-Forsythe variance test",
                    "p_adjustment": "",
                    "n": int(sum(len(groups[c]) for c in CONDITION_ORDER)),
                    "statistic": float(bf.statistic),
                    "df1": math.nan,
                    "df2": math.nan,
                    "p_value": float(bf.pvalue),
                    "p_adjusted": math.nan,
                    "effect": math.nan,
                    "normality_pass": True,
                }
            )
            welch = welch_anova(groups)
            rows.append(
                {
                    "analysis_level": "omnibus",
                    "metric": metric,
                    "condition": "all",
                    "comparison": "all groups",
                    "test": "Welch ANOVA",
                    "p_adjustment": "",
                    "n": int(sum(len(groups[c]) for c in CONDITION_ORDER)),
                    "statistic": welch["statistic"],
                    "df1": welch["df1"],
                    "df2": welch["df2"],
                    "p_value": welch["p_value"],
                    "p_adjusted": math.nan,
                    "effect": math.nan,
                    "normality_pass": True,
                }
            )
            rows.extend(dunnett_t3_vs_control(groups, metric))
        else:
            kw = stats.kruskal(*[groups[c] for c in CONDITION_ORDER])
            rows.append(
                {
                    "analysis_level": "omnibus",
                    "metric": metric,
                    "condition": "all",
                    "comparison": "all groups",
                    "test": "Kruskal-Wallis",
                    "p_adjustment": "",
                    "n": int(sum(len(groups[c]) for c in CONDITION_ORDER)),
                    "statistic": float(kw.statistic),
                    "df1": len(CONDITION_ORDER) - 1,
                    "df2": math.nan,
                    "p_value": float(kw.pvalue),
                    "p_adjusted": math.nan,
                    "effect": math.nan,
                    "normality_pass": False,
                }
            )
            rows.extend(dunn_vs_control(groups, metric))
    return pd.DataFrame(rows)


def resolve_analysis_dir(path: Path) -> Path:
    path = path.resolve()
    sample_file = path / "combined_nadph_sample_differences_all_methods.csv"
    region_file = path / "combined_nadph_region_summary_all_methods.csv"
    if sample_file.exists() and region_file.exists():
        return path
    nested = path / "analysis_TMRM_NADPH_all_conditions"
    if (
        (nested / "combined_nadph_sample_differences_all_methods.csv").exists()
        and (nested / "combined_nadph_region_summary_all_methods.csv").exists()
    ):
        return nested
    raise FileNotFoundError(
        "Could not find combined NAD(P)H CSVs in the provided folder or "
        "analysis_TMRM_NADPH_all_conditions subfolder."
    )


def load_source_data(input_dir: Path, method: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    delta = pd.read_csv(input_dir / "combined_nadph_sample_differences_all_methods.csv")
    region = pd.read_csv(input_dir / "combined_nadph_region_summary_all_methods.csv")
    delta = delta[delta["method"] == method].copy()
    region = region[region["method"] == method].copy()
    missing_conditions = sorted(set(CONDITION_ORDER) - set(delta["condition"].unique()))
    missing_metrics = sorted(set(METRIC_ORDER) - set(delta["metric"].unique()))
    if missing_conditions:
        raise ValueError(f"Missing conditions in delta data: {missing_conditions}")
    if missing_metrics:
        raise ValueError(f"Missing metrics in delta data: {missing_metrics}")
    delta["condition_label"] = delta["condition"].map(CONDITION_LABELS)
    region["condition_label"] = region["condition"].map(CONDITION_LABELS)
    source = delta[
        [
            "condition",
            "condition_label",
            "sample",
            "method",
            "metric",
            "high_n_pixels",
            "low_n_pixels",
            "high_mean",
            "low_mean",
            "mean_diff_high_minus_low",
            "high_median",
            "low_median",
            "median_diff_high_minus_low",
        ]
    ].copy()
    return delta, region, source


def setup_axes_grid(figsize: tuple[float, float]) -> tuple[plt.Figure, np.ndarray]:
    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    return fig, axes.ravel()


def condition_positions() -> dict[str, float]:
    return {condition: i * X_SPACING for i, condition in enumerate(CONDITION_ORDER)}


def add_control_brackets(ax: plt.Axes, positions: dict[str, float], posthoc: pd.DataFrame) -> None:
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    if span <= 0:
        return
    base_y = ymax + 0.04 * span
    gap = 0.085 * span
    bracket_h = 0.025 * span
    ax.set_ylim(ymin, base_y + gap * 3 + 0.08 * span)
    control_x = positions["control"]
    for idx, condition in enumerate(CONDITION_ORDER[1:]):
        row = posthoc[posthoc["condition"] == condition]
        if row.empty:
            continue
        x1 = positions[condition]
        y = base_y + gap * idx
        ax.plot(
            [control_x, control_x, x1, x1],
            [y, y + bracket_h, y + bracket_h, y],
            color="#202020",
            lw=0.65,
            clip_on=False,
        )
        ax.text(
            (control_x + x1) / 2,
            y + bracket_h + 0.012 * span,
            p_to_star(float(row.iloc[0]["p_adjusted"])),
            ha="center",
            va="bottom",
            fontsize=7,
            color="#202020",
        )


def plot_delta(delta: pd.DataFrame, stats_df: pd.DataFrame, output_base: Path) -> None:
    fig, axes = setup_axes_grid((5.8, 5.2))
    rng = np.random.default_rng(20260520)
    positions = condition_positions()
    for ax, metric in zip(axes, METRIC_ORDER):
        metric_df = delta[delta["metric"] == metric]
        for condition in CONDITION_ORDER:
            values = metric_df.loc[
                metric_df["condition"] == condition, "mean_diff_high_minus_low"
            ].to_numpy(dtype=float)
            x = positions[condition]
            color = CONDITION_COLORS[condition]
            ax.scatter(
                np.full(values.size, x) + rng.normal(0, 0.045, size=values.size),
                values,
                s=18,
                color=color,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
            mean, ci = mean_ci(values)
            ax.errorbar(
                x,
                mean,
                yerr=ci,
                fmt="o",
                color="black",
                markerfacecolor="white",
                markeredgecolor="black",
                markersize=4,
                elinewidth=0.8,
                capsize=2.5,
                zorder=4,
            )
        ax.axhline(0, color="#A6A6A6", lw=0.7, ls="--", zorder=1)
        ax.set_xticks([positions[c] for c in CONDITION_ORDER])
        ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER], rotation=30, ha="right")
        ax.set_ylabel(DELTA_YLABELS[metric])
        ax.margins(x=0.08, y=0.22)
        ax.set_xlim(
            positions[CONDITION_ORDER[0]] - 0.38 * X_SPACING,
            positions[CONDITION_ORDER[-1]] + 0.38 * X_SPACING,
        )
        posthoc = stats_df[
            (stats_df["analysis_level"] == "posthoc_vs_control") & (stats_df["metric"] == metric)
        ]
        add_control_brackets(ax, positions, posthoc)
    save_figure(fig, output_base)
    plt.close(fig)


def plot_high_low(region: pd.DataFrame, output_base: Path) -> None:
    fig, axes = setup_axes_grid((5.8, 5.4))
    region_order = ["low_tmrm_outside_mask", "high_tmrm"]
    offsets = {"low_tmrm_outside_mask": -0.16 * X_SPACING, "high_tmrm": 0.16 * X_SPACING}
    positions = condition_positions()
    for ax, metric in zip(axes, METRIC_ORDER):
        metric_df = region[region["metric"] == metric]
        for condition in CONDITION_ORDER:
            condition_df = metric_df[metric_df["condition"] == condition]
            x = positions[condition]
            color = CONDITION_COLORS[condition]
            for _, sample_df in condition_df.groupby("sample"):
                points = []
                for region_name in region_order:
                    row = sample_df[sample_df["region"] == region_name]
                    if not row.empty:
                        points.append((x + offsets[region_name], float(row.iloc[0]["mean"])))
                if len(points) == 2:
                    ax.plot(
                        [points[0][0], points[1][0]],
                        [points[0][1], points[1][1]],
                        color=color,
                        lw=0.55,
                        alpha=0.42,
                        zorder=1,
                    )
            for region_name in region_order:
                values = condition_df.loc[condition_df["region"] == region_name, "mean"].to_numpy(dtype=float)
                xs = np.full(values.size, x + offsets[region_name])
                face = "white" if region_name == "low_tmrm_outside_mask" else color
                ax.scatter(
                    xs,
                    values,
                    s=16,
                    facecolor=face,
                    edgecolor=color,
                    linewidth=0.7,
                    zorder=3,
                )
                mean, ci = mean_ci(values)
                ax.errorbar(
                    x + offsets[region_name],
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
        ax.set_xticks([positions[c] for c in CONDITION_ORDER])
        ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER], rotation=30, ha="right")
        ax.set_ylabel(REGION_YLABELS[metric])
        ax.margins(x=0.08, y=0.14)
        ax.set_xlim(
            positions[CONDITION_ORDER[0]] - 0.42 * X_SPACING,
            positions[CONDITION_ORDER[-1]] + 0.42 * X_SPACING,
        )
    handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            color="#333333",
            markerfacecolor="white",
            markeredgecolor="#333333",
            lw=0,
            markersize=4,
            label="Low-TMRM intracellular",
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
            label="TMRM-high mitochondrial",
        ),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
    save_figure(fig, output_base)
    plt.close(fig)


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


def parse_key_value(items: list[str], target: dict[str, str]) -> dict[str, str]:
    parsed = dict(target)
    for item in items:
        if "=" not in item:
            raise ValueError(f"Expected KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=Path.cwd(), help="Analysis folder or project root.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output folder. Defaults to <analysis-dir>/figures_percentile_p75.")
    parser.add_argument("--method", default=METHOD, help="Mask method to plot.")
    parser.add_argument("--condition-order", nargs="+", default=CONDITION_ORDER, help="Condition keys in plotting order.")
    parser.add_argument("--condition-label", action="append", default=[], help="Override label as key=Label.")
    parser.add_argument("--condition-color", action="append", default=[], help="Override color as key=#RRGGBB.")
    return parser.parse_args()


def configure_from_args(args: argparse.Namespace) -> tuple[Path, Path, str]:
    global CONDITION_ORDER, CONDITION_LABELS, CONDITION_COLORS, METHOD
    input_dir = resolve_analysis_dir(args.analysis_dir)
    output_dir = args.output_dir if args.output_dir else input_dir / "figures_percentile_p75"
    output_dir = output_dir.resolve()
    METHOD = args.method
    CONDITION_ORDER = list(args.condition_order)
    CONDITION_LABELS = parse_key_value(args.condition_label, CONDITION_LABELS)
    CONDITION_COLORS = parse_key_value(args.condition_color, CONDITION_COLORS)
    missing_label = [condition for condition in CONDITION_ORDER if condition not in CONDITION_LABELS]
    missing_color = [condition for condition in CONDITION_ORDER if condition not in CONDITION_COLORS]
    if missing_label:
        raise ValueError(f"Missing condition labels for: {missing_label}")
    if missing_color:
        raise ValueError(f"Missing condition colors for: {missing_color}")
    if "control" not in CONDITION_ORDER:
        raise ValueError("The statistics and brackets require a condition key named 'control'.")
    if CONDITION_ORDER[0] != "control":
        raise ValueError("Put 'control' first in --condition-order for treatment-vs-control annotations.")
    return input_dir, output_dir, METHOD


def main() -> None:
    args = parse_args()
    input_dir, output_dir, method = configure_from_args(args)
    output_dir.mkdir(parents=True, exist_ok=True)
    delta, region, source = load_source_data(input_dir, method)
    source.to_csv(output_dir / "source_data_nadph_percentile_p75.csv", index=False)
    stats_df = compute_statistics(delta)
    stats_df.to_csv(output_dir / "stats_nadph_percentile_p75.csv", index=False)
    plot_delta(delta, stats_df, output_dir / "nadph_delta_percentile_p75")
    plot_high_low(region, output_dir / "nadph_high_low_paired_percentile_p75")
    print(f"Wrote outputs to {output_dir}")


if __name__ == "__main__":
    main()
