#!/usr/bin/env python
"""Plot NADPH phasor G/S shifts from control to treatments.

Inputs may be either:
1. a project root containing analysis_TMRM_NADPH_aligned_inputs/<condition>/TMRM,NADPH, or
2. a project root containing raw <condition>/TMRM and <condition>/NADPH folders.

The plotted unit is the sample/image centroid. Pixels are only used to compute
per-sample G/S centroids inside percentile-p75 TMRM-high and low-TMRM regions.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


METHOD = "percentile_p75"
PAIR_RE = re.compile(r"^.+?_(?P<sample>\d+)_(?P<acq>\d+)_pxData\.csv$")
OUTPUT_STEM = "nadph_phasor_trajectory_percentile_p75"

KNOWN_CONDITIONS = {
    "control": ("control", "Control", "#6F7378", 0),
    "rotenone": ("rotenone", "Rotenone", "#D8893A", 1),
    "antimycin a": ("AA", "Antimycin A", "#C95A71", 2),
    "antimycin_a": ("AA", "Antimycin A", "#C95A71", 2),
    "aa": ("AA", "Antimycin A", "#C95A71", 2),
    "oligomycin": ("oligomycin", "Oligomycin", "#3B78B8", 3),
}
FALLBACK_COLORS = ["#6F7378", "#D8893A", "#C95A71", "#3B78B8", "#5D8A66", "#8E6BB8"]
REGIONS = {
    "high_tmrm_percentile_p75": "TMRM-high mitochondrial region",
    "low_tmrm_percentile_p75": "Low-TMRM intracellular region",
}


@dataclass(frozen=True)
class Condition:
    folder: str
    key: str
    label: str
    color: str
    order: int


@dataclass(frozen=True)
class Pair:
    sample: int
    tmrm_csv: Path
    nadph_csv: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Project root.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help=(
            "Folder containing condition subfolders. Defaults to "
            "<root>/analysis_TMRM_NADPH_aligned_inputs when it exists, otherwise <root>."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output folder. Defaults to "
            "<root>/analysis_TMRM_NADPH_all_conditions/figures_phasor_percentile_p75."
        ),
    )
    parser.add_argument(
        "--condition-order",
        nargs="+",
        default=None,
        help="Optional condition keys/folder names in plotting order; control must be first.",
    )
    return parser.parse_args()


def condition_from_folder(folder: Path, fallback_index: int) -> Condition:
    normalized = folder.name.strip().lower()
    key, label, color, order = KNOWN_CONDITIONS.get(
        normalized,
        (
            re.sub(r"[^a-z0-9]+", "_", normalized).strip("_"),
            folder.name.strip().title(),
            FALLBACK_COLORS[fallback_index % len(FALLBACK_COLORS)],
            100 + fallback_index,
        ),
    )
    return Condition(folder=folder.name, key=key, label=label, color=color, order=order)


def discover_conditions(input_dir: Path, condition_order: list[str] | None = None) -> list[Condition]:
    folders = [
        path
        for path in input_dir.iterdir()
        if path.is_dir() and (path / "TMRM").is_dir() and (path / "NADPH").is_dir()
    ]
    if not folders:
        raise FileNotFoundError(f"No condition folders with TMRM/ and NADPH/ found under {input_dir}")
    conditions = [condition_from_folder(path, idx) for idx, path in enumerate(folders)]
    if condition_order:
        order_lookup = {item.lower(): idx for idx, item in enumerate(condition_order)}
        keyed = {condition.key.lower(): condition for condition in conditions}
        foldered = {condition.folder.lower(): condition for condition in conditions}
        ordered: list[Condition] = []
        for item in condition_order:
            item_key = item.lower()
            match = keyed.get(item_key) or foldered.get(item_key)
            if not match:
                raise ValueError(f"Condition order item {item!r} was not found under {input_dir}")
            ordered.append(match)
        remaining = [
            condition
            for condition in conditions
            if condition.key.lower() not in order_lookup and condition.folder.lower() not in order_lookup
        ]
        return ordered + sorted(remaining, key=lambda c: (c.order, c.label))
    conditions = sorted(conditions, key=lambda c: (c.order, c.label))
    if conditions[0].key != "control":
        raise ValueError("A control condition is required and must sort first; use --condition-order if needed.")
    return conditions


def index_pxdata(folder: Path) -> dict[int, Path]:
    files: dict[int, Path] = {}
    for path in sorted(folder.glob("*_pxData.csv")):
        match = PAIR_RE.match(path.name)
        if not match:
            continue
        sample = int(match.group("sample"))
        if sample in files:
            raise ValueError(f"Multiple pxData CSV files for sample {sample} in {folder}")
        files[sample] = path
    return files


def discover_pairs(root: Path) -> list[Pair]:
    tmrm_files = index_pxdata(root / "TMRM")
    nadph_files = index_pxdata(root / "NADPH")
    samples = sorted(set(tmrm_files) & set(nadph_files))
    missing_tmrm = sorted(set(nadph_files) - set(tmrm_files))
    missing_nadph = sorted(set(tmrm_files) - set(nadph_files))
    if missing_tmrm or missing_nadph:
        raise ValueError(
            "TMRM/NADPH sample mismatch. "
            f"Missing TMRM for {missing_tmrm}; missing NADPH for {missing_nadph}."
        )
    if not samples:
        raise ValueError(f"No paired pxData CSV files found under {root}")
    return [Pair(sample=s, tmrm_csv=tmrm_files[s], nadph_csv=nadph_files[s]) for s in samples]


def read_pair(pair: Pair) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    tmrm = pd.read_csv(pair.tmrm_csv)
    nadph = pd.read_csv(pair.nadph_csv)
    required_tmrm = {"Int", "Ycoord", "Xcoord"}
    required_nadph = {"G", "S", "Int", "Ycoord", "Xcoord"}
    missing_tmrm = required_tmrm - set(tmrm.columns)
    missing_nadph = required_nadph - set(nadph.columns)
    if missing_tmrm or missing_nadph:
        raise ValueError(
            f"Sample {pair.sample} missing columns: "
            f"TMRM={sorted(missing_tmrm)}, NADPH={sorted(missing_nadph)}"
        )
    if tmrm.duplicated(["Xcoord", "Ycoord"]).any() or nadph.duplicated(["Xcoord", "Ycoord"]).any():
        raise ValueError(f"Sample {pair.sample} contains duplicated Xcoord/Ycoord values")

    tmrm = tmrm.copy()
    nadph = nadph.copy()
    tmrm["_tmrm_order"] = np.arange(len(tmrm))
    shared = (
        tmrm[["Xcoord", "Ycoord", "_tmrm_order"]]
        .merge(nadph[["Xcoord", "Ycoord"]], on=["Xcoord", "Ycoord"], how="inner")
        .sort_values("_tmrm_order")
    )
    aligned_tmrm = shared[["Xcoord", "Ycoord"]].merge(tmrm.drop(columns=["_tmrm_order"]), on=["Xcoord", "Ycoord"], how="left")
    aligned_nadph = shared[["Xcoord", "Ycoord"]].merge(nadph, on=["Xcoord", "Ycoord"], how="left")
    alignment = {
        "tmrm_rows_original": int(len(tmrm)),
        "nadph_rows_original": int(len(nadph)),
        "shared_rows_used": int(len(shared)),
        "tmrm_rows_dropped": int(len(tmrm) - len(shared)),
        "nadph_rows_dropped": int(len(nadph) - len(shared)),
    }
    return aligned_tmrm, aligned_nadph, alignment


def percentile_p75_selectors(tmrm: pd.DataFrame) -> dict[str, np.ndarray]:
    values = tmrm["Int"].to_numpy(dtype=float)
    p75 = float(np.nanpercentile(values, 75))
    high = values >= p75
    return {
        "high_tmrm_percentile_p75": high,
        "low_tmrm_percentile_p75": ~high,
    }


def centroid_row(
    condition: Condition,
    pair: Pair,
    region_key: str,
    nadph: pd.DataFrame,
    selector: np.ndarray,
    alignment: dict[str, int],
) -> dict[str, object]:
    sub = nadph.loc[selector, ["G", "S", "Int"]].copy()
    if sub.empty:
        raise ValueError(f"No pixels for {condition.key} sample {pair.sample} region {region_key}")
    row: dict[str, object] = {
        "condition": condition.key,
        "condition_label": condition.label,
        "condition_folder": condition.folder,
        "sample": pair.sample,
        "tmrm_csv": pair.tmrm_csv.name,
        "nadph_csv": pair.nadph_csv.name,
        "method": METHOD,
        "region": region_key,
        "region_label": REGIONS[region_key],
        "n_pixels": int(len(sub)),
        "G_mean": float(sub["G"].mean()),
        "S_mean": float(sub["S"].mean()),
        "G_median": float(sub["G"].median()),
        "S_median": float(sub["S"].median()),
        "G_sd_pixels": float(sub["G"].std(ddof=1)),
        "S_sd_pixels": float(sub["S"].std(ddof=1)),
        "nadph_int_mean": float(sub["Int"].mean()),
    }
    row.update(alignment)
    return row


def load_source_data(input_dir: Path, conditions: list[Condition]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for condition in conditions:
        condition_root = input_dir / condition.folder
        for pair in discover_pairs(condition_root):
            tmrm, nadph, alignment = read_pair(pair)
            selectors = percentile_p75_selectors(tmrm)
            for region_key, selector in selectors.items():
                rows.append(centroid_row(condition, pair, region_key, nadph, selector, alignment))
    return pd.DataFrame(rows)


def sem(values: pd.Series) -> float:
    values = values.astype(float)
    if len(values) < 2:
        return math.nan
    return float(values.std(ddof=1) / math.sqrt(len(values)))


def summarize_groups(source: pd.DataFrame, conditions: list[Condition]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for region_key, region_label in REGIONS.items():
        for condition in conditions:
            sub = source[(source["region"] == region_key) & (source["condition"] == condition.key)]
            rows.append(
                {
                    "condition": condition.key,
                    "condition_label": condition.label,
                    "region": region_key,
                    "region_label": region_label,
                    "n_samples": int(sub["sample"].nunique()),
                    "G_mean": float(sub["G_mean"].mean()),
                    "S_mean": float(sub["S_mean"].mean()),
                    "G_sem": sem(sub["G_mean"]),
                    "S_sem": sem(sub["S_mean"]),
                }
            )
    return pd.DataFrame(rows)


def summarize_trajectories(group_summary: pd.DataFrame, conditions: list[Condition]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    control_key = conditions[0].key
    if control_key != "control":
        raise ValueError("The first condition must be control for control-to-treatment arrows.")
    for region_key in REGIONS:
        control = group_summary[
            (group_summary["region"] == region_key) & (group_summary["condition"] == control_key)
        ].iloc[0]
        for condition in conditions[1:]:
            treatment = group_summary[
                (group_summary["region"] == region_key) & (group_summary["condition"] == condition.key)
            ].iloc[0]
            d_g = float(treatment["G_mean"] - control["G_mean"])
            d_s = float(treatment["S_mean"] - control["S_mean"])
            rows.append(
                {
                    "region": region_key,
                    "region_label": REGIONS[region_key],
                    "from_condition": control_key,
                    "to_condition": condition.key,
                    "to_condition_label": condition.label,
                    "G_start": float(control["G_mean"]),
                    "S_start": float(control["S_mean"]),
                    "G_end": float(treatment["G_mean"]),
                    "S_end": float(treatment["S_mean"]),
                    "delta_G": d_g,
                    "delta_S": d_s,
                    "vector_length": math.hypot(d_g, d_s),
                }
            )
    return pd.DataFrame(rows)


def configure_matplotlib() -> None:
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


def phasor_limits(source: pd.DataFrame) -> tuple[tuple[float, float], tuple[float, float]]:
    x_min = max(0.0, float(source["G_mean"].min()) - 0.030)
    x_max = min(1.0, float(source["G_mean"].max()) + 0.030)
    y_min = max(0.0, float(source["S_mean"].min()) - 0.020)
    y_max = 0.510
    return (x_min, x_max), (y_min, y_max)


def draw_universal_circle(ax: plt.Axes) -> None:
    g = np.linspace(0.0, 1.0, 800)
    s = np.sqrt(np.maximum(0.0, 0.25 - (g - 0.5) ** 2))
    ax.plot(g, s, color="#C5CAD0", lw=0.8, zorder=0)
    ax.text(0.501, 0.496, "universal circle", ha="right", va="bottom", color="#8C939B", fontsize=6.0)


def plot_panel(
    ax: plt.Axes,
    source: pd.DataFrame,
    group_summary: pd.DataFrame,
    trajectories: pd.DataFrame,
    region_key: str,
    conditions: list[Condition],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    draw_universal_circle(ax)
    region_source = source[source["region"] == region_key]
    region_groups = group_summary[group_summary["region"] == region_key]
    region_vectors = trajectories[trajectories["region"] == region_key]

    for condition in conditions:
        sample_df = region_source[region_source["condition"] == condition.key]
        group = region_groups[region_groups["condition"] == condition.key].iloc[0]
        ax.scatter(
            sample_df["G_mean"],
            sample_df["S_mean"],
            s=12,
            color=condition.color,
            alpha=0.28,
            edgecolor="none",
            zorder=2,
        )
        ax.errorbar(
            group["G_mean"],
            group["S_mean"],
            xerr=group["G_sem"],
            yerr=group["S_sem"],
            fmt="o",
            ms=4.6,
            mfc=condition.color,
            mec="white",
            mew=0.55,
            ecolor=condition.color,
            elinewidth=0.75,
            capsize=2.0,
            alpha=0.98,
            zorder=4,
        )

    for _, vector in region_vectors.iterrows():
        color = next(c.color for c in conditions if c.key == vector["to_condition"])
        ax.annotate(
            "",
            xy=(vector["G_end"], vector["S_end"]),
            xytext=(vector["G_start"], vector["S_start"]),
            arrowprops={
                "arrowstyle": "->",
                "color": color,
                "lw": 1.0,
                "alpha": 0.85,
                "shrinkA": 5,
                "shrinkB": 6,
                "mutation_scale": 8,
            },
            zorder=3,
        )

    ax.set_title(REGIONS[region_key], fontsize=7.4, pad=4)
    ax.set_xlabel("G")
    ax.set_ylabel("S")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#ECEFF3", lw=0.45, zorder=0)
    ax.tick_params(labelsize=6.4, length=2.4)


def make_figure(
    source: pd.DataFrame,
    group_summary: pd.DataFrame,
    trajectories: pd.DataFrame,
    conditions: list[Condition],
) -> plt.Figure:
    configure_matplotlib()
    xlim, ylim = phasor_limits(source)
    fig, axes = plt.subplots(1, 2, figsize=(7.05, 3.25), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.83, bottom=0.25, wspace=0.18)

    for ax, region_key in zip(axes, REGIONS):
        plot_panel(ax, source, group_summary, trajectories, region_key, conditions, xlim, ylim)

    fig.suptitle(
        "NADPH phasor shift trajectories from control",
        x=0.075,
        y=0.94,
        ha="left",
        fontsize=8.2,
        fontweight="bold",
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            label=condition.label,
            markerfacecolor=condition.color,
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=5,
        )
        for condition in conditions
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.055),
        ncol=min(len(conditions), 4),
        columnspacing=1.35,
        handletextpad=0.35,
    )
    fig.text(
        0.075,
        0.13,
        "Points are sample centroids; crosshairs show SEM. Arrows connect each treatment mean to the matching control mean.",
        ha="left",
        va="center",
        fontsize=6.2,
        color="#4B5563",
    )
    return fig


def save_figure(fig: plt.Figure, output_dir: Path) -> None:
    base = output_dir / OUTPUT_STEM
    fig.savefig(f"{base}.svg", bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    fig.savefig(f"{base}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{base}.tiff", dpi=600, bbox_inches="tight")


def validate_outputs(source: pd.DataFrame, group_summary: pd.DataFrame, trajectories: pd.DataFrame, conditions: list[Condition]) -> None:
    n_samples = source[["condition", "sample"]].drop_duplicates().shape[0]
    expected_source_rows = n_samples * len(REGIONS)
    expected_group_rows = len(conditions) * len(REGIONS)
    expected_arrows = (len(conditions) - 1) * len(REGIONS)
    if len(source) != expected_source_rows:
        raise ValueError(f"Expected {expected_source_rows} source rows, found {len(source)}")
    if len(group_summary) != expected_group_rows:
        raise ValueError(f"Expected {expected_group_rows} group rows, found {len(group_summary)}")
    if len(trajectories) != expected_arrows:
        raise ValueError(f"Expected {expected_arrows} trajectory vectors, found {len(trajectories)}")


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    default_aligned = root / "analysis_TMRM_NADPH_aligned_inputs"
    input_dir = args.input_dir.resolve() if args.input_dir else (default_aligned if default_aligned.exists() else root)
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else root / "analysis_TMRM_NADPH_all_conditions" / "figures_phasor_percentile_p75"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    conditions = discover_conditions(input_dir, args.condition_order)
    source = load_source_data(input_dir, conditions)
    group_summary = summarize_groups(source, conditions)
    trajectories = summarize_trajectories(group_summary, conditions)
    validate_outputs(source, group_summary, trajectories, conditions)

    source.to_csv(output_dir / f"source_data_{OUTPUT_STEM}.csv", index=False)
    group_summary.to_csv(output_dir / f"group_summary_{OUTPUT_STEM}.csv", index=False)
    trajectories.to_csv(output_dir / f"trajectory_vectors_{OUTPUT_STEM}.csv", index=False)

    fig = make_figure(source, group_summary, trajectories, conditions)
    save_figure(fig, output_dir)
    plt.close(fig)

    print(f"Input directory: {input_dir}")
    print(f"Wrote outputs to {output_dir}")
    print(f"Sample-region centroid rows: {len(source)}")
    print(f"Group centroids: {len(group_summary)}")
    print(f"Trajectory arrows: {len(trajectories)}")


if __name__ == "__main__":
    main()
