#!/usr/bin/env python
"""Generate reusable style variants for percentile-p75 TMRM/NAD(P)H plots."""

from __future__ import annotations

import argparse
import html
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_nadph_percentile_p75 as base


@dataclass(frozen=True)
class StyleSpec:
    name: str
    label: str
    description: str
    delta_size: tuple[float, float]
    paired_size: tuple[float, float]
    font_size: float
    title_size: float
    axis_lw: float
    tick_lw: float
    dot_size: float
    dot_alpha: float
    ci_lw: float
    capsize: float
    mean_marker_size: float
    show_omnibus: bool
    show_panel_labels: bool
    show_metric_titles: bool
    metric_labels_on_y: bool
    show_control_brackets: bool
    x_spacing: float
    palette: dict[str, str]
    markers: dict[str, str]


NATURE_PALETTE = {
    "control": "#767676",
    "rotenone": "#D89042",
    "AA": "#B64342",
    "oligomycin": "#3775BA",
}

GRAPH_PAD_PALETTE = {
    "control": "#5F6368",
    "rotenone": "#E28E2C",
    "AA": "#C95A71",
    "oligomycin": "#2F75B5",
}

GRAYSCALE_PALETTE = {
    "control": "#202020",
    "rotenone": "#5A5A5A",
    "AA": "#8A8A8A",
    "oligomycin": "#B0B0B0",
}

MARKERS = {
    "control": "o",
    "rotenone": "o",
    "AA": "o",
    "oligomycin": "o",
}

GRAYSCALE_MARKERS = {
    "control": "o",
    "rotenone": "s",
    "AA": "^",
    "oligomycin": "D",
}

STYLES = [
    StyleSpec(
        name="graphpad_prism",
        label="GraphPad Prism-like",
        description=(
            "Larger GraphPad-like panel proportions, strong black axes, larger dots, "
            "mean with 95% CI, and clean treatment-vs-control brackets."
        ),
        delta_size=(5.8, 5.2),
        paired_size=(5.8, 5.4),
        font_size=8.2,
        title_size=9.0,
        axis_lw=1.0,
        tick_lw=0.9,
        dot_size=28,
        dot_alpha=0.92,
        ci_lw=1.15,
        capsize=3.0,
        mean_marker_size=4.8,
        show_omnibus=False,
        show_panel_labels=False,
        show_metric_titles=False,
        metric_labels_on_y=True,
        show_control_brackets=True,
        x_spacing=0.72,
        palette=GRAPH_PAD_PALETTE,
        markers=MARKERS,
    ),
    StyleSpec(
        name="nature_compact",
        label="Nature compact",
        description=(
            "Journal-width compact layout, restrained treatment palette, panel labels, "
            "and exact omnibus p-values."
        ),
        delta_size=(7.2, 4.7),
        paired_size=(7.2, 5.0),
        font_size=7.0,
        title_size=7.7,
        axis_lw=0.75,
        tick_lw=0.65,
        dot_size=18,
        dot_alpha=0.9,
        ci_lw=0.85,
        capsize=2.3,
        mean_marker_size=4.0,
        show_omnibus=True,
        show_panel_labels=True,
        show_metric_titles=True,
        metric_labels_on_y=False,
        show_control_brackets=False,
        x_spacing=1.0,
        palette=NATURE_PALETTE,
        markers=MARKERS,
    ),
    StyleSpec(
        name="nature_minimal",
        label="Nature minimal",
        description=(
            "Cleaner manuscript-style view: raw dots and mean/CI only, panel labels, "
            "stars for treatment-vs-control, no omnibus text inside panels."
        ),
        delta_size=(6.6, 4.7),
        paired_size=(6.6, 5.0),
        font_size=7.0,
        title_size=7.5,
        axis_lw=0.7,
        tick_lw=0.6,
        dot_size=18,
        dot_alpha=0.88,
        ci_lw=0.8,
        capsize=2.2,
        mean_marker_size=3.9,
        show_omnibus=False,
        show_panel_labels=True,
        show_metric_titles=True,
        metric_labels_on_y=False,
        show_control_brackets=False,
        x_spacing=1.0,
        palette=NATURE_PALETTE,
        markers=MARKERS,
    ),
    StyleSpec(
        name="print_safe_gray",
        label="Print-safe grayscale",
        description=(
            "Black-and-white friendly version with marker shapes carrying group identity. "
            "Useful for methods supplements or reviewer PDFs."
        ),
        delta_size=(6.8, 5.0),
        paired_size=(6.8, 5.2),
        font_size=7.2,
        title_size=7.8,
        axis_lw=0.8,
        tick_lw=0.7,
        dot_size=22,
        dot_alpha=0.9,
        ci_lw=0.9,
        capsize=2.5,
        mean_marker_size=4.2,
        show_omnibus=True,
        show_panel_labels=True,
        show_metric_titles=True,
        metric_labels_on_y=False,
        show_control_brackets=False,
        x_spacing=1.0,
        palette=GRAYSCALE_PALETTE,
        markers=GRAYSCALE_MARKERS,
    ),
]


def apply_style(style: StyleSpec) -> None:
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": style.font_size,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": style.axis_lw,
            "xtick.major.width": style.tick_lw,
            "ytick.major.width": style.tick_lw,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="#272727",
    )


def condition_positions(style: StyleSpec) -> dict[str, float]:
    return {
        condition: i * style.x_spacing
        for i, condition in enumerate(base.CONDITION_ORDER)
    }


def add_control_brackets(
    ax: plt.Axes,
    positions: dict[str, float],
    posthoc: pd.DataFrame,
    style: StyleSpec,
) -> None:
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    if span <= 0:
        return
    base_y = ymax + 0.04 * span
    gap = 0.085 * span
    bracket_h = 0.025 * span
    ax.set_ylim(ymin, base_y + gap * 3 + 0.08 * span)
    control_x = positions["control"]
    for idx, condition in enumerate(base.CONDITION_ORDER[1:]):
        row = posthoc[posthoc["condition"] == condition]
        if row.empty:
            continue
        x1 = positions[condition]
        y = base_y + gap * idx
        ax.plot(
            [control_x, control_x, x1, x1],
            [y, y + bracket_h, y + bracket_h, y],
            color="#202020",
            lw=style.ci_lw * 0.8,
            clip_on=False,
        )
        ax.text(
            (control_x + x1) / 2,
            y + bracket_h + 0.012 * span,
            base.p_to_star(float(row.iloc[0]["p_adjusted"])),
            ha="center",
            va="bottom",
            fontsize=style.font_size,
            color="#202020",
        )


def save_gallery_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".png"), dpi=600, bbox_inches="tight")


def plot_delta_variant(
    delta: pd.DataFrame, stats_df: pd.DataFrame, style: StyleSpec, output_base: Path
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=style.delta_size, constrained_layout=True)
    axes = axes.ravel()
    rng = np.random.default_rng(20260520)
    positions = condition_positions(style)
    for panel_label, ax, metric in zip(["a", "b", "c", "d"], axes, base.METRIC_ORDER):
        if style.show_panel_labels:
            add_panel_label(ax, panel_label)
        metric_df = delta[delta["metric"] == metric]
        for condition in base.CONDITION_ORDER:
            values = metric_df.loc[
                metric_df["condition"] == condition, "mean_diff_high_minus_low"
            ].to_numpy(dtype=float)
            x = positions[condition]
            color = style.palette[condition]
            marker = style.markers[condition]
            ax.scatter(
                np.full(values.size, x) + rng.normal(0, 0.045, size=values.size),
                values,
                s=style.dot_size,
                marker=marker,
                color=color,
                alpha=style.dot_alpha,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
            mean, ci = base.mean_ci(values)
            ax.errorbar(
                x,
                mean,
                yerr=ci,
                fmt="o",
                color="#202020",
                markerfacecolor="white",
                markeredgecolor="#202020",
                markersize=style.mean_marker_size,
                elinewidth=style.ci_lw,
                capsize=style.capsize,
                zorder=4,
            )
        ax.axhline(0, color="#A6A6A6", lw=0.7, ls="--", zorder=1)
        if style.show_metric_titles:
            ax.set_title(base.METRIC_LABELS[metric], fontweight="bold", fontsize=style.title_size)
        ax.set_xticks([positions[c] for c in base.CONDITION_ORDER])
        ax.set_xticklabels(
            [base.CONDITION_LABELS[c] for c in base.CONDITION_ORDER],
            rotation=30,
            ha="right",
        )
        ax.set_ylabel(base.DELTA_YLABELS[metric] if style.metric_labels_on_y else "High TMRM - Low TMRM")
        ax.margins(x=0.08, y=0.22)
        ax.set_xlim(
            positions[base.CONDITION_ORDER[0]] - 0.38 * style.x_spacing,
            positions[base.CONDITION_ORDER[-1]] + 0.38 * style.x_spacing,
        )

        posthoc = stats_df[
            (stats_df["analysis_level"] == "posthoc_vs_control")
            & (stats_df["metric"] == metric)
        ]
        if style.show_control_brackets:
            add_control_brackets(ax, positions, posthoc, style)
        else:
            ylim = ax.get_ylim()
            y_span = ylim[1] - ylim[0]
            annotation_y = ylim[1] - 0.10 * y_span
            for condition in base.CONDITION_ORDER[1:]:
                row = posthoc[posthoc["condition"] == condition]
                if row.empty:
                    continue
                ax.text(
                    positions[condition],
                    annotation_y,
                    base.p_to_star(float(row.iloc[0]["p_adjusted"])),
                    ha="center",
                    va="bottom",
                    fontsize=style.font_size,
                    color="#202020",
                )
        if style.show_omnibus:
            omnibus = stats_df[
                (stats_df["analysis_level"] == "omnibus") & (stats_df["metric"] == metric)
            ].iloc[0]
            ax.text(
                0.02,
                0.98,
                f"{omnibus['test']}: {base.p_to_text(float(omnibus['p_value']))}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=style.font_size * 0.78,
            )
    save_gallery_figure(fig, output_base)
    plt.close(fig)


def plot_paired_variant(region: pd.DataFrame, style: StyleSpec, output_base: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=style.paired_size, constrained_layout=True)
    axes = axes.ravel()
    region_order = ["low_tmrm_outside_mask", "high_tmrm"]
    offsets = {
        "low_tmrm_outside_mask": -0.16 * style.x_spacing,
        "high_tmrm": 0.16 * style.x_spacing,
    }
    positions = condition_positions(style)
    for panel_label, ax, metric in zip(["a", "b", "c", "d"], axes, base.METRIC_ORDER):
        if style.show_panel_labels:
            add_panel_label(ax, panel_label)
        metric_df = region[region["metric"] == metric]
        for condition in base.CONDITION_ORDER:
            condition_df = metric_df[metric_df["condition"] == condition]
            x = positions[condition]
            color = style.palette[condition]
            marker = style.markers[condition]
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
                        alpha=0.35,
                        zorder=1,
                    )
            for region_name in region_order:
                values = condition_df.loc[condition_df["region"] == region_name, "mean"].to_numpy(dtype=float)
                xs = np.full(values.size, x + offsets[region_name])
                face = "white" if region_name == "low_tmrm_outside_mask" else color
                ax.scatter(
                    xs,
                    values,
                    s=style.dot_size * 0.78,
                    marker=marker,
                    facecolor=face,
                    edgecolor=color,
                    linewidth=0.75,
                    alpha=style.dot_alpha,
                    zorder=3,
                )
                mean, ci = base.mean_ci(values)
                ax.errorbar(
                    x + offsets[region_name],
                    mean,
                    yerr=ci,
                    fmt="o",
                    color="#202020",
                    markerfacecolor=face,
                    markeredgecolor="#202020",
                    markersize=style.mean_marker_size * 0.9,
                    elinewidth=style.ci_lw,
                    capsize=style.capsize,
                    zorder=4,
                )
        if style.show_metric_titles:
            ax.set_title(base.METRIC_LABELS[metric], fontweight="bold", fontsize=style.title_size)
        ax.set_xticks([positions[c] for c in base.CONDITION_ORDER])
        ax.set_xticklabels(
            [base.CONDITION_LABELS[c] for c in base.CONDITION_ORDER],
            rotation=30,
            ha="right",
        )
        ax.set_ylabel(base.REGION_YLABELS[metric] if style.metric_labels_on_y else "NADPH region mean")
        ax.margins(x=0.08, y=0.14)
        ax.set_xlim(
            positions[base.CONDITION_ORDER[0]] - 0.42 * style.x_spacing,
            positions[base.CONDITION_ORDER[-1]] + 0.42 * style.x_spacing,
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
    save_gallery_figure(fig, output_base)
    plt.close(fig)


def write_gallery_index(gallery_dir: Path, style_rows: list[dict[str, str]]) -> None:
    cards = []
    for row in style_rows:
        cards.append(
            f"""
            <section>
              <h2>{html.escape(row['label'])}</h2>
              <p>{html.escape(row['description'])}</p>
              <div class="grid">
                <figure>
                  <img src="{html.escape(row['delta_png'])}" alt="{html.escape(row['label'])} delta">
                  <figcaption>Delta plot</figcaption>
                </figure>
                <figure>
                  <img src="{html.escape(row['paired_png'])}" alt="{html.escape(row['label'])} paired">
                  <figcaption>Paired high/low plot</figcaption>
                </figure>
              </div>
            </section>
            """
        )
    condition_text = " -> ".join(base.CONDITION_LABELS[c] for c in base.CONDITION_ORDER)
    page = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>NAD(P)H percentile-p75 style gallery</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #202020; line-height: 1.45; }}
    h1 {{ margin-bottom: 0.2rem; }}
    h2 {{ margin-top: 2rem; }}
    .grid {{ display: grid; grid-template-columns: repeat(2, minmax(280px, 1fr)); gap: 18px; }}
    img {{ width: 100%; border: 1px solid #d7d7d7; }}
    figcaption {{ font-size: 13px; color: #555; margin-top: 4px; }}
    code {{ background: #f3f3f3; padding: 2px 4px; border-radius: 3px; }}
  </style>
</head>
<body>
  <h1>NAD(P)H percentile-p75 style gallery</h1>
  <p>All variants use the same source data, statistics, and condition order:
  <code>{html.escape(condition_text)}</code>.</p>
  {''.join(cards)}
</body>
</html>
"""
    (gallery_dir / "style_gallery_index.html").write_text(page, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=Path.cwd(), help="Analysis folder or project root.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Gallery output folder. Defaults to <analysis-dir>/figures_percentile_p75/style_gallery.",
    )
    parser.add_argument("--method", default=base.METHOD, help="Mask method to plot.")
    parser.add_argument("--condition-order", nargs="+", default=base.CONDITION_ORDER, help="Condition keys in plotting order.")
    parser.add_argument("--condition-label", action="append", default=[], help="Override label as key=Label.")
    parser.add_argument("--condition-color", action="append", default=[], help="Accepted for parity with the main script.")
    parser.add_argument(
        "--styles",
        nargs="+",
        default=[style.name for style in STYLES],
        help="Style names to export.",
    )
    return parser.parse_args()


def configure(args: argparse.Namespace) -> tuple[Path, Path, str, list[StyleSpec]]:
    base_args = argparse.Namespace(
        analysis_dir=args.analysis_dir,
        output_dir=None,
        method=args.method,
        condition_order=args.condition_order,
        condition_label=args.condition_label,
        condition_color=args.condition_color,
    )
    input_dir, base_output_dir, method = base.configure_from_args(base_args)
    gallery_dir = args.output_dir.resolve() if args.output_dir else base_output_dir / "style_gallery"
    style_by_name = {style.name: style for style in STYLES}
    missing = sorted(set(args.styles) - set(style_by_name))
    if missing:
        raise ValueError(f"Unknown style names: {missing}. Available: {sorted(style_by_name)}")
    selected_styles = [style_by_name[name] for name in args.styles]
    return input_dir, gallery_dir, method, selected_styles


def main() -> None:
    args = parse_args()
    input_dir, gallery_dir, method, selected_styles = configure(args)
    gallery_dir.mkdir(parents=True, exist_ok=True)
    delta, region, source = base.load_source_data(input_dir, method)
    stats_df = base.compute_statistics(delta)
    source.to_csv(gallery_dir / "source_data_nadph_percentile_p75.csv", index=False)
    stats_df.to_csv(gallery_dir / "stats_nadph_percentile_p75.csv", index=False)

    style_rows = []
    for style in selected_styles:
        apply_style(style)
        delta_base = gallery_dir / f"{style.name}_delta_percentile_p75"
        paired_base = gallery_dir / f"{style.name}_paired_percentile_p75"
        plot_delta_variant(delta, stats_df, style, delta_base)
        plot_paired_variant(region, style, paired_base)
        style_rows.append(
            {
                "name": style.name,
                "label": style.label,
                "description": style.description,
                "delta_png": delta_base.with_suffix(".png").name,
                "paired_png": paired_base.with_suffix(".png").name,
            }
        )

    pd.DataFrame(style_rows).to_csv(gallery_dir / "style_gallery_manifest.csv", index=False)
    write_gallery_index(gallery_dir, style_rows)
    print(f"Wrote style gallery to {gallery_dir}")


if __name__ == "__main__":
    main()
