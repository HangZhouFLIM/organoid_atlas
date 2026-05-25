#!/usr/bin/env python
"""Analyze NADPH intensity/lifetime in high- and low-TMRM regions.

Default run:
    python analyze_tmrm_nadph.py

Optional final selected-mask run:
    python analyze_tmrm_nadph.py --final-method otsu_clean
"""

from __future__ import annotations

import argparse
import html
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage, stats
from skimage import filters, measure, morphology, segmentation


IMAGE_SHAPE = (246, 246)
NADPH_METRICS = ["Int", "TauPhase", "TauModulation", "TauNormalized"]
MASK_METHODS = ["percentile_p75", "otsu_clean", "adaptive_clean"]
PAIR_RE = re.compile(r"^.+?_(?P<sample>\d+)_(?P<acq>\d+)_pxData\.csv$")


@dataclass(frozen=True)
class Pair:
    sample: int
    tmrm_csv: Path
    nadph_csv: Path


@dataclass
class MaskResult:
    method: str
    mask: np.ndarray
    threshold: str
    initial_pixels: int
    final_pixels: int
    initial_objects: int
    final_objects: int
    removed_small_objects: int
    removed_small_pixels: int
    filled_hole_pixels: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create TMRM mask QC and NADPH high-vs-low region summaries."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Experiment root containing TMRM/ and NADPH/ folders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output folder for QC report and CSV tables.",
    )
    parser.add_argument(
        "--final-method",
        choices=MASK_METHODS,
        default=None,
        help="Optional selected mask method for final chosen-mask statistics.",
    )
    parser.add_argument(
        "--min-object-size",
        type=int,
        default=20,
        help="Remove high-TMRM objects smaller than this many pixels.",
    )
    parser.add_argument(
        "--hole-size",
        type=int,
        default=20,
        help="Fill mask holes smaller than this many pixels.",
    )
    parser.add_argument(
        "--adaptive-sigma",
        type=float,
        default=8.0,
        help="Gaussian sigma for valid-pixel-aware adaptive threshold background.",
    )
    parser.add_argument(
        "--adaptive-std-factor",
        type=float,
        default=0.35,
        help="Adaptive threshold = local mean + factor * local std.",
    )
    return parser.parse_args()


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


def read_pair(pair: Pair) -> tuple[pd.DataFrame, pd.DataFrame]:
    tmrm = pd.read_csv(pair.tmrm_csv)
    nadph = pd.read_csv(pair.nadph_csv)
    required = {"Int", "Ycoord", "Xcoord"}
    required_nadph = required | set(NADPH_METRICS)
    missing_tmrm = required - set(tmrm.columns)
    missing_nadph = required_nadph - set(nadph.columns)
    if missing_tmrm or missing_nadph:
        raise ValueError(
            f"Sample {pair.sample} missing columns: "
            f"TMRM={sorted(missing_tmrm)}, NADPH={sorted(missing_nadph)}"
        )
    if len(tmrm) != len(nadph):
        raise ValueError(
            f"Sample {pair.sample} row count mismatch: TMRM={len(tmrm)}, NADPH={len(nadph)}"
        )
    if not (
        tmrm["Xcoord"].to_numpy().astype(int).tolist()
        == nadph["Xcoord"].to_numpy().astype(int).tolist()
        and tmrm["Ycoord"].to_numpy().astype(int).tolist()
        == nadph["Ycoord"].to_numpy().astype(int).tolist()
    ):
        raise ValueError(f"Sample {pair.sample} TMRM/NADPH coordinates are not row-aligned")
    return tmrm, nadph


def reconstruct_image(df: pd.DataFrame, column: str, shape: tuple[int, int] = IMAGE_SHAPE) -> tuple[np.ndarray, np.ndarray]:
    y = df["Ycoord"].to_numpy(dtype=int) - 1
    x = df["Xcoord"].to_numpy(dtype=int) - 1
    if y.min() < 0 or x.min() < 0 or y.max() >= shape[0] or x.max() >= shape[1]:
        raise ValueError(
            f"Coordinates exceed image shape {shape}: "
            f"X={x.min() + 1}..{x.max() + 1}, Y={y.min() + 1}..{y.max() + 1}"
        )
    image = np.full(shape, np.nan, dtype=float)
    valid = np.zeros(shape, dtype=bool)
    image[y, x] = df[column].to_numpy(dtype=float)
    valid[y, x] = True
    return image, valid


def valid_gaussian(values: np.ndarray, valid: np.ndarray, sigma: float) -> np.ndarray:
    filled = np.where(valid, values, 0.0)
    weights = valid.astype(float)
    numerator = ndimage.gaussian_filter(filled, sigma=sigma, mode="nearest")
    denominator = ndimage.gaussian_filter(weights, sigma=sigma, mode="nearest")
    out = np.full(values.shape, np.nan, dtype=float)
    np.divide(numerator, denominator, out=out, where=denominator > 1e-8)
    return out


def object_summary(mask: np.ndarray) -> tuple[int, np.ndarray]:
    labeled = measure.label(mask, connectivity=2)
    sizes = np.bincount(labeled.ravel())
    object_sizes = sizes[1:] if len(sizes) > 1 else np.array([], dtype=int)
    return int(object_sizes.size), object_sizes


def cleanup_mask(
    initial: np.ndarray,
    valid: np.ndarray,
    min_object_size: int,
    hole_size: int,
) -> tuple[np.ndarray, dict[str, int]]:
    initial = (initial & valid).astype(bool)
    initial_objects, object_sizes = object_summary(initial)
    small_sizes = object_sizes[object_sizes < min_object_size]
    without_small = morphology.remove_small_objects(
        initial, min_size=min_object_size, connectivity=2
    )
    without_small &= valid
    filled = morphology.remove_small_holes(
        without_small, area_threshold=hole_size, connectivity=2
    )
    filled &= valid
    filled_hole_pixels = int(filled.sum() - without_small.sum())
    closed = morphology.binary_closing(filled, morphology.disk(1))
    closed &= valid
    final_objects, _ = object_summary(closed)
    stats_dict = {
        "initial_objects": int(initial_objects),
        "final_objects": int(final_objects),
        "removed_small_objects": int(small_sizes.size),
        "removed_small_pixels": int(small_sizes.sum()) if small_sizes.size else 0,
        "filled_hole_pixels": max(filled_hole_pixels, 0),
    }
    return closed, stats_dict


def make_masks(
    tmrm_image: np.ndarray,
    valid: np.ndarray,
    min_object_size: int,
    hole_size: int,
    adaptive_sigma: float,
    adaptive_std_factor: float,
) -> dict[str, MaskResult]:
    values = tmrm_image[valid]
    masks: dict[str, MaskResult] = {}

    p75 = float(np.nanpercentile(values, 75))
    raw = valid & (tmrm_image >= p75)
    initial_objects, _ = object_summary(raw)
    masks["percentile_p75"] = MaskResult(
        method="percentile_p75",
        mask=raw,
        threshold=f"per-image p75 = {p75:.4g}",
        initial_pixels=int(raw.sum()),
        final_pixels=int(raw.sum()),
        initial_objects=initial_objects,
        final_objects=initial_objects,
        removed_small_objects=0,
        removed_small_pixels=0,
        filled_hole_pixels=0,
    )

    smooth = valid_gaussian(tmrm_image, valid, sigma=1.0)
    otsu = float(filters.threshold_otsu(smooth[valid]))
    otsu_initial = valid & (smooth >= otsu)
    otsu_mask, otsu_stats = cleanup_mask(
        otsu_initial, valid, min_object_size=min_object_size, hole_size=hole_size
    )
    masks["otsu_clean"] = MaskResult(
        method="otsu_clean",
        mask=otsu_mask,
        threshold=f"smoothed Otsu = {otsu:.4g}",
        initial_pixels=int(otsu_initial.sum()),
        final_pixels=int(otsu_mask.sum()),
        **otsu_stats,
    )

    local_mean = valid_gaussian(tmrm_image, valid, sigma=adaptive_sigma)
    local_sq_mean = valid_gaussian(tmrm_image * tmrm_image, valid, sigma=adaptive_sigma)
    local_std = np.sqrt(np.maximum(local_sq_mean - local_mean * local_mean, 0.0))
    adaptive_threshold = local_mean + adaptive_std_factor * local_std
    adaptive_initial = valid & (tmrm_image >= adaptive_threshold)
    adaptive_mask, adaptive_stats = cleanup_mask(
        adaptive_initial, valid, min_object_size=min_object_size, hole_size=hole_size
    )
    masks["adaptive_clean"] = MaskResult(
        method="adaptive_clean",
        mask=adaptive_mask,
        threshold=f"local mean + {adaptive_std_factor:.2f}*local std, sigma={adaptive_sigma:g}",
        initial_pixels=int(adaptive_initial.sum()),
        final_pixels=int(adaptive_mask.sum()),
        **adaptive_stats,
    )
    return masks


def summarize_values(values: np.ndarray) -> dict[str, float | int]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "n_pixels": 0,
            "mean": math.nan,
            "median": math.nan,
            "sd": math.nan,
            "q1": math.nan,
            "q3": math.nan,
            "iqr": math.nan,
        }
    q1, median, q3 = np.nanpercentile(values, [25, 50, 75])
    return {
        "n_pixels": int(values.size),
        "mean": float(np.nanmean(values)),
        "median": float(median),
        "sd": float(np.nanstd(values, ddof=1)) if values.size > 1 else math.nan,
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(q3 - q1),
    }


def summarize_regions(
    sample: int,
    nadph: pd.DataFrame,
    masks: dict[str, MaskResult],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    y = nadph["Ycoord"].to_numpy(dtype=int) - 1
    x = nadph["Xcoord"].to_numpy(dtype=int) - 1
    region_rows: list[dict[str, object]] = []
    diff_rows: list[dict[str, object]] = []
    for method, result in masks.items():
        high_selector = result.mask[y, x]
        low_selector = ~high_selector
        for metric in NADPH_METRICS:
            metric_values = nadph[metric].to_numpy(dtype=float)
            high_summary = summarize_values(metric_values[high_selector])
            low_summary = summarize_values(metric_values[low_selector])
            for region, summary in [("high_tmrm", high_summary), ("low_tmrm_outside_mask", low_summary)]:
                row = {
                    "sample": sample,
                    "method": method,
                    "region": region,
                    "metric": metric,
                }
                row.update(summary)
                region_rows.append(row)
            diff_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "metric": metric,
                    "high_n_pixels": high_summary["n_pixels"],
                    "low_n_pixels": low_summary["n_pixels"],
                    "high_mean": high_summary["mean"],
                    "low_mean": low_summary["mean"],
                    "mean_diff_high_minus_low": high_summary["mean"] - low_summary["mean"],
                    "high_median": high_summary["median"],
                    "low_median": low_summary["median"],
                    "median_diff_high_minus_low": high_summary["median"] - low_summary["median"],
                }
            )
    return region_rows, diff_rows


def paired_tests(diff_df: pd.DataFrame, method: str | None = None) -> pd.DataFrame:
    df = diff_df if method is None else diff_df[diff_df["method"] == method]
    rows: list[dict[str, object]] = []
    for (mask_method, metric), group in df.groupby(["method", "metric"], sort=True):
        diffs = group.sort_values("sample")["mean_diff_high_minus_low"].to_numpy(dtype=float)
        diffs = diffs[np.isfinite(diffs)]
        if diffs.size == 0:
            continue
        t_stat = t_p = w_stat = w_p = math.nan
        if diffs.size >= 2:
            t_res = stats.ttest_1samp(diffs, popmean=0.0, nan_policy="omit")
            t_stat = float(t_res.statistic)
            t_p = float(t_res.pvalue)
        if diffs.size >= 2 and np.any(np.abs(diffs) > 0):
            try:
                w_res = stats.wilcoxon(diffs)
                w_stat = float(w_res.statistic)
                w_p = float(w_res.pvalue)
            except ValueError:
                pass
        rows.append(
            {
                "method": mask_method,
                "metric": metric,
                "n_samples": int(diffs.size),
                "mean_diff_high_minus_low": float(np.mean(diffs)),
                "median_diff_high_minus_low": float(np.median(diffs)),
                "sd_diff": float(np.std(diffs, ddof=1)) if diffs.size > 1 else math.nan,
                "paired_t_stat": t_stat,
                "paired_t_pvalue": t_p,
                "wilcoxon_stat": w_stat,
                "wilcoxon_pvalue": w_p,
            }
        )
    return pd.DataFrame(rows)


def make_qc_figure(
    sample: int,
    tmrm_image: np.ndarray,
    valid: np.ndarray,
    masks: dict[str, MaskResult],
    figure_path: Path,
) -> None:
    masked_tmrm = np.ma.array(tmrm_image, mask=~valid | ~np.isfinite(tmrm_image))
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad(color="#111111")
    vmax = float(np.nanpercentile(tmrm_image[valid], 99))
    vmin = float(np.nanpercentile(tmrm_image[valid], 1))
    fig, axes = plt.subplots(1, 4, figsize=(15, 4), constrained_layout=True)
    axes[0].imshow(masked_tmrm, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    axes[0].set_title(f"Sample {sample}\nTMRM Int")
    axes[0].axis("off")
    for ax, method in zip(axes[1:], MASK_METHODS):
        ax.imshow(masked_tmrm, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
        mask = masks[method].mask
        boundaries = segmentation.find_boundaries(mask, mode="outer")
        overlay = np.zeros((*mask.shape, 4), dtype=float)
        overlay[mask] = [0.0, 1.0, 0.2, 0.22]
        overlay[boundaries] = [0.0, 0.95, 1.0, 1.0]
        ax.imshow(overlay, interpolation="nearest")
        area = 100.0 * mask.sum() / valid.sum()
        ax.set_title(f"{method}\n{area:.1f}% valid pixels")
        ax.axis("off")
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)


def write_html_report(
    output: Path,
    pairs: list[Pair],
    mask_qc_df: pd.DataFrame,
    tests_df: pd.DataFrame,
    final_method: str | None,
) -> None:
    report_path = output / "qc_report.html"
    rows = []
    for pair in pairs:
        sample_qc = mask_qc_df[mask_qc_df["sample"] == pair.sample]
        table_html = sample_qc[
            [
                "method",
                "threshold",
                "valid_pixels",
                "final_high_pixels",
                "final_high_fraction",
                "initial_objects",
                "final_objects",
                "removed_small_objects",
                "removed_small_pixels",
                "filled_hole_pixels",
            ]
        ].to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
        rows.append(
            f"""
            <section class="sample">
              <h2>Sample {pair.sample}</h2>
              <p><code>{html.escape(pair.tmrm_csv.name)}</code> paired with <code>{html.escape(pair.nadph_csv.name)}</code></p>
              <img src="figures/sample_{pair.sample:02d}_mask_qc.png" alt="Sample {pair.sample} mask QC">
              {table_html}
            </section>
            """
        )
    selected = (
        f"<p>Selected final method for final-stat CSVs: <strong>{html.escape(final_method)}</strong>.</p>"
        if final_method
        else "<p>No final method selected yet. After visual QC, rerun with <code>--final-method METHOD</code>.</p>"
    )
    tests_html = tests_df.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>TMRM/NADPH Mask QC</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; line-height: 1.4; color: #1f2933; }}
    h1 {{ margin-bottom: 0.2rem; }}
    h2 {{ margin-top: 2rem; border-top: 1px solid #d9dee7; padding-top: 1rem; }}
    code {{ background: #f3f5f8; padding: 0.1rem 0.25rem; border-radius: 3px; }}
    img {{ max-width: 100%; border: 1px solid #d9dee7; }}
    table {{ border-collapse: collapse; margin: 1rem 0 1.5rem; font-size: 13px; }}
    th, td {{ border: 1px solid #d9dee7; padding: 0.35rem 0.55rem; text-align: right; }}
    th:first-child, td:first-child, td:nth-child(2) {{ text-align: left; }}
    th {{ background: #eef2f7; }}
    .note {{ max-width: 960px; }}
  </style>
</head>
<body>
  <h1>TMRM/NADPH Mask QC</h1>
  <div class="note">
    <p>This report compares candidate high-TMRM masks. Low-TMRM is defined as exported/valid pixels outside each high-TMRM mask.</p>
    {selected}
    <p>CSV outputs are in this folder: <code>mask_qc_summary.csv</code>, <code>nadph_region_summary_all_methods.csv</code>, <code>nadph_sample_differences_all_methods.csv</code>, and <code>paired_tests_all_methods.csv</code>.</p>
  </div>
  <h2>Exploratory Paired Tests Across All Candidate Methods</h2>
  {tests_html}
  {''.join(rows)}
</body>
</html>
"""
    report_path.write_text(report, encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    pairs = discover_pairs(root)
    mask_qc_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    diff_rows: list[dict[str, object]] = []

    for pair in pairs:
        tmrm, nadph = read_pair(pair)
        tmrm_image, valid = reconstruct_image(tmrm, "Int")
        masks = make_masks(
            tmrm_image=tmrm_image,
            valid=valid,
            min_object_size=args.min_object_size,
            hole_size=args.hole_size,
            adaptive_sigma=args.adaptive_sigma,
            adaptive_std_factor=args.adaptive_std_factor,
        )
        valid_pixels = int(valid.sum())
        for method, result in masks.items():
            mask_qc_rows.append(
                {
                    "sample": pair.sample,
                    "method": method,
                    "threshold": result.threshold,
                    "valid_pixels": valid_pixels,
                    "initial_high_pixels": result.initial_pixels,
                    "final_high_pixels": result.final_pixels,
                    "final_high_fraction": result.final_pixels / valid_pixels,
                    "initial_objects": result.initial_objects,
                    "final_objects": result.final_objects,
                    "removed_small_objects": result.removed_small_objects,
                    "removed_small_pixels": result.removed_small_pixels,
                    "filled_hole_pixels": result.filled_hole_pixels,
                    "tmrm_csv": pair.tmrm_csv.name,
                    "nadph_csv": pair.nadph_csv.name,
                }
            )
        sample_region_rows, sample_diff_rows = summarize_regions(pair.sample, nadph, masks)
        region_rows.extend(sample_region_rows)
        diff_rows.extend(sample_diff_rows)
        make_qc_figure(
            sample=pair.sample,
            tmrm_image=tmrm_image,
            valid=valid,
            masks=masks,
            figure_path=figures / f"sample_{pair.sample:02d}_mask_qc.png",
        )

    mask_qc_df = pd.DataFrame(mask_qc_rows)
    region_df = pd.DataFrame(region_rows)
    diff_df = pd.DataFrame(diff_rows)
    tests_df = paired_tests(diff_df)

    mask_qc_df.to_csv(output / "mask_qc_summary.csv", index=False)
    region_df.to_csv(output / "nadph_region_summary_all_methods.csv", index=False)
    diff_df.to_csv(output / "nadph_sample_differences_all_methods.csv", index=False)
    tests_df.to_csv(output / "paired_tests_all_methods.csv", index=False)

    if args.final_method:
        final_region = region_df[region_df["method"] == args.final_method]
        final_diff = diff_df[diff_df["method"] == args.final_method]
        final_tests = paired_tests(diff_df, method=args.final_method)
        final_region.to_csv(output / f"final_region_summary_{args.final_method}.csv", index=False)
        final_diff.to_csv(output / f"final_sample_differences_{args.final_method}.csv", index=False)
        final_tests.to_csv(output / f"final_paired_tests_{args.final_method}.csv", index=False)

    write_html_report(
        output=output,
        pairs=pairs,
        mask_qc_df=mask_qc_df,
        tests_df=tests_df,
        final_method=args.final_method,
    )

    print(f"Processed {len(pairs)} paired samples.")
    print(f"Wrote QC report: {output / 'qc_report.html'}")
    print("Candidate methods:", ", ".join(MASK_METHODS))


if __name__ == "__main__":
    main()
