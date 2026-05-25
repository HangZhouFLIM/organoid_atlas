#!/usr/bin/env python
"""Analyze TMRM, LGR5-GFP, and NADPH FLIM exports.

This is a project-local extension of the two-channel TMRM/NADPH workflow.
It pairs lgr5_gfp, TMRM, and NADPH files by sample number, aligns channels by
X/Y coordinates, builds TMRM masks, allocates pixels into LGR5-GFP low,
medium, and high groups, and separates high-lifetime NADPH pixels that are
consistent with oxidized-lipid-like signal.
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
TMRM_METRICS = ["Int", "TauPhase", "TauModulation", "TauNormalized"]
MASK_METHODS = ["percentile_p75", "otsu_clean", "adaptive_clean"]
PAIR_RE = re.compile(r"^.+?_(?P<sample>\d+)_(?P<acq>\d+)_pxData\.csv$")
CHANNEL_DIRS = {"gfp": "lgr5_gfp", "tmrm": "TMRM", "nadph": "NADPH"}
CONTRASTS = [
    (
        "within_gfp_high_high_tmrm_minus_low_tmrm",
        "high_tmrm_lgr5_high",
        "low_tmrm_lgr5_high",
    ),
    (
        "within_gfp_medium_high_tmrm_minus_low_tmrm",
        "high_tmrm_lgr5_medium",
        "low_tmrm_lgr5_medium",
    ),
    (
        "within_gfp_low_high_tmrm_minus_low_tmrm",
        "high_tmrm_lgr5_low",
        "low_tmrm_lgr5_low",
    ),
]


@dataclass(frozen=True)
class Triad:
    condition: str
    sample: int
    gfp_csv: Path
    tmrm_csv: Path
    nadph_csv: Path


@dataclass
class MaskResult:
    method: str
    mask: np.ndarray
    threshold: str
    threshold_value: float
    initial_pixels: int
    final_pixels: int
    initial_objects: int
    final_objects: int
    removed_small_objects: int
    removed_small_pixels: int
    filled_hole_pixels: int


@dataclass
class GfpGroupResult:
    method: str
    low_mask: np.ndarray
    intermediate_mask: np.ndarray
    high_mask: np.ndarray
    p25: float
    p75: float
    raw_low_pixels: int
    raw_high_pixels: int
    low_pixels: int
    intermediate_pixels: int
    high_pixels: int
    high_initial_objects: int
    high_final_objects: int
    high_removed_small_objects: int
    high_removed_small_pixels: int
    high_filled_hole_pixels: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze tri-channel TMRM/LGR5-GFP/NADPH pxData exports."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Root containing condition folders with lgr5_gfp, TMRM, and NADPH subfolders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output folder. Defaults to ROOT/analysis_TMRM_LGR5_NADPH.",
    )
    parser.add_argument(
        "--final-tmrm-method",
        choices=MASK_METHODS,
        default="percentile_p75",
        help="TMRM mask used for final reports and figures.",
    )
    parser.add_argument(
        "--tmrm-min-object-size",
        type=int,
        default=20,
        help="Remove high-TMRM objects smaller than this many pixels.",
    )
    parser.add_argument(
        "--tmrm-hole-size",
        type=int,
        default=20,
        help="Fill high-TMRM mask holes smaller than this many pixels.",
    )
    parser.add_argument(
        "--gfp-min-object-size",
        type=int,
        default=20,
        help="Remove LGR5-GFP high seed objects smaller than this many pixels.",
    )
    parser.add_argument(
        "--gfp-hole-size",
        type=int,
        default=20,
        help="Fill LGR5-GFP high cell-shape mask holes smaller than this many pixels.",
    )
    parser.add_argument(
        "--gfp-grouping",
        choices=["p25_p75", "p75_only"],
        default="p25_p75",
        help="Use p25/p75 low-medium-high GFP groups, or p75-only high/rest grouping.",
    )
    parser.add_argument(
        "--adaptive-sigma",
        type=float,
        default=8.0,
        help="Gaussian sigma for valid-pixel-aware adaptive TMRM threshold background.",
    )
    parser.add_argument(
        "--adaptive-std-factor",
        type=float,
        default=0.35,
        help="Adaptive TMRM threshold = local mean + factor * local std.",
    )
    parser.add_argument(
        "--lipid-taunormalized-threshold",
        type=float,
        default=4.0,
        help="NADPH TauNormalized above this value is reported as oxidized-lipid-like.",
    )
    parser.add_argument(
        "--lipid-min-object-size",
        type=int,
        default=5,
        help="Minimum connected TauNormalized-high object area, in pixels, for lipid-droplet QC.",
    )
    parser.add_argument(
        "--skip-lipid-analysis",
        action="store_true",
        help="Do not separate TauNormalized-high oxidized-lipid-like pixels or count lipid droplets.",
    )
    return parser.parse_args()


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


def discover_condition_folders(root: Path) -> list[Path]:
    if all((root / folder).is_dir() for folder in CHANNEL_DIRS.values()):
        return [root]
    condition_folders = [
        path
        for path in sorted(root.iterdir())
        if path.is_dir() and all((path / folder).is_dir() for folder in CHANNEL_DIRS.values())
    ]
    if not condition_folders:
        expected = ", ".join(CHANNEL_DIRS.values())
        raise ValueError(f"No condition folders with {expected} found under {root}")
    return condition_folders


def discover_triads(root: Path) -> list[Triad]:
    triads: list[Triad] = []
    messages: list[str] = []
    for condition_path in discover_condition_folders(root):
        indexed = {
            channel: index_pxdata(condition_path / folder)
            for channel, folder in CHANNEL_DIRS.items()
        }
        sample_sets = {channel: set(paths) for channel, paths in indexed.items()}
        common = set.intersection(*sample_sets.values())
        missing_messages = []
        for channel, samples in sample_sets.items():
            missing = sorted(set.union(*sample_sets.values()) - samples)
            if missing:
                missing_messages.append(f"missing {channel} for {missing}")
        if missing_messages:
            messages.append(f"{condition_path.name}: " + "; ".join(missing_messages))
        for sample in sorted(common):
            triads.append(
                Triad(
                    condition=condition_path.name,
                    sample=sample,
                    gfp_csv=indexed["gfp"][sample],
                    tmrm_csv=indexed["tmrm"][sample],
                    nadph_csv=indexed["nadph"][sample],
                )
            )
    if messages:
        raise ValueError("Tri-channel sample mismatch. " + " | ".join(messages))
    if not triads:
        raise ValueError(f"No tri-channel pxData CSV files found under {root}")
    return triads


def read_channel(path: Path, required_metrics: list[str]) -> pd.DataFrame:
    required = set(required_metrics) | {"Ycoord", "Xcoord"}
    df = pd.read_csv(path)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    if df[["Ycoord", "Xcoord"]].duplicated().any():
        raise ValueError(f"{path} contains duplicate X/Y coordinates")
    return df[list(required)]


def reconstruct_images(
    df: pd.DataFrame,
    columns: list[str],
    shape: tuple[int, int] = IMAGE_SHAPE,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    y = df["Ycoord"].to_numpy(dtype=int) - 1
    x = df["Xcoord"].to_numpy(dtype=int) - 1
    if y.min() < 0 or x.min() < 0 or y.max() >= shape[0] or x.max() >= shape[1]:
        raise ValueError(
            f"Coordinates exceed image shape {shape}: "
            f"X={x.min() + 1}..{x.max() + 1}, Y={y.min() + 1}..{y.max() + 1}"
        )
    valid = np.zeros(shape, dtype=bool)
    valid[y, x] = True
    images: dict[str, np.ndarray] = {}
    for column in columns:
        image = np.full(shape, np.nan, dtype=float)
        image[y, x] = df[column].to_numpy(dtype=float)
        images[column] = image
    return images, valid


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
    return closed, {
        "initial_objects": int(initial_objects),
        "final_objects": int(final_objects),
        "removed_small_objects": int(small_sizes.size),
        "removed_small_pixels": int(small_sizes.sum()) if small_sizes.size else 0,
        "filled_hole_pixels": max(filled_hole_pixels, 0),
    }


def make_tmrm_masks(
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
        threshold=f"per-sample p75 = {p75:.4g}",
        threshold_value=p75,
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
        threshold_value=otsu,
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
        threshold_value=math.nan,
        initial_pixels=int(adaptive_initial.sum()),
        final_pixels=int(adaptive_mask.sum()),
        **adaptive_stats,
    )
    return masks


def make_gfp_groups(
    gfp_image: np.ndarray,
    valid: np.ndarray,
    min_object_size: int,
    hole_size: int,
    grouping: str = "p25_p75",
) -> GfpGroupResult:
    values = gfp_image[valid]
    p25 = float(np.nanpercentile(values, 25))
    p75 = float(np.nanpercentile(values, 75))
    raw_low = valid & (gfp_image <= p25)
    raw_high = valid & (gfp_image >= p75)
    high_mask, high_stats = cleanup_mask(
        raw_high, valid, min_object_size=min_object_size, hole_size=hole_size
    )
    if grouping == "p75_only":
        low_mask = valid & ~high_mask
        intermediate_mask = np.zeros_like(valid, dtype=bool)
        method = "gfp_high_p75_rest"
        raw_low_pixels = 0
    else:
        low_mask = raw_low & ~high_mask
        intermediate_mask = valid & ~(low_mask | high_mask)
        method = "gfp_low_medium_high_p25_p75_high_clean"
        raw_low_pixels = int(raw_low.sum())
    return GfpGroupResult(
        method=method,
        low_mask=low_mask,
        intermediate_mask=intermediate_mask,
        high_mask=high_mask,
        p25=p25,
        p75=p75,
        raw_low_pixels=raw_low_pixels,
        raw_high_pixels=int(raw_high.sum()),
        low_pixels=int(low_mask.sum()),
        intermediate_pixels=int(intermediate_mask.sum()),
        high_pixels=int(high_mask.sum()),
        high_initial_objects=high_stats["initial_objects"],
        high_final_objects=high_stats["final_objects"],
        high_removed_small_objects=high_stats["removed_small_objects"],
        high_removed_small_pixels=high_stats["removed_small_pixels"],
        high_filled_hole_pixels=high_stats["filled_hole_pixels"],
    )


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
            "p95": math.nan,
            "max": math.nan,
        }
    q1, median, q3, p95 = np.nanpercentile(values, [25, 50, 75, 95])
    return {
        "n_pixels": int(values.size),
        "mean": float(np.nanmean(values)),
        "median": float(median),
        "sd": float(np.nanstd(values, ddof=1)) if values.size > 1 else math.nan,
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(q3 - q1),
        "p95": float(p95),
        "max": float(np.nanmax(values)),
    }


def sample_rows_for_value(
    condition: str,
    sample: int,
    channel: str,
    metric: str,
    values: np.ndarray,
) -> dict[str, object]:
    row: dict[str, object] = {
        "condition": condition,
        "sample": sample,
        "channel": channel,
        "metric": metric,
    }
    row.update(summarize_values(values))
    return row


def summarize_tmrm_gfp_groups(
    condition: str,
    sample: int,
    tmrm_images: dict[str, np.ndarray],
    gfp_groups: GfpGroupResult,
    valid: np.ndarray,
) -> list[dict[str, object]]:
    masks = {
        "low": valid & gfp_groups.low_mask,
        "medium": valid & gfp_groups.intermediate_mask,
        "high": valid & gfp_groups.high_mask,
        "stem_like": valid & (gfp_groups.intermediate_mask | gfp_groups.high_mask),
    }
    rows: list[dict[str, object]] = []
    for gfp_group, selector in masks.items():
        group_pixels = int(selector.sum())
        for metric in TMRM_METRICS:
            row: dict[str, object] = {
                "condition": condition,
                "sample": sample,
                "gfp_group": gfp_group,
                "metric": metric,
                "group_pixels": group_pixels,
            }
            row.update(summarize_values(tmrm_images[metric][selector]))
            rows.append(row)
    return rows


def build_compartments(
    tmrm_mask: np.ndarray,
    gfp_groups: GfpGroupResult,
    valid: np.ndarray,
) -> dict[str, np.ndarray]:
    high_tmrm = valid & tmrm_mask
    low_tmrm = valid & ~tmrm_mask
    lgr5_high = valid & gfp_groups.high_mask
    lgr5_medium = valid & gfp_groups.intermediate_mask
    lgr5_low = valid & gfp_groups.low_mask
    return {
        "high_tmrm_lgr5_high": high_tmrm & lgr5_high,
        "high_tmrm_lgr5_medium": high_tmrm & lgr5_medium,
        "high_tmrm_lgr5_low": high_tmrm & lgr5_low,
        "low_tmrm_lgr5_high": low_tmrm & lgr5_high,
        "low_tmrm_lgr5_medium": low_tmrm & lgr5_medium,
        "low_tmrm_lgr5_low": low_tmrm & lgr5_low,
    }


def summarize_compartments(
    condition: str,
    sample: int,
    nadph_images: dict[str, np.ndarray],
    tmrm_masks: dict[str, MaskResult],
    gfp_groups: GfpGroupResult,
    lipid_mask: np.ndarray,
    valid: np.ndarray,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    region_rows: list[dict[str, object]] = []
    lipid_rows: list[dict[str, object]] = []
    for method, tmrm_result in tmrm_masks.items():
        compartments = build_compartments(tmrm_result.mask, gfp_groups, valid)
        for compartment, compartment_mask in compartments.items():
            compartment_pixels = int(compartment_mask.sum())
            lipid_selector = compartment_mask & lipid_mask
            real_selector = compartment_mask & ~lipid_mask
            lipid_pixels = int(lipid_selector.sum())
            lipid_rows.append(
                {
                    "condition": condition,
                    "sample": sample,
                    "method": method,
                    "compartment": compartment,
                    "compartment_pixels": compartment_pixels,
                    "real_nadph_pixels": int(real_selector.sum()),
                    "oxidized_lipid_like_pixels": lipid_pixels,
                    "oxidized_lipid_like_fraction": (
                        lipid_pixels / compartment_pixels if compartment_pixels else math.nan
                    ),
                    "lipid_tau_normalized_threshold": math.nan,
                }
            )
            for nadph_class, selector in [
                ("real_nadph", real_selector),
                ("oxidized_lipid_like", lipid_selector),
            ]:
                for metric in NADPH_METRICS:
                    row: dict[str, object] = {
                        "condition": condition,
                        "sample": sample,
                        "method": method,
                        "compartment": compartment,
                        "nadph_class": nadph_class,
                        "metric": metric,
                        "compartment_pixels": compartment_pixels,
                    }
                    row.update(summarize_values(nadph_images[metric][selector]))
                    region_rows.append(row)
    return region_rows, lipid_rows


def summarize_lipid_droplets(
    condition: str,
    sample: int,
    nadph_images: dict[str, np.ndarray],
    valid: np.ndarray,
    lipid_mask: np.ndarray,
    min_object_size: int,
) -> tuple[list[dict[str, object]], dict[str, object], np.ndarray]:
    labeled = measure.label(lipid_mask, connectivity=2)
    droplet_rows: list[dict[str, object]] = []
    kept_mask = np.zeros(lipid_mask.shape, dtype=bool)
    droplet_id = 0

    for prop in measure.regionprops(labeled):
        area = int(prop.area)
        if area < min_object_size:
            continue
        droplet_id += 1
        object_mask = labeled == prop.label
        kept_mask |= object_mask
        min_row, min_col, max_row, max_col = prop.bbox
        row: dict[str, object] = {
            "condition": condition,
            "sample": sample,
            "droplet_id": droplet_id,
            "label_id": int(prop.label),
            "area_pixels": area,
            "centroid_y": float(prop.centroid[0]),
            "centroid_x": float(prop.centroid[1]),
            "bbox_min_y": int(min_row),
            "bbox_min_x": int(min_col),
            "bbox_max_y": int(max_row),
            "bbox_max_x": int(max_col),
        }
        for metric in NADPH_METRICS:
            values = nadph_images[metric][object_mask]
            row[f"mean_{metric}"] = float(np.nanmean(values)) if values.size else math.nan
            row[f"median_{metric}"] = float(np.nanmedian(values)) if values.size else math.nan
        droplet_rows.append(row)

    areas = np.array([row["area_pixels"] for row in droplet_rows], dtype=float)
    valid_pixels = int(valid.sum())
    candidate_pixels = int(lipid_mask.sum())
    droplet_pixels = int(kept_mask.sum())
    sample_row: dict[str, object] = {
        "condition": condition,
        "sample": sample,
        "valid_common_pixels": valid_pixels,
        "lipid_tau_normalized_threshold": math.nan,
        "lipid_min_object_size": min_object_size,
        "candidate_lipid_pixels": candidate_pixels,
        "candidate_lipid_fraction": candidate_pixels / valid_pixels if valid_pixels else math.nan,
        "lipid_droplet_pixels": droplet_pixels,
        "lipid_droplet_area_fraction": droplet_pixels / valid_pixels if valid_pixels else math.nan,
        "n_lipid_droplets": len(droplet_rows),
        "lipid_droplet_density_per_10000_valid_pixels": (
            len(droplet_rows) / valid_pixels * 10000 if valid_pixels else math.nan
        ),
        "mean_droplet_area_pixels": float(np.nanmean(areas)) if areas.size else math.nan,
        "median_droplet_area_pixels": float(np.nanmedian(areas)) if areas.size else math.nan,
        "max_droplet_area_pixels": float(np.nanmax(areas)) if areas.size else math.nan,
    }
    if droplet_rows:
        droplet_df = pd.DataFrame(droplet_rows)
        for metric in NADPH_METRICS:
            sample_row[f"mean_droplet_{metric}"] = float(droplet_df[f"mean_{metric}"].mean())
            sample_row[f"median_droplet_{metric}"] = float(droplet_df[f"median_{metric}"].median())
    else:
        for metric in NADPH_METRICS:
            sample_row[f"mean_droplet_{metric}"] = math.nan
            sample_row[f"median_droplet_{metric}"] = math.nan
    return droplet_rows, sample_row, kept_mask


def summarize_sample_contrasts(region_df: pd.DataFrame) -> pd.DataFrame:
    real = region_df[region_df["nadph_class"] == "real_nadph"]
    rows: list[dict[str, object]] = []
    group_cols = ["condition", "sample", "method", "metric"]
    for (condition, sample, method, metric), group in real.groupby(group_cols, sort=True):
        by_compartment = group.set_index("compartment")
        for contrast, group_a, group_b in CONTRASTS:
            if group_a not in by_compartment.index or group_b not in by_compartment.index:
                continue
            a = by_compartment.loc[group_a]
            b = by_compartment.loc[group_b]
            rows.append(
                {
                    "condition": condition,
                    "sample": sample,
                    "method": method,
                    "metric": metric,
                    "contrast": contrast,
                    "group_a": group_a,
                    "group_b": group_b,
                    "group_a_n_pixels": int(a["n_pixels"]),
                    "group_b_n_pixels": int(b["n_pixels"]),
                    "group_a_mean": float(a["mean"]),
                    "group_b_mean": float(b["mean"]),
                    "mean_diff_group_a_minus_b": float(a["mean"] - b["mean"]),
                    "group_a_median": float(a["median"]),
                    "group_b_median": float(b["median"]),
                    "median_diff_group_a_minus_b": float(a["median"] - b["median"]),
                }
            )
    return pd.DataFrame(rows)


def paired_tests(contrast_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    test_frames = [(condition, group) for condition, group in contrast_df.groupby("condition")]
    test_frames.append(("ALL", contrast_df))
    for condition, frame in test_frames:
        for (method, metric, contrast), group in frame.groupby(["method", "metric", "contrast"], sort=True):
            diffs = group.sort_values(["condition", "sample"])[
                "mean_diff_group_a_minus_b"
            ].to_numpy(dtype=float)
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
                    "condition": condition,
                    "method": method,
                    "metric": metric,
                    "contrast": contrast,
                    "n_samples": int(diffs.size),
                    "mean_diff_group_a_minus_b": float(np.mean(diffs)),
                    "median_diff_group_a_minus_b": float(np.median(diffs)),
                    "sd_diff": float(np.std(diffs, ddof=1)) if diffs.size > 1 else math.nan,
                    "paired_t_stat": t_stat,
                    "paired_t_pvalue": t_p,
                    "wilcoxon_stat": w_stat,
                    "wilcoxon_pvalue": w_p,
                }
            )
    return pd.DataFrame(rows)


def rgba_overlay(mask: np.ndarray, rgba: tuple[float, float, float, float]) -> np.ndarray:
    out = np.zeros((*mask.shape, 4), dtype=float)
    out[mask] = rgba
    return out


def percentile_limits(image: np.ndarray, valid: np.ndarray) -> tuple[float, float]:
    values = image[valid & np.isfinite(image)]
    if values.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.nanpercentile(values, [1, 99])
    if vmin == vmax:
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


def make_qc_figure(
    triad: Triad,
    tmrm_image: np.ndarray,
    gfp_image: np.ndarray,
    nadph_images: dict[str, np.ndarray],
    valid: np.ndarray,
    tmrm_masks: dict[str, MaskResult],
    gfp_groups: GfpGroupResult,
    lipid_mask: np.ndarray,
    lipid_droplet_mask: np.ndarray,
    final_method: str,
    figure_path: Path,
) -> None:
    tmrm_mask = tmrm_masks[final_method].mask
    compartments = build_compartments(tmrm_mask, gfp_groups, valid)
    raw_low = valid & (gfp_image <= gfp_groups.p25)
    raw_high = valid & (gfp_image >= gfp_groups.p75)
    removed_high = raw_high & ~gfp_groups.high_mask
    filled_high = gfp_groups.high_mask & ~raw_high
    stem_like = gfp_groups.intermediate_mask | gfp_groups.high_mask
    p75_only = gfp_groups.method == "gfp_high_p75_rest"

    fig, axes = plt.subplots(3, 4, figsize=(18, 12), constrained_layout=True)
    cmap_tmrm = plt.get_cmap("magma").copy()
    cmap_gfp = plt.get_cmap("viridis").copy()
    cmap_nadph = plt.get_cmap("plasma").copy()
    for cmap in [cmap_tmrm, cmap_gfp, cmap_nadph]:
        cmap.set_bad(color="#111111")

    tmrm_vmin, tmrm_vmax = percentile_limits(tmrm_image, valid)
    gfp_vmin, gfp_vmax = percentile_limits(gfp_image, valid)
    tau_vmin, tau_vmax = percentile_limits(nadph_images["TauNormalized"], valid)

    tmrm_ma = np.ma.array(tmrm_image, mask=~valid | ~np.isfinite(tmrm_image))
    gfp_ma = np.ma.array(gfp_image, mask=~valid | ~np.isfinite(gfp_image))
    tau_ma = np.ma.array(nadph_images["TauNormalized"], mask=~valid | ~np.isfinite(nadph_images["TauNormalized"]))

    axes[0, 0].imshow(tmrm_ma, cmap=cmap_tmrm, vmin=tmrm_vmin, vmax=tmrm_vmax, interpolation="nearest")
    axes[0, 0].set_title("TMRM Int")
    axes[0, 0].axis("off")

    axes[0, 1].imshow(tmrm_ma, cmap=cmap_tmrm, vmin=tmrm_vmin, vmax=tmrm_vmax, interpolation="nearest")
    axes[0, 1].imshow(rgba_overlay(tmrm_mask, (0.0, 1.0, 0.2, 0.22)), interpolation="nearest")
    axes[0, 1].imshow(rgba_overlay(segmentation.find_boundaries(tmrm_mask, mode="outer"), (0.0, 0.95, 1.0, 1.0)), interpolation="nearest")
    axes[0, 1].set_title(f"TMRM high: {final_method}")
    axes[0, 1].axis("off")

    axes[0, 2].imshow(gfp_ma, cmap=cmap_gfp, vmin=gfp_vmin, vmax=gfp_vmax, interpolation="nearest")
    axes[0, 2].imshow(rgba_overlay(segmentation.find_boundaries(gfp_groups.high_mask, mode="outer"), (1.0, 0.95, 0.0, 1.0)), interpolation="nearest")
    axes[0, 2].set_title("LGR5-GFP Int + cleaned high outline")
    axes[0, 2].axis("off")

    gfp_group_rgb = np.zeros((*valid.shape, 4), dtype=float)
    gfp_group_rgb[valid] = (0.0, 0.0, 0.0, 0.15)
    gfp_group_rgb[gfp_groups.low_mask] = (0.16, 0.38, 0.76, 0.78)
    gfp_group_rgb[gfp_groups.intermediate_mask] = (0.28, 0.62, 0.35, 0.78)
    gfp_group_rgb[gfp_groups.high_mask] = (0.95, 0.63, 0.13, 0.88)
    axes[0, 3].imshow(gfp_group_rgb, interpolation="nearest")
    axes[0, 3].imshow(rgba_overlay(segmentation.find_boundaries(gfp_groups.high_mask, mode="outer"), (1.0, 1.0, 1.0, 1.0)), interpolation="nearest")
    axes[0, 3].set_title(
        "GFP groups: rest / cleaned high" if p75_only else "GFP groups: low / medium / cleaned high"
    )
    axes[0, 3].axis("off")

    raw_seed_rgb = np.zeros((*valid.shape, 4), dtype=float)
    raw_seed_rgb[valid] = (0.0, 0.0, 0.0, 0.12)
    raw_seed_rgb[raw_low] = (0.16, 0.38, 0.76, 0.35)
    raw_seed_rgb[raw_high] = (0.92, 0.22, 0.18, 0.78)
    axes[1, 0].imshow(raw_seed_rgb, interpolation="nearest")
    axes[1, 0].set_title("Raw GFP high seed: >=p75" if p75_only else "Raw GFP tails: <=p25 and >=p75 seed")
    axes[1, 0].axis("off")

    cleaned_rgb = np.zeros((*valid.shape, 4), dtype=float)
    cleaned_rgb[valid] = (0.0, 0.0, 0.0, 0.12)
    cleaned_rgb[gfp_groups.high_mask] = (0.95, 0.63, 0.13, 0.88)
    cleaned_rgb[removed_high] = (0.92, 0.22, 0.18, 0.70)
    cleaned_rgb[filled_high] = (1.0, 0.95, 0.0, 0.75)
    axes[1, 1].imshow(cleaned_rgb, interpolation="nearest")
    axes[1, 1].set_title("Cleaned high GFP: kept / removed / filled")
    axes[1, 1].axis("off")

    axes[1, 2].imshow(gfp_ma, cmap=cmap_gfp, vmin=gfp_vmin, vmax=gfp_vmax, interpolation="nearest")
    axes[1, 2].imshow(rgba_overlay(gfp_groups.high_mask, (0.95, 0.63, 0.13, 0.35)), interpolation="nearest")
    axes[1, 2].imshow(rgba_overlay(segmentation.find_boundaries(gfp_groups.high_mask, mode="outer"), (1.0, 1.0, 1.0, 1.0)), interpolation="nearest")
    axes[1, 2].set_title("Cleaned high GFP cell-shape check")
    axes[1, 2].axis("off")

    stem_rgb = np.zeros((*valid.shape, 4), dtype=float)
    stem_rgb[valid] = (0.0, 0.0, 0.0, 0.12)
    stem_rgb[gfp_groups.low_mask] = (0.16, 0.38, 0.76, 0.45)
    stem_rgb[gfp_groups.intermediate_mask] = (0.28, 0.62, 0.35, 0.72)
    stem_rgb[gfp_groups.high_mask] = (0.95, 0.63, 0.13, 0.85)
    axes[1, 3].imshow(stem_rgb, interpolation="nearest")
    axes[1, 3].imshow(rgba_overlay(segmentation.find_boundaries(stem_like, mode="outer"), (1.0, 1.0, 1.0, 1.0)), interpolation="nearest")
    axes[1, 3].set_title("Stem-like GFP for plots: high only" if p75_only else "Stem-like GFP for plots: medium + high")
    axes[1, 3].axis("off")

    axes[2, 0].imshow(tau_ma, cmap=cmap_nadph, vmin=tau_vmin, vmax=tau_vmax, interpolation="nearest")
    axes[2, 0].set_title("NADPH TauNormalized")
    axes[2, 0].axis("off")

    lipid_candidate_only = lipid_mask & ~lipid_droplet_mask
    axes[2, 1].imshow(tau_ma, cmap=cmap_nadph, vmin=tau_vmin, vmax=tau_vmax, interpolation="nearest")
    axes[2, 1].imshow(rgba_overlay(lipid_candidate_only, (0.75, 0.75, 0.75, 0.35)), interpolation="nearest")
    axes[2, 1].imshow(rgba_overlay(lipid_droplet_mask, (0.0, 0.85, 1.0, 0.52)), interpolation="nearest")
    axes[2, 1].imshow(rgba_overlay(segmentation.find_boundaries(lipid_droplet_mask, mode="outer"), (1.0, 1.0, 1.0, 1.0)), interpolation="nearest")
    axes[2, 1].set_title("Lipid-like droplets: TauN > 4, objects >=5 px")
    axes[2, 1].axis("off")

    compartment_rgb = np.zeros((*valid.shape, 4), dtype=float)
    compartment_rgb[compartments["high_tmrm_lgr5_high"]] = (0.95, 0.63, 0.13, 0.95)
    compartment_rgb[compartments["high_tmrm_lgr5_medium"]] = (0.28, 0.62, 0.35, 0.90)
    compartment_rgb[compartments["high_tmrm_lgr5_low"]] = (0.16, 0.38, 0.76, 0.90)
    compartment_rgb[compartments["low_tmrm_lgr5_high"]] = (0.95, 0.63, 0.13, 0.42)
    compartment_rgb[compartments["low_tmrm_lgr5_medium"]] = (0.28, 0.62, 0.35, 0.42)
    compartment_rgb[compartments["low_tmrm_lgr5_low"]] = (0.16, 0.38, 0.76, 0.42)
    axes[2, 2].imshow(np.ma.array(np.zeros(valid.shape), mask=~valid), cmap="gray", interpolation="nearest")
    axes[2, 2].imshow(compartment_rgb, interpolation="nearest")
    axes[2, 2].imshow(rgba_overlay(segmentation.find_boundaries(tmrm_mask, mode="outer"), (1.0, 1.0, 1.0, 0.85)), interpolation="nearest")
    axes[2, 2].set_title("Final compartments: bright = high TMRM")
    axes[2, 2].axis("off")

    values = gfp_image[valid & np.isfinite(gfp_image)]
    axes[2, 3].hist(values, bins=60, color="#317873", alpha=0.85)
    axes[2, 3].axvline(gfp_groups.p25, color="#3366cc", linewidth=2, label="p25")
    axes[2, 3].axvline(gfp_groups.p75, color="#d8893a", linewidth=2, label="p75")
    axes[2, 3].set_title("GFP Int distribution")
    axes[2, 3].set_xlabel("GFP Int")
    axes[2, 3].set_ylabel("Pixels")
    axes[2, 3].legend(frameon=False)

    for ax in axes.ravel():
        ax.tick_params(labelsize=8)

    fig.text(
        0.015,
        0.01,
        ("GFP colors: rest=blue, cleaned high=orange. " if p75_only else "GFP colors: low=blue, medium=green, cleaned high=orange. ")
        + "High-GFP cleaning: orange=kept, red=removed, yellow=filled. "
        "Lipid droplets: cyan=kept object, gray=small TauNormalized-high candidate.",
        fontsize=9,
        ha="left",
        va="bottom",
    )

    fig.suptitle(
        f"{triad.condition} sample {triad.sample}: coordinate-joined common pixels",
        fontsize=14,
    )
    fig.savefig(figure_path, dpi=170)
    plt.close(fig)


def write_by_condition_tables(output: Path, tables: dict[str, pd.DataFrame], condition_col: str = "condition") -> None:
    by_condition = output / "by_condition"
    by_condition.mkdir(exist_ok=True)
    conditions = sorted(
        {
            str(value)
            for df in tables.values()
            if condition_col in df.columns
            for value in df[condition_col].dropna().unique()
            if str(value) != "ALL"
        }
    )
    for condition in conditions:
        condition_dir = by_condition / condition
        condition_dir.mkdir(exist_ok=True)
        for name, df in tables.items():
            if condition_col in df.columns:
                df[df[condition_col] == condition].to_csv(condition_dir / name, index=False)


def write_html_report(
    output: Path,
    triads: list[Triad],
    join_qc_df: pd.DataFrame,
    tmrm_qc_df: pd.DataFrame,
    gfp_qc_df: pd.DataFrame,
    lipid_summary_df: pd.DataFrame,
    lipid_droplet_sample_df: pd.DataFrame,
    paired_tests_df: pd.DataFrame,
    final_method: str,
    lipid_threshold: float,
    lipid_min_object_size: int,
    lipid_analysis_enabled: bool,
) -> None:
    sections = []
    for triad in triads:
        selector = (join_qc_df["condition"] == triad.condition) & (join_qc_df["sample"] == triad.sample)
        join_html = join_qc_df[selector].to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
        gfp_html = gfp_qc_df[selector].to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
        tmrm_selector = selector & (tmrm_qc_df["method"] == final_method)
        tmrm_html = tmrm_qc_df[tmrm_selector].to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
        lipid_selector = (
            (lipid_summary_df["condition"] == triad.condition)
            & (lipid_summary_df["sample"] == triad.sample)
            & (lipid_summary_df["method"] == final_method)
        )
        lipid_html = lipid_summary_df[lipid_selector][
            [
                "compartment",
                "compartment_pixels",
                "real_nadph_pixels",
                "oxidized_lipid_like_pixels",
                "oxidized_lipid_like_fraction",
            ]
        ].to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
        droplet_selector = (
            (lipid_droplet_sample_df["condition"] == triad.condition)
            & (lipid_droplet_sample_df["sample"] == triad.sample)
        )
        droplet_columns = [
            "valid_common_pixels",
            "candidate_lipid_pixels",
            "candidate_lipid_fraction",
            "lipid_droplet_pixels",
            "lipid_droplet_area_fraction",
            "n_lipid_droplets",
            "lipid_droplet_density_per_10000_valid_pixels",
            "mean_droplet_area_pixels",
            "median_droplet_area_pixels",
            "max_droplet_area_pixels",
        ]
        droplet_html = lipid_droplet_sample_df[droplet_selector][droplet_columns].to_html(
            index=False,
            escape=True,
            float_format=lambda x: f"{x:.4g}",
        )
        image_name = f"figures/{triad.condition}_sample_{triad.sample:02d}_qc.png"
        sections.append(
            f"""
            <section class="sample">
              <h2>{html.escape(triad.condition)} sample {triad.sample}</h2>
              <p>
                GFP: <code>{html.escape(triad.gfp_csv.name)}</code>;
                TMRM: <code>{html.escape(triad.tmrm_csv.name)}</code>;
                NADPH: <code>{html.escape(triad.nadph_csv.name)}</code>
              </p>
              <img src="{html.escape(image_name)}" alt="{html.escape(triad.condition)} sample {triad.sample} QC">
              <h3>Coordinate Join QC</h3>
              {join_html}
              <h3>LGR5-GFP Low/Medium/High QC</h3>
              {gfp_html}
              <h3>Final TMRM Mask QC</h3>
              {tmrm_html}
              <h3>Oxidized-Lipid-Like NADPH Pixels</h3>
              {lipid_html}
              <h3>Oxidized-Lipid-Like Droplet Objects</h3>
              {droplet_html}
            </section>
            """
        )

    p75_only_report = set(gfp_qc_df["method"].astype(str).unique()) == {"gfp_high_p75_rest"}
    if p75_only_report:
        gfp_grouping_text = (
            "LGR5-GFP statistics use p75-only grouping: cleaned high seeded from "
            "<code>GFP Int &gt;= p75</code>, and all valid non-high pixels are treated as rest/low."
        )
        stem_text = "For the main plots, cleaned high GFP is treated as the stem-like compartment."
    else:
        gfp_grouping_text = (
            "LGR5-GFP statistics use three groups for plotting: low "
            "<code>GFP Int &lt;= p25</code>, medium <code>p25 &lt; GFP Int &lt; p75</code>, "
            "and cleaned high seeded from <code>GFP Int &gt;= p75</code>."
        )
        stem_text = (
            "For the main plots, medium + cleaned high GFP is treated as the stem-like compartment, "
            "and low GFP is the low-GFP/non-stem compartment."
        )
    if lipid_analysis_enabled:
        lipid_text = (
            f"Real NADPH summaries exclude all pixels with <code>TauNormalized &gt; {lipid_threshold:.4g}</code>. "
            f"The droplet QC additionally counts connected lipid-like objects from the same mask after keeping objects "
            f"with area <code>&gt;= {lipid_min_object_size}</code> pixels."
        )
    else:
        lipid_text = (
            "Oxidized-lipid-like separation was disabled for this run; "
            "real NADPH summaries include all valid NADPH pixels."
        )

    tests_preview = paired_tests_df[
        (paired_tests_df["condition"] != "ALL")
        & (paired_tests_df["method"] == final_method)
    ].to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}")
    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>TMRM/LGR5-GFP/NADPH QC</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; line-height: 1.4; color: #1f2933; }}
    h1 {{ margin-bottom: 0.2rem; }}
    h2 {{ margin-top: 2rem; border-top: 1px solid #d9dee7; padding-top: 1rem; }}
    h3 {{ margin-top: 1.2rem; }}
    code {{ background: #f3f5f8; padding: 0.1rem 0.25rem; border-radius: 3px; }}
    img {{ max-width: 100%; border: 1px solid #d9dee7; }}
    table {{ border-collapse: collapse; margin: 0.8rem 0 1.4rem; font-size: 12px; }}
    th, td {{ border: 1px solid #d9dee7; padding: 0.32rem 0.48rem; text-align: right; }}
    th:first-child, td:first-child, td:nth-child(2), td:nth-child(3) {{ text-align: left; }}
    th {{ background: #eef2f7; }}
    .note {{ max-width: 980px; }}
  </style>
</head>
<body>
  <h1>TMRM/LGR5-GFP/NADPH QC</h1>
  <div class="note">
    <p>Pixels are aligned by <code>Xcoord,Ycoord</code>. The analysis-valid area is the coordinate intersection across LGR5-GFP, TMRM, and NADPH.</p>
    <p>{gfp_grouping_text}</p>
    <p>{stem_text}</p>
    <p>{lipid_text}</p>
    <p>Final report figures use TMRM method <strong>{html.escape(final_method)}</strong>. All candidate TMRM methods are included in CSV outputs.</p>
  </div>
  <h2>Final-Method Paired Tests by Condition</h2>
  {tests_preview}
  {''.join(sections)}
</body>
</html>
"""
    (output / "qc_report.html").write_text(report, encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    output = (args.output or (root / "analysis_TMRM_LGR5_NADPH")).resolve()
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    triads = discover_triads(root)
    join_qc_rows: list[dict[str, object]] = []
    tmrm_qc_rows: list[dict[str, object]] = []
    gfp_qc_rows: list[dict[str, object]] = []
    sample_distribution_rows: list[dict[str, object]] = []
    tmrm_gfp_group_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    lipid_rows: list[dict[str, object]] = []
    lipid_droplet_rows: list[dict[str, object]] = []
    lipid_droplet_sample_rows: list[dict[str, object]] = []

    for triad in triads:
        gfp_df = read_channel(triad.gfp_csv, ["Int"])
        tmrm_df = read_channel(triad.tmrm_csv, TMRM_METRICS)
        nadph_df = read_channel(triad.nadph_csv, NADPH_METRICS)

        gfp_images, gfp_valid = reconstruct_images(gfp_df, ["Int"])
        tmrm_images, tmrm_valid = reconstruct_images(tmrm_df, TMRM_METRICS)
        nadph_images, nadph_valid = reconstruct_images(nadph_df, NADPH_METRICS)

        common_valid = gfp_valid & tmrm_valid & nadph_valid
        valid_pixels = int(common_valid.sum())
        if valid_pixels == 0:
            raise ValueError(f"{triad.condition} sample {triad.sample} has no common channel pixels")

        join_qc_rows.append(
            {
                "condition": triad.condition,
                "sample": triad.sample,
                "gfp_rows": int(gfp_valid.sum()),
                "tmrm_rows": int(tmrm_valid.sum()),
                "nadph_rows": int(nadph_valid.sum()),
                "common_pixels": valid_pixels,
                "common_fraction_vs_gfp": valid_pixels / int(gfp_valid.sum()),
                "common_fraction_vs_tmrm": valid_pixels / int(tmrm_valid.sum()),
                "common_fraction_vs_nadph": valid_pixels / int(nadph_valid.sum()),
                "gfp_csv": triad.gfp_csv.name,
                "tmrm_csv": triad.tmrm_csv.name,
                "nadph_csv": triad.nadph_csv.name,
            }
        )

        for channel, images in [
            ("lgr5_gfp", gfp_images),
            ("tmrm", tmrm_images),
            ("nadph", nadph_images),
        ]:
            for metric, image in images.items():
                sample_distribution_rows.append(
                    sample_rows_for_value(
                        triad.condition,
                        triad.sample,
                        channel,
                        metric,
                        image[common_valid],
                    )
                )

        tmrm_masks = make_tmrm_masks(
            tmrm_image=tmrm_images["Int"],
            valid=common_valid,
            min_object_size=args.tmrm_min_object_size,
            hole_size=args.tmrm_hole_size,
            adaptive_sigma=args.adaptive_sigma,
            adaptive_std_factor=args.adaptive_std_factor,
        )
        gfp_groups = make_gfp_groups(
            gfp_image=gfp_images["Int"],
            valid=common_valid,
            min_object_size=args.gfp_min_object_size,
            hole_size=args.gfp_hole_size,
            grouping=args.gfp_grouping,
        )
        if args.skip_lipid_analysis:
            lipid_mask = np.zeros_like(common_valid, dtype=bool)
        else:
            lipid_mask = common_valid & (nadph_images["TauNormalized"] > args.lipid_taunormalized_threshold)
        sample_lipid_droplets, sample_lipid_droplet_summary, lipid_droplet_mask = summarize_lipid_droplets(
            condition=triad.condition,
            sample=triad.sample,
            nadph_images=nadph_images,
            valid=common_valid,
            lipid_mask=lipid_mask,
            min_object_size=args.lipid_min_object_size,
        )
        for row in sample_lipid_droplets:
            row["lipid_tau_normalized_threshold"] = (
                math.nan if args.skip_lipid_analysis else args.lipid_taunormalized_threshold
            )
            row["lipid_min_object_size"] = args.lipid_min_object_size
        sample_lipid_droplet_summary["lipid_tau_normalized_threshold"] = (
            math.nan if args.skip_lipid_analysis else args.lipid_taunormalized_threshold
        )
        sample_lipid_droplet_summary["lipid_analysis_enabled"] = not args.skip_lipid_analysis
        lipid_droplet_rows.extend(sample_lipid_droplets)
        lipid_droplet_sample_rows.append(sample_lipid_droplet_summary)

        for method, result in tmrm_masks.items():
            tmrm_qc_rows.append(
                {
                    "condition": triad.condition,
                    "sample": triad.sample,
                    "method": method,
                    "threshold": result.threshold,
                    "valid_common_pixels": valid_pixels,
                    "initial_high_pixels": result.initial_pixels,
                    "final_high_pixels": result.final_pixels,
                    "final_high_fraction": result.final_pixels / valid_pixels,
                    "initial_objects": result.initial_objects,
                    "final_objects": result.final_objects,
                    "removed_small_objects": result.removed_small_objects,
                    "removed_small_pixels": result.removed_small_pixels,
                    "filled_hole_pixels": result.filled_hole_pixels,
                }
            )

        gfp_values = gfp_images["Int"][common_valid]
        gfp_summary = summarize_values(gfp_values)
        gfp_qc_row: dict[str, object] = {
            "condition": triad.condition,
            "sample": triad.sample,
            "method": gfp_groups.method,
            "valid_common_pixels": valid_pixels,
            "p25_threshold": gfp_groups.p25,
            "p75_threshold": gfp_groups.p75,
            "raw_low_pixels": gfp_groups.raw_low_pixels,
            "raw_high_pixels": gfp_groups.raw_high_pixels,
            "low_pixels": gfp_groups.low_pixels,
            "medium_pixels": gfp_groups.intermediate_pixels,
            "cleaned_high_pixels": gfp_groups.high_pixels,
            "low_fraction": gfp_groups.low_pixels / valid_pixels,
            "medium_fraction": gfp_groups.intermediate_pixels / valid_pixels,
            "cleaned_high_fraction": gfp_groups.high_pixels / valid_pixels,
            "high_initial_objects": gfp_groups.high_initial_objects,
            "high_final_objects": gfp_groups.high_final_objects,
            "high_removed_small_objects": gfp_groups.high_removed_small_objects,
            "high_removed_small_pixels": gfp_groups.high_removed_small_pixels,
            "high_filled_hole_pixels": gfp_groups.high_filled_hole_pixels,
            "gfp_min": float(np.nanmin(gfp_values)),
            "gfp_p25": float(np.nanpercentile(gfp_values, 25)),
            "gfp_p50": float(np.nanpercentile(gfp_values, 50)),
            "gfp_p75": float(np.nanpercentile(gfp_values, 75)),
            "gfp_p90": float(np.nanpercentile(gfp_values, 90)),
            "gfp_p95": float(np.nanpercentile(gfp_values, 95)),
            "gfp_p99": float(np.nanpercentile(gfp_values, 99)),
            "gfp_max": float(np.nanmax(gfp_values)),
        }
        gfp_qc_row.update({f"gfp_{key}": value for key, value in gfp_summary.items() if key == "sd"})
        gfp_qc_rows.append(gfp_qc_row)
        tmrm_gfp_group_rows.extend(
            summarize_tmrm_gfp_groups(
                condition=triad.condition,
                sample=triad.sample,
                tmrm_images=tmrm_images,
                gfp_groups=gfp_groups,
                valid=common_valid,
            )
        )

        sample_region_rows, sample_lipid_rows = summarize_compartments(
            condition=triad.condition,
            sample=triad.sample,
            nadph_images=nadph_images,
            tmrm_masks=tmrm_masks,
            gfp_groups=gfp_groups,
            lipid_mask=lipid_mask,
            valid=common_valid,
        )
        for row in sample_lipid_rows:
            row["lipid_tau_normalized_threshold"] = (
                math.nan if args.skip_lipid_analysis else args.lipid_taunormalized_threshold
            )
            row["lipid_analysis_enabled"] = not args.skip_lipid_analysis
        region_rows.extend(sample_region_rows)
        lipid_rows.extend(sample_lipid_rows)

        make_qc_figure(
            triad=triad,
            tmrm_image=tmrm_images["Int"],
            gfp_image=gfp_images["Int"],
            nadph_images=nadph_images,
            valid=common_valid,
            tmrm_masks=tmrm_masks,
            gfp_groups=gfp_groups,
            lipid_mask=lipid_mask,
            lipid_droplet_mask=lipid_droplet_mask,
            final_method=args.final_tmrm_method,
            figure_path=figures / f"{triad.condition}_sample_{triad.sample:02d}_qc.png",
        )

    join_qc_df = pd.DataFrame(join_qc_rows)
    tmrm_qc_df = pd.DataFrame(tmrm_qc_rows)
    gfp_qc_df = pd.DataFrame(gfp_qc_rows)
    sample_distribution_df = pd.DataFrame(sample_distribution_rows)
    tmrm_gfp_group_df = pd.DataFrame(tmrm_gfp_group_rows)
    region_df = pd.DataFrame(region_rows)
    lipid_summary_df = pd.DataFrame(lipid_rows)
    lipid_droplet_df = pd.DataFrame(lipid_droplet_rows)
    lipid_droplet_sample_df = pd.DataFrame(lipid_droplet_sample_rows)
    contrast_df = summarize_sample_contrasts(region_df)
    tests_df = paired_tests(contrast_df)

    join_qc_df.to_csv(output / "triad_join_qc.csv", index=False)
    tmrm_qc_df.to_csv(output / "tmrm_mask_qc_all_methods.csv", index=False)
    gfp_qc_df.to_csv(output / "lgr5_gfp_mask_qc.csv", index=False)
    sample_distribution_df.to_csv(output / "sample_channel_distributions.csv", index=False)
    tmrm_gfp_group_df.to_csv(output / "tmrm_gfp_group_summary.csv", index=False)
    region_df.to_csv(output / "nadph_compartment_summary_all_methods.csv", index=False)
    lipid_summary_df.to_csv(output / "oxidized_lipid_summary_all_methods.csv", index=False)
    lipid_droplet_df.to_csv(output / "oxidized_lipid_droplet_objects.csv", index=False)
    lipid_droplet_sample_df.to_csv(output / "oxidized_lipid_droplet_sample_summary.csv", index=False)
    contrast_df.to_csv(output / "nadph_sample_contrasts_all_methods.csv", index=False)
    tests_df.to_csv(output / "paired_tests_all_methods.csv", index=False)

    final_selector = region_df["method"] == args.final_tmrm_method
    final_region_df = region_df[final_selector]
    final_lipid_df = lipid_summary_df[lipid_summary_df["method"] == args.final_tmrm_method]
    final_contrast_df = contrast_df[contrast_df["method"] == args.final_tmrm_method]
    final_tests_df = tests_df[tests_df["method"] == args.final_tmrm_method]
    final_region_df.to_csv(
        output / f"final_nadph_compartment_summary_{args.final_tmrm_method}.csv",
        index=False,
    )
    final_lipid_df.to_csv(
        output / f"final_oxidized_lipid_summary_{args.final_tmrm_method}.csv",
        index=False,
    )
    final_contrast_df.to_csv(
        output / f"final_nadph_sample_contrasts_{args.final_tmrm_method}.csv",
        index=False,
    )
    final_tests_df.to_csv(
        output / f"final_paired_tests_{args.final_tmrm_method}.csv",
        index=False,
    )

    write_by_condition_tables(
        output,
        {
            "triad_join_qc.csv": join_qc_df,
            "tmrm_mask_qc_all_methods.csv": tmrm_qc_df,
            "lgr5_gfp_mask_qc.csv": gfp_qc_df,
            "sample_channel_distributions.csv": sample_distribution_df,
            "tmrm_gfp_group_summary.csv": tmrm_gfp_group_df,
            "nadph_compartment_summary_all_methods.csv": region_df,
            "oxidized_lipid_summary_all_methods.csv": lipid_summary_df,
            "oxidized_lipid_droplet_objects.csv": lipid_droplet_df,
            "oxidized_lipid_droplet_sample_summary.csv": lipid_droplet_sample_df,
            "nadph_sample_contrasts_all_methods.csv": contrast_df,
            "paired_tests_all_methods.csv": tests_df,
        },
    )

    write_html_report(
        output=output,
        triads=triads,
        join_qc_df=join_qc_df,
        tmrm_qc_df=tmrm_qc_df,
        gfp_qc_df=gfp_qc_df,
        lipid_summary_df=lipid_summary_df,
        lipid_droplet_sample_df=lipid_droplet_sample_df,
        paired_tests_df=tests_df,
        final_method=args.final_tmrm_method,
        lipid_threshold=args.lipid_taunormalized_threshold,
        lipid_min_object_size=args.lipid_min_object_size,
        lipid_analysis_enabled=not args.skip_lipid_analysis,
    )

    print(f"Processed {len(triads)} tri-channel samples.")
    print(f"Wrote output folder: {output}")
    print(f"Wrote QC report: {output / 'qc_report.html'}")
    print(f"Final TMRM method: {args.final_tmrm_method}")
    if args.skip_lipid_analysis:
        print("Oxidized-lipid-like gate: disabled")
    else:
        print(f"Oxidized-lipid-like gate: NADPH TauNormalized > {args.lipid_taunormalized_threshold:g}")
        print(f"Lipid-like droplet object filter: area >= {args.lipid_min_object_size} pixels")


if __name__ == "__main__":
    main()
