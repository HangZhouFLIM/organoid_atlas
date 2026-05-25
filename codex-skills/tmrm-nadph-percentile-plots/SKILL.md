---
name: tmrm-nadph-percentile-plots
description: Generate reusable publication-style percentile-p75 TMRM/NAD(P)H plots and statistics from combined analysis CSV folders. Use when the user asks to reuse the TMRM-high versus low-TMRM NAD(P)H plotting workflow, make delta plots, paired high/low plots, GraphPad-like Nature-style figures, or run the existing NADPH/TMRM plotting code on another folder containing combined_nadph_sample_differences_all_methods.csv and combined_nadph_region_summary_all_methods.csv.
---

# TMRM NAD(P)H Percentile Plots

## Workflow

Use the bundled script for repeatable plotting instead of rewriting the analysis code.

1. Locate the analysis folder. It should either be the folder containing:
   - `combined_nadph_sample_differences_all_methods.csv`
   - `combined_nadph_region_summary_all_methods.csv`

   or a project root containing `analysis_TMRM_NADPH_all_conditions/` with those files.

2. Run:

```powershell
& "<python.exe>" "<skill>/scripts/plot_nadph_percentile_p75.py" --analysis-dir "<analysis-or-project-folder>"
```

Use the existing project Python environment when available. If Python packages are missing, install or request permission for `pandas`, `numpy`, `scipy`, and `matplotlib`.

3. Check outputs in the script-reported output directory. By default it writes:

```text
figures_percentile_p75/
  nadph_delta_percentile_p75.svg/pdf/png/tiff
  nadph_high_low_paired_percentile_p75.svg/pdf/png/tiff
  source_data_nadph_percentile_p75.csv
  stats_nadph_percentile_p75.csv
```

4. To generate the style gallery variants, run:

```powershell
& "<python.exe>" "<skill>/scripts/plot_nadph_percentile_p75_style_gallery.py" --analysis-dir "<analysis-or-project-folder>"
```

By default the gallery writes:

```text
figures_percentile_p75/style_gallery/
  graphpad_prism_delta_percentile_p75.svg/pdf/png
  graphpad_prism_paired_percentile_p75.svg/pdf/png
  nature_compact_delta_percentile_p75.svg/pdf/png
  nature_compact_paired_percentile_p75.svg/pdf/png
  nature_minimal_delta_percentile_p75.svg/pdf/png
  nature_minimal_paired_percentile_p75.svg/pdf/png
  print_safe_gray_delta_percentile_p75.svg/pdf/png
  print_safe_gray_paired_percentile_p75.svg/pdf/png
  style_gallery_index.html
  style_gallery_manifest.csv
```

## Defaults

- Method: `percentile_p75`
- Condition order: `Control`, `Rotenone`, `Antimycin A`, `Oligomycin`
- Statistical unit: sample/image, not pixels
- Main delta: `high-TMRM region - low-TMRM region`
- Delta y labels:
  - `Δ NAD(P)H photon/pixel`
  - `Δ NAD(P)H τ(φ) ns`
  - `Δ NAD(P)H τ(m) ns`
  - `Δ NAD(P)H τ(N) ns`
- Region labels:
  - `Low-TMRM intracellular`
  - `TMRM-high mitochondrial`

## Statistics

Use the same decisions as the original workflow:

- Shapiro-Wilk normality per condition and metric.
- Within each condition, test delta versus zero:
  - normal: one-sample two-tailed t-test
  - non-normal: two-tailed Wilcoxon signed-rank test
- Between conditions:
  - all groups normal: Brown-Forsythe variance test, Welch ANOVA, then Dunnett T3 treatment-vs-control comparisons
  - any group non-normal: Kruskal-Wallis, then Dunn treatment-vs-control comparisons with Holm adjustment
- Annotate delta plots with treatment-vs-control adjusted stars only; keep exact statistics in CSV.

## Common Options

Use these when a new dataset has different folders, labels, or methods:

```powershell
& "<python.exe>" "<skill>/scripts/plot_nadph_percentile_p75.py" `
  --analysis-dir "<folder>" `
  --method percentile_p75 `
  --output-dir "<folder>\figures_percentile_p75" `
  --condition-order control rotenone AA oligomycin `
  --condition-label AA="Antimycin A" `
  --condition-color AA="#C95A71"
```

Only change labels/colors/order if the source CSV uses different condition keys.
