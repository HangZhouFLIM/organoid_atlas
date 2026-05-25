---
name: tmrm-lgr5-nadph-analysis
description: Run QC-first tri-channel LGR5-GFP/TMRM/NADPH FLIM pxData analysis for folders containing lgr5_gfp, TMRM, and NADPH channels. Use when the user asks to analyze stem-like LGR5-GFP low/medium/high or high/rest compartments, TMRM-high mitochondrial allocation, real NADPH shifts, optional oxidized-lipid-like droplet objects, QC reports, and publication-style plots across conditions.
---

# TMRM/LGR5-GFP/NADPH FLIM Analysis

Use this skill for tri-channel FLIM export folders shaped like:

```text
project_root/
  condition_a/
    lgr5_gfp/*_pxData.csv
    TMRM/*_pxData.csv
    NADPH/*_pxData.csv
  condition_b/
    lgr5_gfp/*_pxData.csv
    TMRM/*_pxData.csv
    NADPH/*_pxData.csv
```

The scripts align pixels by `Xcoord,Ycoord`, not row order. Treat samples/images as replicates; never treat pixels or droplets as biological replicates for final tests.

## Required Gates

1. **Ask before lipid-droplet analysis.** Before final analysis, ask whether this NADPH dataset needs oxidized-lipid-like separation/droplet analysis. Do not assume yes. If yes, use `TauNormalized > 4.0` by default and count connected droplet objects with minimum area `>= 5 px`.
2. **QC before choosing GFP grouping.** Inspect GFP intensity ranges, histograms, and QC masks before final plots.
   - Use `p25/p75` three groups when low, medium, and high GFP are interpretable: `low <= p25`, `medium p25-p75`, cleaned `high >= p75`.
   - Use `p75-only` high/rest interpretation when GFP is weak or the low tail is not biologically meaningful. In that case, treat cleaned `>=p75` as high/stem-like and all valid non-high pixels as rest/low.
   - If ambiguous, show the QC report and ask the user to choose before plotting.

## Runtime

Prefer the shared environment:

```powershell
$PY = Join-Path $env:USERPROFILE "analysis_envs\tmrm_nadph\.venv\Scripts\python.exe"
```

If packages are missing, install/request: `pandas`, `numpy`, `scipy`, `matplotlib`, `scikit-image`.

## Workflow

1. Confirm the folder contains condition folders with `lgr5_gfp`, `TMRM`, and `NADPH`.
2. Inspect one CSV header per channel. TMRM and NADPH should ideally include `Int`, `TauPhase`, `TauModulation`, `TauNormalized`, `Ycoord`, `Xcoord`.
3. Ask the lipid-droplet gate question.
4. Run QC/analysis with the bundled script:

```powershell
$SCRIPT = Join-Path $env:USERPROFILE ".codex\skills\tmrm-lgr5-nadph-analysis\scripts\analyze_tmrm_lgr5_nadph.py"
& $PY $SCRIPT --root "PROJECT_ROOT" --output "PROJECT_ROOT\analysis_TMRM_LGR5_NADPH" --final-tmrm-method percentile_p75 --gfp-grouping p25_p75
```

If the user says the dataset does **not** need oxidized-lipid-like separation, add:

```powershell
--skip-lipid-analysis
```

If QC shows p75-only high/rest is better than low/medium/high, rerun final analysis with:

```powershell
--gfp-grouping p75_only
```

5. Open `analysis_TMRM_LGR5_NADPH/qc_report.html`. Check:
   - GFP raw intensity and histogram p25/p75
   - low/medium/cleaned-high masks
   - whether cleaned high GFP follows cell-like regions
   - TMRM-high mask quality
   - NADPH real-signal exclusion and optional lipid-like droplets
6. Decide final GFP interpretation (`p25/p75` three groups or `p75-only` high/rest). If the bundled script needs adjustment for strict p75-only plots, adapt the plotting-level grouping before final export and state the assumption.
7. Generate plots:

```powershell
$PLOT = Join-Path $env:USERPROFILE ".codex\skills\tmrm-lgr5-nadph-analysis\scripts\plot_lgr5_stem_like_percentile_p75.py"
& $PY $PLOT --analysis-dir "PROJECT_ROOT\analysis_TMRM_LGR5_NADPH"
```

## Expected Outputs

Analysis folder:

- `qc_report.html`
- `lgr5_gfp_mask_qc.csv`
- `tmrm_gfp_group_summary.csv`
- `final_nadph_compartment_summary_percentile_p75.csv`
- `final_nadph_sample_contrasts_percentile_p75.csv`
- optional `oxidized_lipid_droplet_objects.csv`
- optional `oxidized_lipid_droplet_sample_summary.csv`

Figure folder:

- `tmrm_features_stem_like_percentile_p75`
- `real_nadph_delta_stem_like_percentile_p75`
- `real_nadph_high_low_paired_stem_like_percentile_p75`
- support `real_nadph_delta_three_gfp_groups_percentile_p75`
- optional oxidized-lipid droplet density and droplet metric figures

Export SVG/PDF/PNG/TIFF plus source, summary, and stats CSVs.

## Validation Checklist

- All expected samples are included.
- GFP groups sum to valid common pixels.
- Stem-like pixels equal medium + cleaned high for three-group runs.
- Real NADPH excludes oxidized-lipid-like pixels when the lipid gate is enabled.
- Optional lipid droplets are condition-only, not GFP-grouped.
- Plots use sample-level dots/means and match the percentile-p75 plotting style.
