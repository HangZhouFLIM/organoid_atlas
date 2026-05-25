---
name: tmrm-nadph-analysis
description: Analyze paired TMRM and NADPH FLIM pxData CSV exports from folders containing TMRM/ and NADPH/ subfolders. Use when the user asks to run TMRM/NADPH analysis, create high-TMRM mitochondrial region QC masks, compare high- versus low-TMRM NADPH intensity or lifetime metrics, choose percentile_p75/otsu_clean/adaptive_clean masks, or generate final per-sample and paired-test CSV outputs from raw FLIM exports.
---

# TMRM/NADPH FLIM Analysis

Use this skill for raw paired TMRM/NADPH FLIM export folders with this shape:

```text
experiment_or_condition/
  TMRM/
    *_pxData.csv
  NADPH/
    *_pxData.csv
```

The bundled script expects pixel-level CSVs with `Int`, `TauPhase`, `TauModulation`, `TauNormalized`, `Ycoord`, and `Xcoord`. It pairs files by sample number in filenames shaped like `prefix_sample_acquisition_pxData.csv`, for example `Oligomycin_1_01_pxData.csv` paired with `Oligomycin_1_02_pxData.csv`.

Use `tmrm-nadph-percentile-plots` instead when the inputs are already combined analysis CSVs and the user only wants the reusable percentile/TMRM-high plotting workflow.

## Preferred Runtime

Use the shared environment first:

```powershell
$PY = Join-Path $env:USERPROFILE "analysis_envs\tmrm_nadph\.venv\Scripts\python.exe"
```

If it is missing, create it with the user's approval:

```powershell
$VENV = Join-Path $env:USERPROFILE "analysis_envs\tmrm_nadph\.venv"
python -m venv $VENV
$PY = Join-Path $VENV "Scripts\python.exe"
$REQUIREMENTS = Join-Path $env:USERPROFILE ".codex\skills\tmrm-nadph-analysis\scripts\requirements.txt"
& $PY -m pip install -r $REQUIREMENTS
```

## Workflow

1. Confirm the target folder contains `TMRM/` and `NADPH/` subfolders and inspect `*_pxData.csv` names.
2. Ensure the runtime has the packages in `scripts/requirements.txt`.
3. Run QC first:

```powershell
$SCRIPT = Join-Path $env:USERPROFILE ".codex\skills\tmrm-nadph-analysis\scripts\analyze_tmrm_nadph.py"
& $PY $SCRIPT --root "PATH_TO_CONDITION_FOLDER" --output "PATH_TO_CONDITION_FOLDER\analysis_TMRM_NADPH"
```

4. Inspect `qc_report.html` and the per-sample images in `figures/`. If visual QC is ambiguous, ask the user to choose one final mask method: `percentile_p75`, `otsu_clean`, or `adaptive_clean`.
5. Run final chosen-mask statistics:

```powershell
& $PY $SCRIPT --root "PATH_TO_CONDITION_FOLDER" --output "PATH_TO_CONDITION_FOLDER\analysis_TMRM_NADPH" --final-method otsu_clean
```

## Analysis Defaults

- High-TMRM candidate masks: `percentile_p75`, `otsu_clean`, and `adaptive_clean`.
- Low-TMRM region: exported/valid pixels outside the selected high-TMRM mask.
- NADPH metrics: `Int`, `TauPhase`, `TauModulation`, and `TauNormalized`.
- Treat samples as independent replicates; do not treat pixels as independent samples for final statistics.
- Default image reconstruction shape is `246 x 246`; adjust `IMAGE_SHAPE` in the script only if coordinates exceed that shape.
- The script requires row-aligned TMRM and NADPH coordinates for each paired sample.

## Outputs

The script writes `qc_report.html`, per-sample mask figures, `mask_qc_summary.csv`, `nadph_region_summary_all_methods.csv`, `nadph_sample_differences_all_methods.csv`, and `paired_tests_all_methods.csv`. When `--final-method` is supplied, it also writes `final_region_summary_METHOD.csv`, `final_sample_differences_METHOD.csv`, and `final_paired_tests_METHOD.csv`.
