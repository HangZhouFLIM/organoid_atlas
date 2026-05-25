---
name: nadph-phasor-trajectory
description: Generate publication-style NADPH phasor G/S shift trajectory figures from paired TMRM and NADPH FLIM pxData CSV folders. Use when the user asks to plot NADPH phasor shifts, G/S trajectories, control-to-treatment phasor arrows, or TMRM-high versus low-TMRM NADPH phasor centroids across conditions such as control, rotenone, antimycin A, and oligomycin.
---

# NADPH Phasor Trajectory

Use this skill to draw a two-panel phasor trajectory figure from paired TMRM/NADPH FLIM CSV exports. The script uses raw NADPH `G` and `S` columns directly and defines regions with the same `percentile_p75` TMRM-high logic used in the TMRM/NADPH analysis workflow.

## Input Shape

Preferred project shape:

```text
project_root/
  analysis_TMRM_NADPH_aligned_inputs/
    control/
      TMRM/*_pxData.csv
      NADPH/*_pxData.csv
    rotenone/
      TMRM/*_pxData.csv
      NADPH/*_pxData.csv
```

If `analysis_TMRM_NADPH_aligned_inputs/` is absent, the script reads raw condition folders directly under `project_root/` and aligns TMRM/NADPH rows in memory by shared `Xcoord,Ycoord`.

Required NADPH columns: `G`, `S`, `Int`, `Ycoord`, `Xcoord`. Required TMRM columns: `Int`, `Ycoord`, `Xcoord`.

## Workflow

1. Use Python/matplotlib only for plotting.
2. Prefer the shared runtime when available:

```powershell
$PY = Join-Path $env:USERPROFILE "analysis_envs\tmrm_nadph\.venv\Scripts\python.exe"
```

3. Run the bundled script:

```powershell
$SCRIPT = Join-Path $env:USERPROFILE ".codex\skills\nadph-phasor-trajectory\scripts\plot_nadph_phasor_trajectory.py"
& $PY $SCRIPT --root "PATH_TO_PROJECT_ROOT"
```

4. If the dataset uses nonstandard condition names or order, pass control first:

```powershell
& $PY $SCRIPT --root "PATH_TO_PROJECT_ROOT" --condition-order control rotenone AA oligomycin
```

## Outputs

Default output folder:

```text
project_root/analysis_TMRM_NADPH_all_conditions/figures_phasor_percentile_p75/
```

The script writes:

```text
nadph_phasor_trajectory_percentile_p75.svg
nadph_phasor_trajectory_percentile_p75.pdf
nadph_phasor_trajectory_percentile_p75.png
nadph_phasor_trajectory_percentile_p75.tiff
source_data_nadph_phasor_trajectory_percentile_p75.csv
group_summary_nadph_phasor_trajectory_percentile_p75.csv
trajectory_vectors_nadph_phasor_trajectory_percentile_p75.csv
```

## Figure Contract

- Statistical unit: sample/image centroid, not pixels.
- Regions: `TMRM-high mitochondrial region` and `Low-TMRM intracellular region`.
- Trajectory: radial arrows from control mean to each treatment mean within each region.
- Uncertainty: sample centroids plus SEM crosshairs.
- Reference: faint universal phasor semicircle; do not add free/bound NADH labels unless trusted reference coordinates are provided.

## Validation

After running, check the script log:

- source rows should equal `number_of_samples * 2 regions`.
- group centroids should equal `number_of_conditions * 2 regions`.
- trajectory arrows should equal `(number_of_conditions - 1) * 2 regions`.

Preview the PNG or SVG to confirm the universal semicircle, points, SEM crosshairs, arrows, and labels are readable.
