# organoid_atlas

This repository contains a MATLAB pipeline for analysing crypt-to-villus
cell gradients. The combined script `crypt_to_villus_gradients_combined.m`
allows interactive processing of CSV files with cell measurements and
produces gradient summaries for the metrics **MeanInt**, **MeanTauPhase**
and **MeanTauModulation**.

Run the script in MATLAB and follow on-screen prompts to select your files
and define the crypt and villus reference cells. Results and figures are
saved under `./output/<Metric>/` for each metric, including a pooled mean ± SEM
curve and a combined single-cell heatmap.

## Binary mask utility for napari

`batch_make_binary_from_masks.py` provides a small helper for generating
binary masks from grayscale mask images within a napari session. Launch napari
with an existing `viewer` instance and run the script in the console:

```python
from batch_make_binary_from_masks import batch_make_binary_from_masks

batch_make_binary_from_masks()
```

Select one or more mask files when prompted. Each mask is binarized (non-zero
pixels become white), added to the viewer, and saved alongside the source
image with a `_binaryMask` suffix.
