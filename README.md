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
