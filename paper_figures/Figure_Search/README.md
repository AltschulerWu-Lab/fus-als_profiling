# Marker search screens
Scripts in this directory are used to train and evalaute models used in marker
set search screens and generate figures summarizing model performance. Data are
derived from three screens:

- `data_profiles/110420_BiomarkerScreen/Cell_Filtered`
- `data_profiles/121820_SporadicScreen/Cell_Filtered`
- `data_profiles/030321_FUSWT/Cell_Filtered`

## Model fitting
Models for the first two screens can be fit by running the
`train_marker_screens.R` script. Select the screen to analyze by setting the
`screen` variable as one of '110420_BiomarkerScreen' or '121820_SporadicScreen'.
Models for the third screen can be trained using the `train_fus.R` script. Each
of these scripts saves the output model in `output.fig.dir/model_<screen>.Rdata`

### Parameters in modeling scripts

- `genetics.train`: cell line genetics to be used for model training. Default is
  to use all genetics available for a given marker set. (only defined for
  `train_marker_screens.R`)

- `model`: modeling function that takes arguments `x` (feature matrix) and `y`
  (response vector) and returns a fitted model.

- `model_predict: model prediction function that takes arguments `model` (a
  fitted mode) and `x` (feature matrix) and returns predictions on `x`.

- `markers`: marker sets to be used for model. Single marker models can be
  defined from a marker pair as <marker1>_<marker2>:<marker*>. (only defined for
  `train_fus_screen.R`)

## Figures
Scripts to generate figures for the marker search screens require fitted models,
generated using the scripts described above. `fig_heatmap.R` produces a heatmap
for a selected `screen` reporting classification accuracy by genetics / marker.
`fig_marker_tuning.R` produces figures summarizing marker set vs. individual
marker performance for the FUS screen.
