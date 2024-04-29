# FUS-ALS imaging scores
Scripts in this directory train and evaluate imaging-based FUS-ALS vs. Healthy
classifiers and generate figures summarizing model performance. Data are derived
from:

- `data_profiles/032422_WTFUSSporadics/Cell_Filtered`
- `data_profiles/032422_WTFUSSporadics/Well`

## Modeling
Models can be fit using the `train_fus_models.R` script.

### Parameters in modeling scripts
- `cytosolic`: TRUE/FALSE, indicating whether models should be restricted to
  cytosolic FUS features.

- `model`: modeling function that takes arguments `x` (feature matrix) and `y`
  (response vector) and returns a fitted model.

- `model_predict: model prediction function that takes arguments `model` (a
  fitted mode) and `x` (feature matrix) and returns predictions on `x`.

## Figures
Scripts to generate figures for the marker search screens require fitted models,
generated using the script described above. `fig_scores.R` generates summary
plots based on the FUS-ALS vs. WT classifier that appear in the main text of the
manuscript. To generate summaries for different models (e.g., irf, logistic
regression, neural network), change the `model.str` argument in this script.
`fig_coordinate_space.R` generates plots of the cell lines in 2D coordinate
space based on raw features. `fig_qc_pca.R` plots points in PCA space.
`fig_qc_plates.R` runs models and evaluates performance on hold-out plates.
`fig_model_comparisons.R` generates plots comparing predictions across logistic
regression, irf, and neural networks. `run_models.sh` provides a shell wrapper
to iterate over models and cytosolic / full cell analyses.
