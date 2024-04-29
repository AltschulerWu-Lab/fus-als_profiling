# FUS-ALS transcriptomic scores
Scripts in this directory train and evaluate transcriptomic-based FUS-ALS vs.
Healthy classifiers and generate figures summarizing model performance. Data for
both spinal cord and fibroblast classifiers are derived from:

- **fibroblast raw read counts:** `data/050322_RNAseq/read_counts.Rdata`
- **spinal cord raw read counts:** `data/050322_RNAseq/nygc/read_counts.Rdata`

## Preprocessing
Spinal cord and fibroblast RNA-seq datasets are preprocessed from raw read
counts using the `load_spine.R` and `load_fibro.R` scripts. Running these
scripts generates the `processed_*.Rdata` files contained in 
`data/050322_RNAseq`. The resulting files contain sparse matrices with
both (i) raw read counts and (ii) TPM for genes that are both annotated in 
the reactome pathway database and expressed in fibroblasts / spinal cord 
samples.

## Modeling
Models can be fit using the `fit_*.R` scripts. The spinal cord models train ALS
vs. healthy control classifiers on patient-derived spinal cord samples, using
subsets of genes corresponding to reactome pathways. The fibroblast models train
classifiers using the subset of reactome pathways with better than random
prediction in the spinal cord data.

### Parameters in modeling scripts
#### Spinal cord model
- `min.healthy`: minimum number of healthy controls for a site to be included in
  analysis. I.e., data from sites with fewer than `min.healthy` HCs are dropped.

- `n.null`: number of bootstrap samples to use in generating null distribution
  on metadata features.

- `min.size`, `max.size`: minimum and maximum pathway sizes included in
  analysis.

#### Fibroblast model
- `n.rep`: number of training model replicates to fit for each hold-out cell
  line. Upsampling for class balance, and optional bootstrap samples, result in
  variation across model replicates.
- `bootstrap`: T/F, should models be trained on bootstrap samples of the
  training data. 
- `gene_transform`: transformation applied to each feature. Default is rank
  transform.


## Figures
`fig_*.R` scripts generate figures reported in the manuscript. Models must be
trained first, using the `fit_*.R` scripts described above.
