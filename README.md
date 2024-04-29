# FUS-ALS Molecular ALS Phenotype Scores
This directory contains scripts used to process and analyze screens from the ALS
patient-derived fibroblast study. We provide a Dockerfile to install all
dependencies for running analyses reported in Kumbier et al. 2024. 

```
  docker build -t als/als .
  docker run -it --rm -v <path to project_als directory>:/ALS als/als
```

The core analysis scripts used for each screen / analysis include:

## General scripts
- `preprocessing`: directory containing scripts for processing raw image-derived
features into the selected set of ALS-relevant eigenfeatures. Proecsses image
features are available in the `data_profiles` directory of the zenodo
repository. Raw image features and images are available upon request.

- `utilities.R`: script containg functions for fitting models along with helper
functions for data processing and visualization.

- `color_palette.R`: sets color palette used in figures.

## Analysis specific scripts
Scripts used to generate figures for the paper have been copied to
`paper_figures` and are for the most part organized by figure.

- `Figure_Search`: analyses of imaging marker set search—reported in
  supplemental section. 

- `Figure_Scores`: analyses for image-based scores (i-MAP scores).

- `Figure_Transcriptomics`: analyses for transcriptomic-based scores (t-MAP
scores).

- `Figure_ASO`: analyses from ASO screen. *Note:* some of the sporadic figures
are generated in `Figure_Transcriptomics`

## Data
Data are available through the Zenodo repository:
https://zenodo.org/records/10499037. Data for image-based cell profiles are
contained in `data_profiles`, organized by screen. All other data can be found
in `data`.
