This directory contains scripts used to process and analyze screens from the ALS
patient-derived fibroblast study. Analyses for each screen are organized within
individual directories corresponding to Figures 2-4 of the manuscript
Figure_Scores=Fig. 2, Figure_ASO=Fig. 3, Figure_Transcriptomics=Fig. 4,
Figure_Search=Extended data Figs. E1-E3. 

Data required to run the scripts can be downloaded from Zenodo:
https://zenodo.org/record/7247995. To run scripts, Initialize `.Renviron`
variables:

- `ALS_PAPER`=`<path_to_this_directory>`

- `ALS_DATA`=`<path_to_data_profiles_directory>` 

The `data` directory on Zenodo must be placed in this directory.
