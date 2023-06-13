#' The following script generates heatmaps of prediction accuracy by genetics 
#' and marker set as reported in paper Fig. 1b. Accuracy is measured within each
#' screen individually.
#' 
#' Requires input datasets
#'   `data_profiles/110420_BiomarkerScreen/Cell_Filtered` 
#'   `data_profiles/121820_SporadicScreen/Cell_Filtered`
#' 
#'   Karl Kumbier 2/23/2022`
#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(scales)

# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA') 
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Search')

#screen <- '110420_BiomarkerScreen'
screen <- '121820_SporadicScreen'
data.dir <- file.path(data.base.dir, screen, 'Cell_Filtered')
plate.file <- file.path(fig.base.dir, 'data', 'plate_list.txt')
platemap.dir <- file.path(fig.base.dir, 'data', 'platemaps')

output.fig.dir <- file.path(figure.dir, 'fig')
dir.create(output.fig.dir, recursive=FALSE)
dir.create(platemap.dir)

# Source utility functions
source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Parameters for analysis
n.core <- 16 # number of cores for fitting models
type <- 'Cell' # sample type for modeling
genetics.train <- c('FUS-ALS', 'C9orf72', 'sporadic', 'TBK1') # genetics for training models
save.plates <- TRUE # save plate IDs to experimental summary file

# Initialize plot ordering for genetics
genetic.order <- c('Healthy', 'TBK1', 'sporadic', 'C9orf72', 'FUS-ALS')

# Initialize models to be used for supervised learning
model <- function(x, y) irf(x, y, n.core=n.core) 
model_predict <- predict_irf

#' # Loading data
#' Our dataset consists of phenotypic profiles derived from multiple marker sets. 
#' Each observation in the dataset corresponds to a single `r type`. For each 
#' marker set, we evaluate the predictive strength of a single marker by 
#' ALS genetic subtype.
#+load_data

# Load in datasets for each marker set
input.files <- list.files(data.dir) %>% str_subset('Rdata')

xlist <- lapply(input.files, function(f) {
  load(file.path(data.dir, f))
  x <- filter(x, str_detect(Compounds, '(DMSO|None|Control)')) 
  return(as.data.frame(x))
})

#' ## Predictive analysis - cell line generalization
#' To assess the generalizability of our predictive models to new cell lines, we 
#' consider a leave one cell line out strategy. Specifically, we remove a single 
#' cell line from our dataset, train models on the remaining cell lines from the
#' following genetic backgrounds `r genetics.train`, and evaluate predictions on 
#' the held-out line. The wrapper functions defined below run this analysis.
#+ leave_cell_out

#' Below, we run the analysis described above over each cell line.
#+ cell_model, cache=TRUE
set.seed(47)
fout <- file.path(output.fig.dir, str_c('model_', screen, '.Rdata'))

# Fit model to select data
fits <- lapply(xlist, function(z) {
  print('-----------------------------------------------')
  print(z$Markers[1])
  cell_fit(z, model, model_predict, genetics.train)$ypred
})

fits <- fits[!sapply(fits, is.null)]
save(file=fout, fits)


# Save table of plates used
if (save.plates) {
  plate.list <- lapply(xlist, function(z) unique(z$PlateID)) %>% unlist
  
  write.table(
    file=plate.file, 
    unique(plate.list), 
    append=TRUE, 
    row.names=FALSE, 
    col.names=FALSE, 
    quote=FALSE
  )
  
  # Save plate maps
  plate.layout <- lapply(xlist, function(z) {
    dplyr::select(z, PlateID, CellLine, RowID, ColID) %>% distinct()
  })

  plate.layout <- rbindlist(plate.layout)  
  
  for (p in unique(plate.layout$PlateID)) {
    platemap <- matrix('', nrow=16, ncol=24)
    pm <- filter(plate.layout, PlateID == p)
    for (i in 1:nrow(pm)) platemap[pm$RowID[i], pm$ColID[i]] <- pm$CellLine[i]
    colnames(platemap) <- 1:24
    
    print(table(platemap))
    fout <- file.path(platemap.dir, str_c(p, '.csv'))
    write.csv(file=fout, platemap)
  }
}
