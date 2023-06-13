#' The following script fits models by marker set for the FUS, CD63, EEA1 
#' comparison screen. We consider model performance with respect to both marker 
#' pairs and individual markers. 
#' 
#' Requires input datasets
#'   `data_profiles/030321_FUSWT_Cell`
#'   
#'   Karl Kumbier 2/23/2022
#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(scales)

# Set paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Search')

data.dir <- file.path(data.base.dir, '030321_FUSWT', 'Cell_Filtered')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')
plate.file <- file.path(fig.base.dir, 'data', 'plate_list.txt')
platemap.dir <- file.path(fig.base.dir, 'data', 'platemaps')

output.fig.dir <- file.path(figure.dir, 'fig')
dir.create(output.fig.dir, recursive=FALSE)
dir.create(platemap.dir)

# Source utility funcitons
source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))

# Parameters for analysis
type <- 'Cell'
n.core <- 16

markers <- c(
  'FUS_EEA1', 'CD63_EEA1', 'CD63_FUS',
  'FUS_EEA1:FUS', 'CD63_EEA1:CD63', 'CD63_FUS:CD63',
  'FUS_EEA1:EEA1', 'CD63_EEA1:EEA1', 'CD63_FUS:FUS'
)

save.plates <- FALSE

# Initialize models to be used for supervised learning
model <- function(x, y) irf(x, y, n.core=n.core) 
model_predict <- predict_irf

#' # Loading data
#' Our dataset consists of phenotypic profiles derived from `r marker`. Each 
#' observation in the dataset corresponds to a single `r type`.
#+load_data
# Load in datasets for selected marker
load_marker <- function(marker) {
  # wrapper function for loading data from select marker set
  load(file.path(data.dir, str_c('profiles_', marker, '.Rdata')))
  return(x)
}

#' ## Predictive analysis
#' Here we address the question of which marker is best for separating WT from 
#' FUS cell lines. Toward this end, we train classifiers to discriminate between 
#' disease and healthy `r type`s. Classifiers are trained using features 
#' derived from select markers / marker sets.
#' 
#' To assess the generalizability of our predictive models to new cell lines, we 
#' consider a leave one cell line out strategy. Specifically, we remove a single 
#' cell line, train models on the remaining cell lines, and evaluate predictions 
#' on the held-out line. 
#+ leave_cell_out
marker_fit <- function(marker) {
  # Wrapper function to call cell_fit over each marker set
  
  marker.split <- str_split(marker, ':')[[1]]
  
  # Fit model to predict each cell lines
  x <- load_marker(marker.split[1])

  # Filter to singletons
  if (length(marker.split) == 2) {
    marker.set <- str_split(marker.split[1], '_')[[1]]
    single.marker <- marker.split[2]
    channel <- which(marker.set == single.marker) + 1
    
    xchannel <- sapply(str_split(colnames(x), '\\.\\.'), '[', 3)
    id.meta <- !str_detect(colnames(x), '^X')
    id.marker <- xchannel %in% c(channel, str_c(channel, channel))
    x <- select_if(x, id.marker | id.meta)
  }

  out <- cell_fit(x, model, model_predict)
  out$ypred$Markers <- marker
  return(out)
}

# Fit model to select data
set.seed(47)
fout <- file.path(output.fig.dir, 'model_030321_FUSWT.Rdata')

fits <- lapply(markers, marker_fit)
fits <- lapply(fits, function(f) f$ypred)
save(file=fout, fits)


################################################################################
# Save table of plates used
################################################################################
if (save.plates) {
  predictions <- rbindlist(fits)
  plate.list <- unique(predictions$PlateID)
  
  write.table(
    file=plate.file, 
    plate.list, 
    append=TRUE, 
    row.names=FALSE, 
    col.names=FALSE, 
    quote=FALSE
  )
  
  # Save plate maps
  plate.layout <- dplyr::select(predictions, PlateID, CellLine, WellID) %>%
    distinct() %>%
    mutate(Row=as.numeric(str_remove_all(WellID, '-.*$'))) %>%
    mutate(Col=as.numeric(str_remove_all(WellID, '^.*-')))
  
  for (p in unique(plate.layout$PlateID)) {
    platemap <- matrix('', nrow=16, ncol=24)
    pm <- filter(plate.layout, PlateID == p)
    for (i in 1:nrow(pm)) platemap[pm$Row[i], pm$Col[i]] <- pm$CellLine[i]
    colnames(platemap) <- 1:24
    
    print(table(platemap))
    fout <- file.path(platemap.dir, str_c(p, '.csv'))
    write.csv(file=fout, platemap)
  }
}
