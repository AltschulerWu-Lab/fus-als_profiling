###############################################################################
#+ setup, message=FALSE
###############################################################################
library(tidyverse)
library(data.table)

# Set paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Search")

data.dir <- file.path(data.base.dir, "030321_FUSWT", "Cell_Filtered")
library.file <- file.path(fig.base.dir, "data", "Fibroblasts_Library.csv")
plate.file <- file.path(fig.base.dir, "data", "plate_list.txt")
platemap.dir <- file.path(fig.base.dir, "data", "platemaps")

output.fig.dir <- file.path(figure.dir, "fig")
dir.create(output.fig.dir, recursive = FALSE)
dir.create(platemap.dir)

# Source utility funcitons
source(file.path(fig.base.dir, "paper_figures", "utilities.R"))

# save plate IDs to experiment summary file
save.plates <- FALSE

# Parameters for analysis
n.core <- 16

###############################################################################
#+ parameters
###############################################################################
markers <- c(
  "FUS_EEA1", "CD63_EEA1", "CD63_FUS",
  "FUS_EEA1:FUS", "CD63_EEA1:CD63", "CD63_FUS:CD63",
  "FUS_EEA1:EEA1", "CD63_EEA1:EEA1", "CD63_FUS:FUS"
)

# Initialize models to be used for supervised learning
model <- function(x, y) irf(x, y, n.core = n.core)
model_predict <- predict_irf

###############################################################################
#' # Load data
#+load_data
###############################################################################
# Load in datasets for selected marker
load_marker <- function(marker) {
  # wrapper function for loading data from select marker set
  load(file.path(data.dir, str_c("profiles_", marker, ".Rdata")))
  return(x)
}

################################################################################
#' # Predictive analysis
#+ modeling
################################################################################
marker_fit <- function(marker) {
  # Wrapper function to call cell_fit over each marker set

  marker.split <- str_split(marker, ":")[[1]]

  # Fit model to predict each cell lines
  x <- load_marker(marker.split[1])

  # Filter to singletons
  if (length(marker.split) == 2) {
    marker.set <- str_split(marker.split[1], "_")[[1]]
    single.marker <- marker.split[2]
    channel <- which(marker.set == single.marker) + 1

    xchannel <- sapply(str_split(colnames(x), "\\.\\."), "[", 3)
    id.meta <- !str_detect(colnames(x), "^X")
    id.marker <- xchannel %in% c(channel, str_c(channel, channel))
    x <- select_if(x, id.marker | id.meta)
  }

  out <- cell_fit(x, model, model_predict)
  out$ypred$Markers <- marker
  return(out)
}

# Fit model to select data
set.seed(47)
fout <- file.path(output.fig.dir, "model_030321_FUSWT.Rdata")

fits <- lapply(markers, marker_fit)
fits <- lapply(fits, function(f) f$ypred)
save(file = fout, fits)


################################################################################
# Save table of plates used
################################################################################
if (save.plates) {
  predictions <- rbindlist(fits)
  plate.list <- unique(predictions$PlateID)

  write.table(
    file = plate.file,
    plate.list,
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # Save plate maps
  plate.layout <- dplyr::select(predictions, PlateID, CellLine, WellID) %>%
    distinct() %>%
    mutate(Row = as.numeric(str_remove_all(WellID, "-.*$"))) %>%
    mutate(Col = as.numeric(str_remove_all(WellID, "^.*-")))

  for (p in unique(plate.layout$PlateID)) {
    platemap <- matrix("", nrow = 16, ncol = 24)
    pm <- filter(plate.layout, PlateID == p)
    for (i in 1:nrow(pm)) platemap[pm$Row[i], pm$Col[i]] <- pm$CellLine[i]
    colnames(platemap) <- 1:24

    print(table(platemap))
    fout <- file.path(platemap.dir, str_c(p, ".csv"))
    write.csv(file = fout, platemap)
  }
}
