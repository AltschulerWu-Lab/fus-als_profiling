###############################################################################
#+ setup, message=FALSE, echo=FALSE
###############################################################################
library(tidyverse)
library(data.table)

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Search")

#screen <- "121820_SporadicScreen"
args <- R.utils::commandArgs(asValues=TRUE)
screen <- args$SCREEN
if (is.null(screen)) screen <- "110420_BiomarkerScreen" 

data.dir <- file.path(data.base.dir, screen, "Cell_Filtered")
plate.file <- file.path(fig.base.dir, "data", "plate_list.txt")
platemap.dir <- file.path(fig.base.dir, "data", "platemaps")

output.fig.dir <- file.path(figure.dir, "fig")
dir.create(output.fig.dir, recursive = FALSE)
dir.create(platemap.dir)

# Source utility functions
source(file.path(fig.base.dir, "paper_figures", "utilities.R"))

# save plate IDs to experiment summary file
save.plates <- TRUE

# number of cores to use
n.core <- 16

###############################################################################
#+ parameters
###############################################################################
# genetics used for model training
genetics.train <- c("FUS-ALS")

# models used for supervised learning
model <- function(x, y) irf(x, y, n.core = n.core)
model_predict <- predict_irf

###############################################################################
#' # Load data
#+load_data
###############################################################################
# Load in datasets for each marker set
input.files <- list.files(data.dir) %>% str_subset("Rdata")

xlist <- lapply(input.files, function(f) {
  load(file.path(data.dir, f))
  x <- filter(x, str_detect(Compounds, "(DMSO|None|Control)"))
  x <- filter(x, Genetics %in% c("FUS-ALS", "Healthy"))
  return(as.data.frame(x))
})

###############################################################################
#' # Predictive analysis
#+ modeling, cache=TRUE
###############################################################################
set.seed(47)
fout <- file.path(output.fig.dir, str_c("model_", screen, ".Rdata"))

# Fit model to select data
fits <- lapply(xlist, function(z) {
  print("-----------------------------------------------")
  print(z$Markers[1])
  cell_fit(z, model, model_predict, genetics.train)$ypred
})

fits <- fits[!sapply(fits, is.null)]
save(file = fout, fits)


# Save table of plates used
if (save.plates) {
  plate.list <- lapply(xlist, function(z) unique(z$PlateID)) %>% unlist()

  write.table(
    file = plate.file,
    unique(plate.list),
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # Save plate maps
  plate.layout <- lapply(xlist, function(z) {
    dplyr::select(z, PlateID, CellLine, RowID, ColID) %>% distinct()
  })

  plate.layout <- rbindlist(plate.layout)

  for (p in unique(plate.layout$PlateID)) {
    platemap <- matrix("", nrow = 16, ncol = 24)
    pm <- filter(plate.layout, PlateID == p)
    for (i in 1:nrow(pm)) platemap[pm$RowID[i], pm$ColID[i]] <- pm$CellLine[i]
    colnames(platemap) <- 1:24

    print(table(platemap))
    fout <- file.path(platemap.dir, str_c(p, ".csv"))
    write.csv(file = fout, platemap)
  }
}
