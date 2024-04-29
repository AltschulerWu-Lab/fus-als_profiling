################################################################################
#+ setup, message=FALSE
################################################################################
library(tidyverse)
library(data.table)

# Set user args for model choice
args <- R.utils::commandArgs(asValues = TRUE)
cytosolic <- as.logical(args$CYTOSOLIC) # TRUE/FALSE
model.str <- args$MODEL # 'irf', 'logistic', or 'fcnn'

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Scores")
platemap.dir <- file.path(fig.base.dir, "data", "platemaps")
dir.create(platemap.dir)

source(file.path(fig.base.dir, "paper_figures", "utilities.R"))

# Set paths to data directories
data.dir <- file.path(data.base.dir, "032422_WTFUSSporadics", "Cell_Filtered")
plate.file <- file.path(fig.base.dir, "data", "plate_list.txt")

n.core <- 10
marker <- "FUS_EEA1"
save.fig <- FALSE
save.plate <- FALSE

################################################################################
#+ parameters, message=FALSE
################################################################################

# Initialize models to be used for supervised learning
if (model.str == "fcnn") {
  model <- fcnn
  model_predict <- predict_fcnn
} else if (model.str == "logistic") {
  model <- logistic
  model_predict <- predict_logistic
} else {
  model <- function(x, y) irf(x, y, n.core = n.core)
  model_predict <- predict_irf
}

ext <- ifelse(cytosolic, "_cytosolic", "_full")
output.fig.dir <- file.path(figure.dir, model.str, str_c("fig", ext))
dir.create(output.fig.dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
#' # Load data
#+ parameters, message=FALSE
################################################################################
# Load in datasets for selected marker
input.file <- str_c("profiles_", marker, ".Rdata")
load(file.path(data.dir, input.file))

# Filter to cytosolic features
if (cytosolic) x <- select_cytosolic(x)

################################################################################
#' # Predictive analysis
#+ modeling, message=FALSE
################################################################################
# Fit model to select data
set.seed(47)

if (!file.exists(file.path(output.fig.dir, "model.Rdata"))) {
  cell.models <- cell_fit(x, model, model_predict)
  save(file = file.path(output.fig.dir, "model.Rdata"), cell.models)
} else {
  warning("Model file already exists, model not fit")
}
