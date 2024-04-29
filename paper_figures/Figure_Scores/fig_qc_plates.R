#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(ggsci)
library(gridExtra)
library(caret)

theme_set(theme_bw(base_size = 22))

# Paths to data directories
als.dir <- Sys.getenv("ALS_PAPER")
figure.dir <- file.path(als.dir, "paper_figures", "Figure_Scores", "fig_qc")
dir.create(figure.dir, showWarnings = FALSE)

# Source utility functions
source(file.path(als.dir, "paper_figures", "utilities.R"))
source(file.path(als.dir, "paper_figures", "color_palette.R"))
library.file <- file.path(als.dir, "data", "Fibroblasts_Library.csv")

# Parameters for analysis
save.fig <- TRUE # save figures or generate in notebook
heat.pal <- pal_material("blue-grey")(10)

model <- irf
model_predict <- predict_irf

# Initialize cell line metadata table
xcell <- fread(library.file) %>%
  dplyr::rename(CellLine = `Cell line`) %>%
  mutate(Sex = toupper(Sex)) %>%
  dplyr::rename(Genetics = `Disease Gene`) %>%
  dplyr::rename(Figure = `Figure Names`) %>%
  mutate(Genetics = str_replace_all(Genetics, "WT", "Healthy")) %>%
  mutate(Genetics = str_replace_all(Genetics, "FUS", "FUS-ALS"))


#' # Plate generalization
#+ cell_counts
data.dir <- file.path(als.dir, "data_profiles", "032422_WTFUSSporadics", "Cell_Filtered")
load(file.path(data.dir, "profiles_FUS_EEA1.Rdata"))

plate_fit <- function(x,
                      model = irf,
                      model_predict = predict_irf) {
  # Wrapper function to call plate_fit_single predict over all plates

  # Initialize cell lines to evaluate
  plates <- unique(x$PlateID)

  # Rank-normalize features within plate
  ecdf <- function(x) rank(x) / length(x)
  xfeat <- ungroup(x) %>%
    dplyr::select(matches("(^X|PlateID)")) %>%
    group_by(PlateID) %>%
    mutate_if(is.numeric, ecdf) %>%
    ungroup() %>%
    dplyr::select(-PlateID)

  x <- dplyr::select(x, -matches("^X"))
  x <- cbind(x, xfeat)

  # Drop sporadic lines from model trai
  fits <- lapply(plates, plate_fit_single,
    x = x,
    model = model,
    model_predict = model_predict
  )

  # Initialize importance
  return(rbindlist(fits))
}

plate_fit_single <- function(x,
                             plate.test,
                             model = irf,
                             model_predict = predict_irf) {
  # Check for valid input
  print(plate.test)
  if (nrow(x) == 0) {
    return(NULL)
  }

  # Set training/test indices based on hold-out place
  id.test <- which(x$Plate %in% plate.test)
  id.train <- setdiff(1:nrow(x), id.test)

  if (length(id.test) == 0) {
    return(NULL)
  }

  # Set feature matrix and response vector
  y <- as.numeric(x$Genetics != "Healthy")
  xx <- dplyr::select(x, matches("^X"))

  # Take balanced class sample
  ytab <- data.frame(Y = as.factor(y)) %>%
    mutate(Idx = 1:n()) %>%
    filter(Idx %in% id.train)

  id.downsample <- downSample(ytab, as.factor(ytab$Y))

  # Fit model
  id.train <- intersect(id.train, id.downsample$Idx)
  fit <- fit_model(xx, y, id.train, model, model_predict, test.id = id.test)

  # Aggregate predictions across each model
  ypred <- dplyr::select(x[id.test, ], WellID, PlateID, CellLine, Genetics) %>%
    mutate(Ypred = fit$ypred) %>%
    mutate(Ytest = y[id.test])

  return(ypred)
}

x <- filter(x, Genetics != "sporadic")
plate.models <- plate_fit(x, model, model_predict)

################################################################################
# Plate generalization by cell line
################################################################################
# Aggregate FUS v. WT predictions by cell line
predictions <- plate.models %>%
  group_by(CellLine, Ytest, PlateID) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine")

# Compute test statistic for difference in prediction scores
p <- predictions %>%
  filter(Genetics %in% c("FUS-ALS", "Healthy")) %>%
  mutate(PlateID = str_remove_all(PlateID, "^.*_plate_")) %>%
  ggplot(aes(x = reorder(Genetics, Ypred, median), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics), outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  ylab("i-MAP score") +
  xlab(NULL) +
  facet_wrap(~PlateID, ncol = 1) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = col.pal) +
  scale_color_manual(values = col.pal) +
  ylim(0:1)


if (save.fig) pdf(file = file.path(figure.dir, "plate_generalization_a.pdf"), h = 16, w = 4)
plot(p)
if (save.fig) dev.off()

################################################################################
# Plate generalization by well replicate
################################################################################
# Aggregate predictions by well replicate
predictions <- plate.models %>%
  group_by(CellLine, WellID, Ytest, PlateID) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine")

predictions.agg <- filter(predictions, Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(CellLine, Genetics, PlateID) %>%
  summarize(Ypred = median(Ypred), .groups = "drop") %>%
  arrange(Ypred)

# Plot predictions, FUS v. WT
p <- filter(predictions, Genetics %in% c("FUS-ALS", "Healthy")) %>%
  mutate(PlateID = str_remove_all(PlateID, "^.*_plate_")) %>%
  ggplot(aes(x = reorder(Figure, Ypred, median), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics), outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab("i-MAP score") +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  facet_wrap(~PlateID, ncol = 1) +
  ylim(c(0, 1))

if (save.fig) pdf(file = file.path(figure.dir, "plate_generalization_b.pdf"), h = 16, w = 16)
plot(p)
if (save.fig) dev.off()
