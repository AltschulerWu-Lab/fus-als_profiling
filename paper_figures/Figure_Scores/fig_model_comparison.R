################################################################################
#+ setup, message=FALSE
################################################################################
library(tidyverse)
library(data.table)
library(ggsci)
library(PRROC)

col.pal <- pal_jama()(7)[c(2, 6, 7)]

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Scores")
output.fig.dir <- file.path(figure.dir, "fig_comparisons")
dir.create(output.fig.dir, recursive = FALSE)

source(file.path(fig.base.dir, "paper_figures", "theme.R"))
theme_set(theme_als())

save.fig <- TRUE
models <- c("fcnn", "logistic", "irf")
cytosolic <- FALSE

# Load predictions from each model
xpred <- lapply(models, function(m) {
  subdir <- ifelse(cytosolic, "fig_cytosolic", "fig_full")
  load(file.path(figure.dir, m, subdir, "model.Rdata"))
  return(mutate(cell.models$ypred, Model = m))
})


################################################################################
#' # ROC curves
#+ roc, message=FALSE
################################################################################
roc <- function(ypred, ytrue) {
  roc.curve(scores.class0 = ypred, weights.class0 = ytrue, curve = TRUE)
}

roc.cell <- rbindlist(xpred) %>%
  filter(Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(Model) %>%
  summarize(ROC = list(roc(Ypred, Ytest))) %>%
  mutate(Type = "Single cell")

# Aggregate predictions by single well
roc.well <- rbindlist(xpred) %>%
  filter(Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(CellLine, WellID, Ytest, PlateID, Model) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  group_by(Model) %>%
  summarize(ROC = list(roc(Ypred, Ytest))) %>%
  mutate(Type = "Well replicate")

# Aggregate predictions by cell line
roc.line <- rbindlist(xpred) %>%
  filter(Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(CellLine, Ytest, Model) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  group_by(Model) %>%
  summarize(ROC = list(roc(Ypred, Ytest))) %>%
  mutate(Type = "Cell line")

rocs <- rbind(roc.cell, roc.well, roc.line)

# Initialize type labels
aucs <- mutate(rocs, AUC = sapply(ROC, function(z) z$auc)) %>%
  mutate(AUC = round(AUC, 2)) %>%
  mutate(AUC = str_c("AUROC = ", AUC)) %>%
  mutate(X = 0.75) %>%
  mutate(Y = ifelse(Model == "fcnn", 0.2, 0.15)) %>%
  mutate(Y = ifelse(Model == "logistic", 0.1, Y)) %>%
  mutate(AUC=Model)

rocs <- mutate(rocs, FPR = sapply(ROC, function(z) z$curve[, 1])) %>%
  mutate(TPR = sapply(ROC, function(z) z$curve[, 2]))

xplot <- lapply(1:nrow(rocs), function(i) {
  out <- data.frame(FPR = rocs$FPR[[i]], TPR = rocs$TPR[[i]]) %>%
    mutate(Model = rocs$Model[i]) %>%
    mutate(Type = rocs$Type[i])
  return(out)
})

type.levels <- c("Single cell", "Well replicate", "Cell line")
aucs <- mutate(aucs, Type = factor(Type, levels = type.levels))
xplot <- rbindlist(xplot) %>% mutate(Type = factor(Type, levels = type.levels))

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_supp_roc_comparison.pdf"), h = 8, w = 16)
ggplot(xplot, aes(col = Model)) +
  geom_line(aes(x = FPR, y = TPR), size = 2) +
  geom_abline(col = "grey", lty = 2) +
  scale_color_manual(values = col.pal) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  theme(legend.position = "none") +
  labs(color = NULL) +
  facet_wrap(~Type) +
  geom_text(data = aucs, aes(x = X, y = Y, label = AUC))
if (save.fig) dev.off()


if (save.fig) pdf(file = file.path(output.fig.dir, "fig_roc_comparison_well.pdf"), h = 8, w = 8)
filter(xplot, Type == "Well replicate") %>%
  ggplot(aes(col = Model)) +
  geom_line(aes(x = FPR, y = TPR), size = 2) +
  scale_color_manual(values = col.pal) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  theme(legend.position = c(0.85, 0.15)) +
  labs(col = NULL) +
  theme(legend.text = element_text(size = 24))
if (save.fig) dev.off()
