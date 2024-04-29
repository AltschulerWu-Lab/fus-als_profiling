################################################################################
#+ setup, message=FALSE
################################################################################
library(tidyverse)
library(data.table)
library(tidytext)
library(scales)
library(ggsci)
library(ggrepel)
library(PRROC)

# Set user args for model choice
args <- R.utils::commandArgs(asValues = TRUE)
cytosolic <- as.logical(args$CYTOSOLIC) # FALSE
model.str <- args$MODEL #' irf'

if (!length(cytosolic)) cytosolic <- FALSE
if (!length(model.str)) model.str <- "irf"

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Scores")
platemap.dir <- file.path(fig.base.dir, "data", "platemaps")
dir.create(platemap.dir)

source(file.path(fig.base.dir, "paper_figures", "utilities.R"))
source(file.path(fig.base.dir, "paper_figures", "color_palette.R"))
source(file.path(fig.base.dir, "paper_figures", "theme.R"))
theme_set(theme_als())

# Set paths to data directories
data.dir <- file.path(data.base.dir, "032422_WTFUSSporadics", "Cell_Filtered")
library.file <- file.path(fig.base.dir, "data", "Fibroblasts_Library.csv")
plate.file <- file.path(fig.base.dir, "data", "plate_list.txt")

marker <- "FUS_EEA1"
save.fig <- TRUE
save.plate <- FALSE

importance_norm <- function(x) x / max(x)
ext <- ifelse(cytosolic, "_cytosolic", "_full")
output.fig.dir <- file.path(figure.dir, model.str, str_c("fig", ext))
dir.create(output.fig.dir, showWarnings = FALSE)

################################################################################
#' # Load data
#+ parameters, message=FALSE
################################################################################
# Load in datasets for selected marker
input.file <- str_c("profiles_", marker, ".Rdata")
load(file.path(data.dir, input.file))

# Initialize cell line metadata table
xcell <- fread(library.file) %>% dplyr::rename(CellLine = `Cell line`)

xcell <- dplyr::select(x, CellLine, CellLinePassageNumber, Genetics) %>%
  distinct() %>%
  left_join(xcell, by = "CellLine") %>%
  mutate(Sex = toupper(Sex)) %>%
  dplyr::rename(Passage = CellLinePassageNumber) %>%
  dplyr::rename(Site = `Origin Institution`)

if (!file.exists(file.path(output.fig.dir, "model.Rdata"))) {
  stop("Model not fit, run `train_fus_models.R`")
} else {
  load(file.path(output.fig.dir, "model.Rdata"))
}

################################################################################
#' # Predictions by cell line
#+ predictions_cell_line, message=FALSE
################################################################################
# Aggregate FUS v. Healthy predictions by cell line
predictions <- cell.models$ypred %>%
  group_by(CellLine, Ytest) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine")

################################################################################
# FUS-ALS vs. WT
################################################################################
# Compute test statistic for difference in prediction scores
x.wt.fus <- filter(predictions, Genetics != "sporadic")
fit <- lm(Ypred ~ Genetics + Sex + Site + Passage, data = x.wt.fus)
pval.lm <- summary(fit)$coefficients["GeneticsHealthy", "Pr(>|t|)"]
print(pval.lm)

yfus <- predictions$Ypred[predictions$Genetics == "FUS-ALS"]
ywt <- predictions$Ypred[predictions$Genetics == "Healthy"]
pval.sr <- wilcox.test(yfus, ywt, alternative = "greater")$p.value
print(pval.sr)

pval.lab <- label_pval(pval.sr)
ymax <- max(predictions$Ypred) + 1e-2

p <- predictions %>%
  filter(Genetics != "sporadic") %>%
  ggplot(aes(x = reorder(Genetics, Ypred, median), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics), outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  geom_text(x = 1.5, y = ymax + 5e-2, label = pval.lab, size = 8) +
  geom_segment(x = 1, xend = 2, y = ymax, yend = ymax, col = "#0B4668", linewidth = 1) +
  ylab("i-MAP score") +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = col.pal) +
  scale_color_manual(values = col.pal) +
  ylim(c(0, 1))

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_a.pdf"), h = 8, w = 4)
plot(p)
if (save.fig) dev.off()

################################################################################
# sporadic vs. WT
################################################################################
# Compute test statistic for difference in prediction scores
x.wt.sp <- filter(predictions, Genetics != "FUS-ALS")
fit <- lm(Ypred ~ Genetics + Sex + Site + Passage, data = x.wt.sp)
pval.lm <- summary(fit)$coefficients["Geneticssporadic", "Pr(>|t|)"]

ysporadic <- predictions$Ypred[predictions$Genetics == "sporadic"]
pval.sr <- wilcox.test(ysporadic, ywt, alternative = "greater")$p.value

pval.lab.s <- label_pval(pval.lm)
ymax.s <- max(predictions$Ypred[predictions$Genetics == "sporadic"]) + 5e-2

p <- predictions %>%
  ggplot(aes(x = reorder(Genetics, Ypred, median), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics), outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  geom_text(x = 1.5, y = ymax + 5e-2, label = pval.lab, size = 8) +
  geom_segment(x = 1, xend = 3, y = ymax, yend = ymax, col = "#0B4668", size = 1) +
  geom_text(x = 1.5, y = ymax.s + 5e-2, label = pval.lab.s, size = 8) +
  geom_segment(x = 1, xend = 2, y = ymax.s, yend = ymax.s, col = "#0B4668", size = 1) +
  ylab("i-MAP score") +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = col.pal) +
  scale_color_manual(values = col.pal) +
  ylim(0:1)

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_sa.pdf"), h = 8, w = 4)
plot(p)
if (save.fig) dev.off()

################################################################################
# Predictions vs. metadata features
################################################################################
p1 <- predictions %>%
  filter(Genetics == "Healthy") %>%
  ggplot(aes(x = Sex, y = Ypred, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  xlab("Sex") +
  theme(legend.position = "none") +
  ylab("i-MAP score") +
  scale_fill_hue(l = 55) +
  labs(fill = NULL) +
  ylim(0:1)

p2 <- predictions %>%
  filter(Genetics == "Healthy") %>%
  ggplot(aes(x = Site, y = Ypred, fill = Site)) +
  geom_boxplot(outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  xlab("Origin institution") +
  theme(legend.position = "none") +
  ylab(NULL) +
  scale_fill_jama() +
  labs(fill = NULL) +
  ylim(0:1)

p3 <- predictions %>%
  filter(Genetics == "Healthy") %>%
  ggplot(aes(x = Passage, y = Ypred)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(size = 3) +
  theme(legend.position = "none") +
  ylab(NULL) +
  xlab("Cell Line Passage Number") +
  ylim(0:1)

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_supp_metadata.pdf"), h = 8, w = 24)
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
if (save.fig) dev.off()

# Aggregate by cell line / plate
predictions <- cell.models$ypred %>%
  group_by(CellLine, Ytest, PlateID) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine") %>%
  dplyr::rename(Figure = `Figure Names`) %>%
  mutate(Plate = as.numeric(as.factor(PlateID))) %>%
  mutate(Plate = str_c("Plate ", Plate))

p4 <- predictions %>%
  filter(Genetics == "Healthy") %>%
  ggplot(aes(x = Plate, y = Ypred, fill = PlateID)) +
  geom_boxplot(outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  theme(legend.position = "none") +
  ylab("i-MAP score") +
  xlab(NULL) +
  scale_fill_hue(l = 55) +
  labs(fill = NULL) +
  ylim(0:1)

# Aggregate by cell line / plate
predictions <- cell.models$ypred %>%
  group_by(CellLine, Ytest, PlateID, WellID) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine") %>%
  dplyr::rename(Figure = `Figure Names`) %>%
  mutate(Row = str_remove_all(WellID, "-.*$")) %>%
  mutate(Row = as.numeric(Row)) %>%
  mutate(Col = str_remove_all(WellID, "^.*-")) %>%
  mutate(Col = as.numeric(Col)) %>%
  mutate(Plate = as.numeric(as.factor(PlateID))) %>%
  mutate(Plate = str_c("Plate ", Plate))

p5 <- predictions %>%
  mutate(Ypred = ifelse(Genetics == "Healthy", Ypred, NA)) %>%
  ggplot(aes(x = Col, y = Row, fill = Ypred)) +
  geom_tile() +
  ylab("Plate Column") +
  xlab("Plate Row") +
  scale_fill_viridis_c(limits = 0:1) +
  facet_wrap(~Plate) +
  labs(fill = "i-MAP score")


if (save.fig) pdf(file = file.path(output.fig.dir, "fig_supp_plate.pdf"), h = 8, w = 8)
plot(p4)
if (save.fig) dev.off()

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_supp_well.pdf"), h = 12, w = 24)
plot(p5)
if (save.fig) dev.off()


################################################################################
#' # Predictions by well replicate
#+ predictions_by_well, message=FALSE
################################################################################
# Aggregate predictions by well replicate
predictions <- cell.models$ypred %>%
  group_by(CellLine, WellID, Ytest, PlateID) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine") %>%
  dplyr::rename(Figure = `Figure Names`)

################################################################################
# FUS-ALS vs. WT
################################################################################
predictions.agg <- filter(predictions, Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(CellLine, Genetics) %>%
  summarize(Ypred = median(Ypred), .groups = "drop") %>%
  arrange(Ypred) %>%
  mutate(Face = ifelse(CellLine %in% c("F12", "F1"), "bold", "plain"))

# Plot predictions, FUS v. Healthy
p <- filter(predictions, Genetics %in% c("FUS-ALS", "Healthy")) %>%
  ggplot(aes(x = reorder(Figure, Ypred, median), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics), outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90, face = predictions.agg$Face)) +
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab("i-MAP score") +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  ylim(c(0, 1))

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_b.pdf"), h = 8, w = 16)
plot(p)
if (save.fig) dev.off()

################################################################################
# sporadic vs. WT
################################################################################
n.sporadic.top <- 5

# plot predictions, sporadic v. Healthy
predictions.agg <- filter(predictions, Genetics %in% c("sporadic", "Healthy")) %>%
  group_by(CellLine, Genetics) %>%
  summarize(Ypred = median(Ypred), .groups = "drop") %>%
  arrange(Ypred)

p <- filter(predictions, Genetics %in% c("sporadic", "Healthy")) %>%
  ggplot(aes(x = reorder(Figure, Ypred, median), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics), outlier.shape = NA, coef = NULL) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 90, face = predictions.agg$Face)) +
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab("i-MAP score") +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal) +
  ylim(c(0, 1)) 

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_sb.pdf"), h = 8, w = 16)
plot(p)
if (save.fig) dev.off()

################################################################################
#' # ROC curves
#+ roc, message=FALSE
################################################################################
predictions.cell <- filter(cell.models$ypred, Genetics %in% c("FUS-ALS", "Healthy"))

# Aggregate predictions by single well
predictions.well <- filter(cell.models$ypred, Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(CellLine, WellID, Ytest, PlateID) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine")

# Aggregate predictions by cell line
predictions.line <- filter(cell.models$ypred, Genetics %in% c("FUS-ALS", "Healthy")) %>%
  group_by(CellLine, Ytest) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine")

# Compute ROC curves
roc.cell <- roc.curve(
  scores.class0 = predictions.cell$Ypred,
  weights.class0 = predictions.cell$Ytest,
  curve = TRUE
)

roc.well <- roc.curve(
  scores.class0 = predictions.well$Ypred,
  weights.class0 = predictions.well$Ytest,
  curve = TRUE
)

roc.line <- roc.curve(
  scores.class0 = predictions.line$Ypred,
  weights.class0 = predictions.line$Ytest,
  curve = TRUE
)

# Initialize type labels
type.levels <- c("Single cell", "Well replicate", "Cell line")
aucs <- c(roc.cell$auc, roc.well$auc, roc.line$auc) %>% round(2)
type.levels <- str_c(type.levels, ", AUROC = ", aucs)

# Plot ROC curves
xplot <- rbind(
  data.frame(roc.cell$curve[, 1:2], Type = type.levels[1]),
  data.frame(roc.well$curve[, 1:2], Type = type.levels[2]),
  data.frame(roc.line$curve[, 1:2], Type = type.levels[3])
)

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_c.pdf"), h = 8, w = 8)
mutate(xplot, Type = factor(Type, levels = type.levels)) %>%
  ggplot(aes(x = X1, y = X2, col = Type)) +
  geom_line(size = 2) +
  geom_abline(col = "grey", lty = 2) +
  scale_color_jama() +
  xlab("False positive rate") +
  ylab("True positive rate") +
  theme(legend.position = c(0.75, 0.15)) +
  theme(legend.text = element_text(size = 16)) +
  labs(color = NULL)
if (save.fig) dev.off()


################################################################################
#' # Feature importance
#+ importance, message=FALSE
################################################################################
clean_marker <- function(x, markers = c("2" = "FUS", "3" = "EEA1")) {
  # Map channel indices to marker names
  channels <- str_split(x, "\\.\\.") %>%
    sapply("[", 2) %>%
    str_split("")
  return(sapply(channels, function(ch) str_c(markers[ch], collapse = "/")))
}

clean_region <- function(x) {
  # Map region indicators to full names
  region.key <- c(D = "nucleus", N = "cytoplasm", C = "cell")
  x <- str_remove_all(x, "^.*_")
  regions <- str_split(x, "\\.\\.") %>%
    sapply("[", 1) %>%
    str_split("")
  return(sapply(regions, function(r) str_c(region.key[r], collapse = "/")))
}

clean_feature <- function(x) {
  # Clean feature name, remove region / marker indicators
  x <- str_remove_all(x, "_(D|C|N)*\\.\\.[0-9]*$")
  return(str_replace_all(x, "_", " ") %>% str_c("\n"))
}

# Max-normalize feature importance for each hold-out cell line
importance <- apply(abs(cell.models$importance), MAR = 2, importance_norm)
features <- rownames(cell.models$importance) %>% str_remove_all("^X\\.\\.")

# Compute maximum importance within feature category, average across models
importance.avg <- reshape2::melt(importance) %>%
  mutate(Feature = features[Var1]) %>%
  dplyr::rename(Model = Var2) %>%
  mutate(Feature = str_remove_all(Feature, "\\.\\.[0-9]*$")) %>%
  group_by(Model, Feature) %>%
  summarize(Importance = max(value), .groups = "drop") %>%
  group_by(Feature) %>%
  summarize(SD = sd(Importance), Importance = mean(Importance)) %>%
  mutate(Marker = clean_marker(Feature)) %>%
  mutate(Region = clean_region(Feature)) %>%
  mutate(Feature = clean_feature(Feature)) %>%
  arrange(desc(Importance))

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_d.pdf"), h = 7, w = 8)
top_n(importance.avg, 5, Importance) %>%
  mutate(Feature = str_c(Feature, Region, ", ", Marker)) %>%
  ggplot(aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#0B4668", color = "black", width = 0.5) +
  geom_errorbar(aes(ymin = Importance, ymax = Importance + SD), width = 0.25) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  coord_flip() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle(NULL) +
  theme(text = element_text(size = 22))
if (save.fig) dev.off()

fout <- file.path(output.fig.dir, "importance.csv")
write.csv(file = fout, importance.avg)

################################################################################
#' # Age correlation
#+ age, message=FALSE
################################################################################
# Aggregate model predictions by cell line
predictions <- cell.models$ypred %>%
  dplyr::rename(Plate = PlateID) %>%
  group_by(CellLine, Ytest, Genetics, WellID, Plate) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop")

# Plot age v. predictions by genetics
xcell.s <- dplyr::select(xcell, -Genetics) %>%
  dplyr::rename(Age = `Age of biopsy`) %>%
  dplyr::rename(AgeOfOnset = `Age of onset`) %>%
  dplyr::rename(Figure = `Figure Names`)

xplot <- predictions %>% left_join(xcell.s, by = "CellLine")

xplot.avg <- group_by(predictions, CellLine, Genetics) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell.s, by = "CellLine") %>%
  mutate(AgeOfOnset = ifelse(AgeOfOnset == "Control", Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset = ifelse(AgeOfOnset == "control", Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset = ifelse(AgeOfOnset == "", Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset = as.numeric(AgeOfOnset)) %>%
  mutate(Figure = str_remove_all(Figure, "-.*$")) %>%
  dplyr::rename(Site = Site)

# Compute age/prediction correlations within genetic background
xfus <- filter(xplot.avg, Genetics == "FUS-ALS")
cor.fus <- cor.test(xfus$Ypred, xfus$Age, alternative = "less")

xwt <- filter(xplot.avg, Genetics == "Healthy")
cor.wt <- cor.test(xwt$Ypred, xwt$Age, alternative = "less")

xsporadic <- filter(xplot.avg, Genetics == "sporadic")
cor.spor <- cor.test(xsporadic$Ypred, xsporadic$Age, alternative = "less")

rhos <- round(c(cor.fus$estimate, cor.wt$estimate, cor.spor$estimate), 2)
pvals <- c(cor.fus$p.value, cor.wt$p.value, cor.spor$p.value)
pvals <- p.adjust(pvals, method='bonferroni')
pval.lab <- str_c("(", sapply(pvals, label_pval), ")")

xplot.text <- data.frame(Genetics = c("FUS-ALS", "Healthy", "sporadic")) %>%
  mutate(text = str_c("Cor = ", rhos, pval.lab)) %>%
  mutate(Ypred = c(0.15, 0.1, 0.05)) %>%
  mutate(Age = 55)

# Plot predictions v. age, FUS v. Healthy
if (save.fig) pdf(file = file.path(output.fig.dir, "fig_supp_aob.pdf"), h = 10, w = 10)
filter(xplot.avg, Genetics %in% c("Healthy", "FUS-ALS")) %>%
  ggplot(aes(x = Age, y = Ypred, col = Genetics)) +
  geom_point(size = 3) +
  geom_line(stat = "smooth", method = "lm", formula = y ~ x, alpha = 0.7, size = 1.5) +
  geom_text_repel(aes(label = CellLine), size = 7, min.segment.length = 0.01) +
  geom_text(data = filter(xplot.text, Genetics != "sporadic"), aes(label = text), size = 10) +
  scale_color_manual(values = col.pal) +
  ylim(0:1) +
  ylab("i-MAP score") +
  theme(legend.position = "none")
if (save.fig) dev.off()

if (save.fig) pdf(file = file.path(output.fig.dir, "fig_supp_sporadic_aob.pdf"), h = 10, w = 10)
filter(xplot.avg, Genetics %in% c("Healthy", "sporadic")) %>%
  ggplot(aes(x = Age, y = Ypred, col = Genetics)) +
  geom_point(size = 3) +
  geom_line(stat = "smooth", method = "lm", formula = y ~ x, alpha = 0.7, size = 1.5) +
  geom_text_repel(aes(label = CellLine), size = 7, min.segment.length = 0.01) +
  geom_text(data = filter(xplot.text, Genetics != "FUS-ALS"), aes(label = text), size = 10) +
  scale_color_manual(values = col.pal) +
  ylim(0:1) +
  ylab("i-MAP score") +
  theme(legend.position = "none")
if (save.fig) dev.off()


################################################################################
# Age of onset
################################################################################
# Compute age/prediction correlations within genetic background
cor.fus.aoo <- cor.test(xfus$Ypred, xfus$AgeOfOnset, alternative = "less")

rhos <- round(c(cor.fus.aoo$estimate, cor.wt$estimate, cor.spor$estimate), 2)
pvals <- c(cor.fus.aoo$p.value, cor.wt$p.value, cor.spor$p.value)
pvals <- p.adjust(pvals, method='bonferroni')
pval.lab <- str_c("(", sapply(pvals, label_pval), ")")

xplot.text <- data.frame(Genetics = c("FUS-ALS", "Healthy", "sporadic")) %>%
  mutate(text = str_c("Cor = ", rhos, pval.lab)) %>%
  mutate(Ypred = c(0.15, 0.1, 0.05)) %>%
  mutate(AgeOfOnset = 55)

# Plot predictions v. age of onset, FUS v. Healthy
if (save.fig) pdf(file = file.path(output.fig.dir, "fig_g.pdf"), h = 12, w = 12)
filter(xplot.avg, Genetics %in% c("Healthy", "FUS-ALS")) %>%
  ggplot(aes(x = AgeOfOnset, y = Ypred, col = Genetics)) +
  geom_point(size = 3) +
  geom_line(stat = "smooth", method = "lm", formula = y ~ x, size = 2) +
  geom_text_repel(aes(label = CellLine), size = 7, min.segment.length = 0.01) +
  geom_text(data = filter(xplot.text, Genetics != "sporadic"), aes(label = text), size = 10) +
  scale_color_manual(values = col.pal) +
  ylim(0:1) +
  ylab("i-MAP score") +
  xlab("Age of onset") +
  theme(legend.position = "none")
if (save.fig) dev.off()


################################################################################
# Save table of predictions by cell line
################################################################################
# Predictions by genetics
predictions <- cell.models$ypred %>%
  group_by(CellLine, Ytest) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  left_join(xcell, by = "CellLine")

save(file = file.path(output.fig.dir, "predictions.Rdata"), predictions)

################################################################################
# Save table of plates used
################################################################################
if (save.plate) {
  plate.list <- unique(predictions$Plate)

  write.table(
    file = plate.file,
    plate.list,
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # Save plate maps
  plate.layout <- cell.models$ypred %>%
    dplyr::select(PlateID, CellLine, WellID) %>%
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
