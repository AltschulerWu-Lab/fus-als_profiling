# " # ASO analysis
# " The following script fits and evaluates FUS vs. WT classifiers on ASO-treated
# " cell lines. Models can be fit for the full feature set or cytosolic features
# " only by setting the variable `cytosolic  =  FALSE` or `cytosolic  =  TRUE`
# " respectively. Models can be fit using EEA1/FUS features or only FUS features
# " by setting the variable `fus.only  =  FALSE` or `fus.only  =  TRUE` respectively.
# "
# " Models are trained using data from EP cells and evaluated on EP cells treated
# " with NTC, low ASO and high ASO.
library(data.table)
library(tidyverse)
library(tidytext)
library(ggsci)
library(ggpubr)
library(caret)
library(gridExtra)
library(lme4)
library(lmerTest)

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_ASO")
meta.dir <- file.path(fig.base.dir, "data")

source(file.path(fig.base.dir, "paper_figures", "utilities.R"))
source(file.path(fig.base.dir, "paper_figures", "color_palette.R"))
source(file.path(fig.base.dir, "paper_figures", "theme.R"))
theme_set(theme_als())

# Set paths to data directories
data.dir <- file.path(data.base.dir, "120423_ASO14", "Cell_Filtered")

# Set parameters for analysis
n.core <- 16
marker <- "FUS_EEA1" # marker set used in screen
marker.select <- "FUS_EEA1" # marker(s) used in model
region <- "cytosolic"
drop.sporadic <- FALSE # if TRUE, plots FUS/healthy, else plots sporadic
save.fig <- TRUE
correction <- "bonferroni"
train.condition <- "untreated"
treatments <- c("untreated", "NTC 10uM", "ASO 1uM", "ASO 10uM")

if (drop.sporadic) {
  genetics.keep <- c("Healthy", "FUS-ALS")
} else {
  genetics.keep <- "sporadic"
}

# Set paths to data directories
params <- str_c(region, "_", marker, ":", marker.select)
output.fig.dir <- file.path(figure.dir, "fig")
dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(data.dir, str_c("profiles_", marker, ".Rdata")))

# Initialize models to be used for supervised learning
model <- function(x, y) irf(x, y, n.core = n.core)
model_predict <- predict_irf

# " ## Loading data
# " Our dataset consists of phenotypic profiles. Each
# " observation in the dataset corresponds to a single cell.
#+load_data
# Load in dataset for select marker set
x <- mutate(x, PlateID = str_remove_all(PlateID, ".*plate_")) %>%
  mutate(ID = str_c(ID, seq_len(n()))) %>%
  dplyr::rename(ASO = ASOtreatment) %>%
  mutate(ASO = str_remove_all(ASO, "EP_")) %>%
  mutate(ASO = str_replace_all(ASO, "ASO", " ASO")) %>%
  mutate(ASO = str_replace_all(ASO, "NTC", "NTC 10uM")) %>%
  mutate(ASO = str_replace_all(ASO, "low ASO", "ASO 1uM")) %>%
  mutate(ASO = str_replace_all(ASO, "high ASO", "ASO 10uM")) %>%
  mutate(ASO = factor(ASO, treatments))

# Initialize gene/cell line key
genetic.key <- dplyr::select(x, CellLine, Genetics) %>% distinct()

# Initialize channel / region IDs
channels <- sapply(str_split(colnames(x), "\\.\\."), "[", 3)
channels[is.na(channels)] <- 0

regions <- sapply(str_split(colnames(x), "\\.\\."), "[", 2)
regions <- lapply(regions, str_split, pattern = "_")
regions <- sapply(regions, function(z) tail(z[[1]], 1))
regions[is.na(regions)] <- "X"

# Filter to selected region / marker
id <- rep(TRUE, ncol(x))
if (region == "cytosolic") id <- id & regions %in% c("X", "N", "NN")
if (marker.select == "FUS") id <- id & channels %in% c("0", "2")
if (marker.select == "CD63") id <- id & channels %in% c("0", "2")
if (marker.select == "EEA1") id <- id & channels %in% c("0", "3")
x <- dplyr::select_if(x, id)

# " ## ASO models
# " We evaluate the degree to which ASO treatment "pushes" cell lines (disease
# " to health, health to disease) in the context of supervised models trained
# " to discriminate between ALS/FUS. Models are trained on DMSO treated cell
# " lines and evaluated on lines treated with (i) H2O (ii) NTC (iii) low ASO
# " and (iv) high ASO.
#+ leave_cell_out, fig.height = 12, fig.width = 16
aso_fit <- function(
    x,
    model = irf,
    model_predict = predict_irf) {
  # Check for valid input
  if (nrow(x) == 0) {
    return(NULL)
  }

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

  # Set feature matrix and reponse vector
  y <- as.numeric(x$Genetics != "Healthy")

  # Downsample for class balance
  xds <- downSample(x, as.factor(y))
  xx <- dplyr::select(xds, matches("^X"))
  y <- as.numeric(xds$Class) - 1

  # Set training/test indices
  id.train <- which(xds$ASO == train.condition)
  id.test <- which(xds$ASO != train.condition)
  fit <- fit_model(xx, y, id.train, model, model_predict)

  predicted <- data.frame(
    Ypred = fit$ypred,
    Ytest = y[id.test],
    Genetics = xds$Genetics[id.test],
    CellLine = xds$CellLine[id.test],
    Treatment = xds$ASO[id.test],
    WellID = xds$WellID[id.test],
    ID = xds$ID[id.test],
    Marker = marker,
    Plate = xds$PlateID[id.test]
  )

  return(list(fit = fit, predicted = predicted))
}


# Filter data to selected screen/compound
set.seed(47)
cell.models <- aso_fit(x, model, model_predict)

# " ### Predictions by genetics
# " Model prediction scores by genetics, treatment.
#+ raw_predictions, fig.height = 8, fig.width = 12
main <- str_c(marker.select, ", ", region)

xplot <- filter(cell.models$predicted) %>%
  group_by(Genetics, CellLine, Treatment) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop")

fout <- str_c("pred_", params, ".pdf")
fout <- file.path(output.fig.dir, fout)

if (save.fig) pdf(file = fout, height = 8, width = 8)
ggplot(xplot, aes(x = Treatment, y = Ypred)) +
  geom_boxplot(aes(fill = Genetics)) +
  scale_fill_manual(values = col.pal) +
  ylim(0:1) +
  ggtitle(main)
if (save.fig) dev.off()

#####################
# Genetic predictions
#####################
# Average predictions by well and normalize relative to NTC
xplot <- filter(cell.models$predicted, !Treatment %in% train.condition) %>%
  group_by(Genetics, CellLine, Treatment, WellID, Plate) %>%
  summarize(Ypred = mean(Ypred), .groups = "drop") %>%
  filter(Treatment != "untreated") %>%
  group_by(CellLine) %>%
  mutate(Ypred = Ypred / mean(Ypred[Treatment == "NTC 10uM"])) %>%
  ungroup()

# Initialize plot range
xrange <- group_by(xplot, Genetics, Treatment) %>%
  summarize(Mean = mean(Ypred), SD = sd(Ypred)) %>%
  mutate(MaxRange = Mean + SD) %>%
  mutate(MinRange = Mean - SD)

xrange <- c(min(xrange$MinRange), max(xrange$MaxRange))
eps <- diff(xrange) / 10

# Compute significance levels for treatment effect
lme.test <- lapply(unique(xplot$Genetics), function(g) {
  xg <- filter(xplot, Genetics == g)
  mixed <- lmerTest::lmer(Ypred ~ Treatment + (1 | CellLine), data = xg)
  out <- summary(mixed)$coefficients[3, 5]
  out <- data.frame(Treatment = "ASO 10uM", p = out) %>% 
    mutate(Treatment = str_remove_all(Treatment, "Treatment")) %>%
    mutate(Genetics = g)
  return(out)
})

test.plot <- rbindlist(lme.test) %>%
  mutate(padj = p.adjust(p, method = correction)) %>%
  mutate(padj = signif(padj, 3)) %>%
  mutate(x = ifelse(Treatment == "ASO 1uM", 1, 2)) %>%
  mutate(x = 2) %>%
  mutate(y = max(xrange) - eps) %>%
  mutate(y = ifelse(Treatment == "ASO 1uM", y + eps / 2, y)) %>%
  mutate(xend = x + 1) %>%
  mutate(yend = y) %>%
  mutate(xtext = 2) %>%
  mutate(ytext = y + eps / 4) %>%
  mutate(x = 1) %>%
  mutate(plab = label_pval(padj)) %>%
  filter(Genetics %in% genetics.keep)

p <- filter(xplot, Genetics %in% genetics.keep) %>%
  ggbarplot(
    x = "Treatment",
    fill = "Genetics",
    y = "Ypred",
    facet.by = "Genetics",
    merge = TRUE,
    add = "mean_se",
    error.plot = "upper_errorbar"
  ) +
  geom_text(
    data = filter(test.plot, str_detect(Treatment, "1uM")),
    aes(x = xtext, y = ytext, label = plab), size = 8
  ) +
  geom_text(
    data = filter(test.plot, str_detect(Treatment, "10uM")),
    aes(x = xtext, y = ytext, label = plab), size = 8
  ) +
  geom_segment(
    data = test.plot,
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = col.pal) +
  coord_cartesian(ylim = xrange) +
  ggtitle(main) +
  geom_hline(yintercept = 1, col = "grey", lty = 2) +
  ylab("i-MAP score, relative to NTC") +
  theme(text = element_text(size = 22))

w <- ifelse(drop.sporadic, 8, 5)
fout <- str_c("adj_cell_line_", params, "_", drop.sporadic, ".pdf")
fout <- file.path(output.fig.dir, fout)
if (save.fig) pdf(file = fout, height = 8, width = w)
plot(p)
if (save.fig) dev.off()


#######################
# Cell line predictions
#######################
# Initialize plot range
xrange <- group_by(xplot, CellLine, Treatment) %>%
  summarize(Mean = mean(Ypred), SD = sd(Ypred)) %>%
  mutate(MaxRange = Mean + SD) %>%
  mutate(MinRange = Mean - SD)

xrange <- c(min(xrange$MinRange), max(xrange$MaxRange))
eps <- diff(xrange) / 8

# Compute significance levels for treatment effect
rs.test <- lapply(unique(xplot$CellLine), function(cl) {
  xc <- filter(xplot, CellLine == cl)
  xt2 <- filter(xc, Treatment == "ASO 10uM")
  xt1 <- filter(xc, Treatment == "ASO 1uM")
  xt0 <- filter(xc, Treatment == "NTC 10uM")

  rs.20 <- wilcox.test(xt2$Ypred, xt0$Ypred, alternative = "less")
  rs.10 <- wilcox.test(xt1$Ypred, xt0$Ypred, alternative = "less")

  out <- data.frame(Treatment = c("ASO 10uM")) %>%
    mutate(p = c(rs.20$p.value)) %>%
    mutate(CellLine = cl)

  return(out)
})

test.plot <- rbindlist(rs.test) %>%
  mutate(padj = p.adjust(p, method = correction)) %>%
  mutate(padj = signif(padj, 3)) %>%
  mutate(x = ifelse(Treatment == "ASO 1uM", 1, 2)) %>%
  mutate(y = max(xrange) - eps) %>%
  mutate(y = ifelse(Treatment == "ASO 1uM", y + eps / 2, y)) %>%
  mutate(xend = x + 1) %>%
  mutate(yend = y) %>%
  mutate(xtext = 2) %>% # x + 0.5) %>%
  mutate(ytext = y + eps / 2) %>%
  mutate(x = 1) %>%
  mutate(plab = label_pval(padj))

if (drop.sporadic) {
  test.plot <- filter(test.plot, !str_detect(CellLine, "S"))
} else {
  test.plot <- filter(test.plot, str_detect(CellLine, "S"))
}

# Plot FUS, healthy
p <- filter(xplot, Genetics %in% genetics.keep) %>%
  ggbarplot(
    x = "Treatment",
    fill = "Genetics",
    y = "Ypred",
    facet.by = "CellLine",
    merge = TRUE,
    add = "mean_se",
    error.plot = "upper_errorbar",
    nrow = 1
  ) +
  geom_text(
    data = filter(test.plot, str_detect(Treatment, "1uM")),
    aes(x = xtext, y = ytext, label = plab), size = 7
  ) +
  geom_text(
    data = filter(test.plot, str_detect(Treatment, "10uM")),
    aes(x = xtext, y = ytext, label = plab), size = 7
  ) +
  geom_segment(
    data = test.plot,
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = col.pal) +
  coord_cartesian(ylim = xrange) +
  ggtitle(main) +
  ylab("i-MAP score, relative to NTC") +
  theme(text = element_text(size = 22))

fout <- str_c("well_f-h_", params, "_", drop.sporadic, ".pdf")
fout <- file.path(output.fig.dir, fout)
if (save.fig) pdf(file = fout, height = 8, width = 16)
plot(p)
if (save.fig) dev.off()
