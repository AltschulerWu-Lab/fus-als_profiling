################################################################################
#+ setup, message=FALSE
################################################################################
library(tidyverse)
library(data.table)
library(tidytext)
library(ggrepel)
theme_set(theme_bw(base_size = 22))

save.fig <- TRUE
age.of.onset <- TRUE

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Transcriptomics")
output.fig.dir <- file.path(figure.dir, "fig")

source(file.path(fig.base.dir, "paper_figures", "utilities.R"))
source(file.path(fig.base.dir, "paper_figures", "color_palette.R"))
library.file <- file.path(fig.base.dir, "data", "Fibroblasts_Library.csv")

load(file.path(output.fig.dir, "fibro_pw_models_auroc.Rdata"))
load(file.path(output.fig.dir, "spine_auc_enrich_group.Rdata"))

# Initialize top-leve pathway grouping
groups <- setNames(pw.enriched$RootName, pw.enriched$Name)

# Initialize cell line metadata table
xcell <- fread(library.file) %>%
  dplyr::rename(CellLine = `Cell line`) %>%
  dplyr::rename(Age = `Age of biopsy`) %>%
  dplyr::rename(AgeOfOnset = `Age of onset`) %>%
  dplyr::rename(Figure = `Figure Names`) %>%
  dplyr::select(CellLine, Age, AgeOfOnset, Figure)

# Select top-performaing pathway classifiers from spinal cord data
top.spine.pw <- group_by(top.pw, RootName) %>%
  mutate(Count = n()) %>%
  filter(Count >= 3) %>%
  top_n(1, AUC) %>%
  filter(!is.na(RootName))

# Initialize enzemblized fibroblast predictions
ypred.ens <- filter(ypred, Name %in% top.spine.pw$Name) %>%
  filter(Genetics != "sporadic") %>%
  group_by(CellLine, Genetics, Sex, Site) %>%
  summarize(YpredSeq = mean(YpredSeq), Ypred = mean(Ypred)) %>%
  ungroup()

# Merge predictions with cell metadata
predictions <- left_join(ypred.ens, xcell, by = "CellLine") %>%
  dplyr::select(-Ypred) %>%
  dplyr::rename(Ypred = YpredSeq) %>%
  mutate(AgeOfOnset = ifelse(AgeOfOnset == "Control", Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset = ifelse(AgeOfOnset == "control", Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset = ifelse(AgeOfOnset == "", Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset = as.numeric(AgeOfOnset)) %>%
  mutate(Figure = str_remove_all(Figure, "-.*$"))

if (age.of.onset) predictions$Age <- predictions$AgeOfOnset

# Compute age / prediction correlation by group
xplot.text <- group_by(predictions, Genetics) %>%
  summarize(
    Cor = cor.test(Ypred, Age, alternative = "less")$estimate,
    CorP = cor.test(Ypred, Age, alternative = "less")$p.value
  ) %>%
  mutate(CorP = p.adjust(CorP, method = "bonferroni")) %>%
  mutate(Cor = round(Cor, 3)) %>%
  mutate(CorP = round(CorP, 3)) %>%
  mutate(PLab = sapply(CorP, label_pval)) %>%
  mutate(Ypred = ifelse(Genetics == "Healthy", 0.95, 0.85)) %>%
  mutate(text = str_c("Cor = ", Cor, PLab))

# Plot predictions v. age, FUS v. Healthy
fout <- file.path(output.fig.dir, "fig_supp_rna_aob.pdf")
if (save.fig) pdf(file = fout, h = 8, w = 8)
xplot.text$Age <- 60
filter(predictions, Genetics %in% c("Healthy", "FUS-ALS")) %>%
  ggplot(aes(x = Age, y = Ypred, col = Genetics)) +
  geom_point(size = 3) +
  geom_line(stat = "smooth", method = "lm", formula = y ~ x, size = 2) +
  geom_text(data = filter(xplot.text, Genetics != "sporadic"), aes(label = text), size = 8) +
  geom_text_repel(aes(label = CellLine)) +
  scale_color_manual(values = col.pal) +
  ylim(0:1) +
  ylab("t-MAP score") +
  theme(legend.position = "none")
if (save.fig) dev.off()
