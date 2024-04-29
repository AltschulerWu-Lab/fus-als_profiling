######################################################## i#######################
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
################################################################################
library(data.table)
library(tidyverse)
library(ggsci)
library(Matrix)
library(PRROC)
library(ggrepel)
library(org.Hs.eg.db)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

auroc <- function(yp, yt) roc.curve(scores.class0 = yp, weights.class0 = yt)$auc

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.dir <- file.path(fig.base.dir, "data/050322_RNAseq/")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Transcriptomics")
pathway.file <- file.path(data.dir, "reactome_pathway.Rdata")
gene.file <- file.path(fig.base.dir, "data", "transcriptomic_genes.csv")

output.fig.dir <- file.path(figure.dir, "fig")
dir.create(output.fig.dir, showWarnings = FALSE)

source(file.path(figure.dir, "utilities.R"))
source(file.path(fig.base.dir, "paper_figures", "utilities.R"))
source(file.path(fig.base.dir, "paper_figures", "color_palette.R"))
source(file.path(fig.base.dir, "paper_figures", "theme.R"))
theme_set(theme_als())

# Parameters for analysis
heat.pal <- pal_material("blue-grey")(10)
save.fig <- TRUE

################################################################################
#+ load_data
################################################################################
load(file.path(output.fig.dir, "fibro_models.Rdata"))
load(file.path(output.fig.dir, "spine_models.Rdata"))
load(pathway.file)

xgene <- fread(gene.file) %>%
  mutate(Disease = ifelse(!Disease %in% c("ALS", ""), "ND", Disease))

# Initialize pathway groups
groups <- setNames(pw.enriched$RootName, pw.enriched$Name)

# Initialize keys for converting between ENSEMBL / Symbol gene IDs
genes <- rbindlist(coefs)$Var1 %>%
  as.character() %>%
  str_subset("^ENS")

gene.key.se <- mapIds(
  x = org.Hs.eg.db,
  keys = genes,
  column = "SYMBOL",
  keytype = "ENSEMBL"
)

# Select top-performing pathway classifiers from spinal cord data
top.spine.pw <- filter(pw.enriched, RootName != "Other") %>%
  group_by(RootName) %>%
  top_n(1, AUC)

# Compute fibroblast AUROC for pathways / groups
auc.fibro <- filter(ypred, Genetics != "sporadic") %>%
  group_by(Name) %>%
  summarize(AUC = auroc(YpredSeq, ALS)) %>%
  arrange(desc(AUC))

auc.fibro.group <- mutate(auc.fibro, Group = groups[Name]) %>%
  group_by(Group) %>%
  summarize(AUC = max(AUC), Name = Name[which.max(AUC)]) %>%
  arrange(AUC)

# Generate pathway enrichment table
xtab <- left_join(pw.enriched, auc.fibro, by='Name') %>%
  dplyr::rename(AUROC_Spine=AUC.x) %>%
  dplyr::rename(AUROC_Fibroblast=AUC.y) %>%
  group_by(RootName) %>%
  arrange(desc(AUROC_Fibroblast), .by_group = TRUE) %>%
  ungroup() %>%
  dplyr::select(Name, RootName, Size, AUROC_Fibroblast, AUROC_Spine)

fout <- file.path(output.fig.dir, 'pw_enrichment_table.csv')
write.csv(file = fout, xtab, row.names = FALSE)

################################################################################
#' ## Distribution of predictions by pathway
#+ pw_dist, fig.height=12, fig.width=12
################################################################################
ypred.ens <- filter(ypred, Name %in% top.spine.pw$Name) %>%
  filter(Genetics != "sporadic") %>%
  group_by(CellLine, Genetics, Sex, Site) %>%
  summarize(YpredSeq = mean(YpredSeq), Ypred = mean(Ypred))

fit <- lm(YpredSeq ~ Genetics + Sex + Site, data = ypred.ens)
pval.lm <- summary(fit)$coefficients["GeneticsHealthy", "Pr(>|t|)"]
print(pval.lm)

yfus <- ypred.ens$YpredSeq[ypred.ens$Genetics == "FUS-ALS"]
ywt <- ypred.ens$YpredSeq[ypred.ens$Genetics == "Healthy"]
pval.sr <- wilcox.test(yfus, ywt, alternative = "greater")$p.value
print(pval.sr)

pval.lab <- label_pval(pval.sr)
ymax <- max(ypred.ens$Ypred) + 1e-2

fout <- file.path(output.fig.dir, "predictions_by_genetics.pdf")
if (save.fig) pdf(fout, height = 8, width = 4)
ggplot(ypred.ens, aes(x = Genetics, y = YpredSeq, fill = Genetics)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.6) +
  geom_text(x = 1.5, y = ymax + 5e-2, label = pval.lab, size = 8) +
  geom_segment(x = 1, xend = 2, y = ymax, yend = ymax, col = "#0B4668", linewidth = 1) +
  scale_fill_manual(values = col.pal) +
  ylab("t-MAP score") +
  xlab(NULL) +
  theme(legend.position = "none") +
  ylim(0:1)
if (save.fig) dev.off()

################################################################################
#' ## Imaging vs. sequencing predictions for select pathway
#+ pw_predictions, fig.height=12, fig.width=12
################################################################################
xthr <- filter(ypred.ens, Genetics == "Healthy") %>% arrange(desc(Ypred))
xthr <- xthr$Ypred[1]

ythr <- filter(ypred.ens, Genetics == "Healthy") %>% arrange(desc(YpredSeq))
ythr <- ythr$YpredSeq[1]

fout <- file.path(output.fig.dir, "pathway_predictions.pdf")
if (save.fig) pdf(fout, height = 10, width = 10)
ypred.ens %>%
  ggplot(aes(x = Ypred, y = YpredSeq, col = Genetics)) +
  geom_vline(xintercept = xthr, col = "grey", lty = 2, linewidth = 1.5) +
  geom_hline(yintercept = ythr, col = "grey", lty = 2, linewidth = 1.5) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = CellLine), size = 6) +
  scale_color_manual(values = col.pal) +
  theme(legend.position = "none") +
  xlab("i-MAP score") +
  ylab("t-MAP score") +
  xlim(0:1) +
  ylim(0:1)
if (save.fig) dev.off()


# Plot ensemblized predictions with sporadics
ypred.ens.full <- filter(ypred, Name %in% top.spine.pw$Name) %>%
  group_by(CellLine, Genetics, Sex, Site) %>%
  summarize(YpredSeq = mean(YpredSeq), Ypred = mean(Ypred))

fout <- file.path(output.fig.dir, "pathway_predictions_full.pdf")
if (save.fig) pdf(fout, height = 10, width = 10)
ypred.ens.full %>%
  ggplot(aes(x = Ypred, y = YpredSeq, col = Genetics)) +
  geom_vline(xintercept = xthr, col = "grey", lty = 2, linewidth = 1.5) +
  geom_hline(yintercept = ythr, col = "grey", lty = 2, linewidth = 1.5) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = CellLine), size = 6) +
  scale_color_manual(values = col.pal) +
  theme(legend.position = "none") +
  xlab("i-MAP score") +
  ylab("t-MAP score") +
  xlim(0:1) +
  ylim(0:1)
if (save.fig) dev.off()

################################################################################
#+ ## imaging and pathway predictions by cell line
#+ predictions_by_celline
################################################################################
xplot <- group_by(ypred, Name) %>%
  mutate(SeqThreshold = max(YpredSeq[ALS == 0])) %>%
  mutate(Threshold = max(Ypred[ALS == 0])) %>%
  ungroup() %>%
  mutate(SepSeq = YpredSeq > SeqThreshold) %>%
  mutate(Sep = Ypred > Threshold) %>%
  mutate(Group = groups[Name]) %>%
  group_by(Group, CellLine, Genetics) %>%
  summarize(Sep = any(Sep), SepSeq = any(SepSeq))

xplot.seq <- filter(xplot, SepSeq) %>%
  dplyr::select(CellLine, Group, Genetics)

xplot.img <- filter(xplot, Sep) %>%
  dplyr::select(CellLine, Group, Genetics) %>%
  mutate(Group = "Imaging")

cell.order <- dplyr::select(ypred, CellLine, Ypred, ALS) %>%
  distinct() %>%
  arrange(desc(Ypred))

fout <- file.path(output.fig.dir, "cell_line_sep.pdf")
if (save.fig) pdf(fout, height = 8, width = 10)
rbind(xplot.seq, xplot.img) %>%
  filter(Genetics != "sporadic") %>%
  mutate(Group = factor(Group, levels = c(auc.fibro.group$Group, "Imaging"))) %>%
  mutate(CellLine = factor(CellLine, levels = cell.order$CellLine)) %>%
  ggplot(aes(x = CellLine, y = Group)) +
  geom_point(size = 8) +
  theme(axis.text.x = element_text(color = col.pal["FUS-ALS"])) +
  geom_hline(yintercept = 8.5, lty = 2, col = "#363636") +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position = "none")
if (save.fig) dev.off()

################################################################################
#' ## Number of cell lines classified by modality
#+ num_cell_line, fig.height=12, fig.width=12
################################################################################
xplot <- auc.fibro.group %>% 
  mutate(Group = factor(Group, levels = Group)) %>%
  filter(Group != "Other")

fout <- file.path(output.fig.dir, "num_classified_auc.pdf")
if (save.fig) pdf(fout, height = 8, width = 2)
ggplot(xplot, aes(x = Group, y = AUC)) +
  geom_point(size = 6, color = "#0B4668") +
  geom_segment(aes(xend = Group, yend = 0.7), lwd = 3, color = "#0B4668") +
  coord_flip() +
  theme(legend.position = "none") +
  xlab(NULL) +
  theme(axis.text.y = element_blank()) +
  ylab("AUROC") +
  scale_y_continuous(breaks = c(0.7, 1), limits = c(0.7, 1))
if (save.fig) dev.off()

################################################################################
#' ## Minimal gene set
#' We sought to identify minimal gene sets cable of predicting ALS vs. health
#' within each PW group. This strategy was analogous to imaging: (i) define a
#' refined, pathway focused assay (ii) train classifier (iii) predict hold-out.
#'
#' As with imaging, we focused on sparse (minimal) set of readouts for
#' prediction. Focusing on spinal cord enriched PWs to prioritize ALS relevant
#' pathways.
#+ coeffs, fig.height=8, fig.width=12
################################################################################
#' ## Minimal gene set
coef.select.full <- rbindlist(coefs) %>%
  filter(str_detect(Var1, "^ENSG")) %>%
  mutate(Group = groups[Name]) %>%
  filter(Group != "Other") %>%
  group_by(Var1, Group, Name) %>%
  summarize(Stability = mean(AbsBeta != 0), AbsBeta = mean(AbsBeta)) %>%
  mutate(Var1 = as.character(Var1)) %>%
  group_by(Group, Var1) %>%
  summarize(Stability = max(Stability), AbsBeta = max(AbsBeta)) %>%
  mutate(SYMBOL = gene.key.se[Var1]) %>%
  filter(Stability >= 0.5)

# Initialize stability score / coefficient matrices
coef.select <- group_by(coef.select.full, Group) %>%
  top_n(10, Stability) %>%
  mutate(Group = factor(Group, levels = auc.fibro.group$Group))

xgene <- dplyr::rename(xgene, SYMBOL = Gene)[-14, ]
xgene <- mutate(xgene, Disease = factor(Disease, levels = c("ALS", "ND", "")))

# Plot gene fingerprints by pathway
gene.order <- left_join(coef.select, xgene, by = "SYMBOL") %>%
  group_by(SYMBOL) %>%
  summarize(Stability = max(Stability)) %>%
  arrange(desc(Stability))

fout <- file.path(output.fig.dir, "gene_fingerprint.pdf")
if (save.fig) pdf(fout, height = 8, width = 24)
mutate(coef.select, SYMBOL = factor(SYMBOL, levels = gene.order$SYMBOL)) %>%
  ggplot(aes(x = SYMBOL, y = Group)) +
  geom_point(aes(col = Stability), size=8) +
  geom_point(shape = 1, size=8) +
  scale_color_material("blue-grey", breaks=c(0.5, 1), limits=c(0.5, 1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.key.width = unit(2, "cm"))
if (save.fig) dev.off()


fout <- file.path(output.fig.dir, "gene_refs.pdf")
if (save.fig) pdf(fout, height = 3, width = 24)
mutate(xgene, SYMBOL = factor(SYMBOL, levels = gene.order$SYMBOL)) %>%
  ggplot(aes(x = SYMBOL, y = 1, shape = Disease)) +
  geom_point(size = 8) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_shape_manual(guide = "none", values = c(8, 10, 32)) +
  xlab(NULL) +
  ylab(NULL) +
  ylim(0.95, 1.05)
if (save.fig) dev.off()

# Plot gene fingerprints aggregated
fout <- file.path(output.fig.dir, "gene_fingerprint_agg.pdf")
if (save.fig) pdf(fout, height = 7, width = 20)
left_join(gene.order, xgene, by = "SYMBOL") %>%
  mutate(Disease = as.character(Disease)) %>%
  mutate(Disease = ifelse(Disease == "", "No reported ND", Disease)) %>%
  mutate(Disease = ifelse(is.na(Disease), "No reported ND", Disease)) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = rev(gene.order$SYMBOL))) %>%
  ggplot(aes(y = SYMBOL, x = Stability, col = Disease)) +
  geom_point(size=5) +
  geom_segment(aes(yend=SYMBOL), xend = 0.5, linewidth = 3) +
  scale_color_jama() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Stability score') +
  ylab(NULL) +
  labs(col=NULL) +
  theme(legend.position = c(0.1, 0.85)) +
  theme(legend.text = element_text(size = 18)) +
  coord_flip()
if (save.fig) dev.off()

# Write gene fingerprint table
xtab <- left_join(coef.select.full, xgene, by = "SYMBOL") %>%
  dplyr::rename(RootName = Group) %>%
  group_by(RootName) %>%
  arrange(desc(Stability), .by_group = TRUE) %>%
  dplyr::select(RootName, SYMBOL, Stability, `Protein name`, Function, Disease, Reference, `Author, Year`)

fout <- file.path(output.fig.dir, 'gene_enrichment_table.csv')
write.csv(file = fout, xtab, row.names = FALSE)