#' The following script fits models to predict ALS vs. health from sequencing 
#' data.
#' 
#' Requires input datasets
#'   `data/050322_RNAseq`
#'   `scripts/Figure1/models/predictions.Rdata`
#' 
#' **Note:** RNA-seq read count datasets for fibroblasts and spinal cord samples
#' can be processed from raw fastq files using the `als_rnaseq_fibro` and 
#' `als_rnaseq_spine` pipelines respectively. Data are loaded from calls to 
#' `load_fibro.R` and `load_spine.R`.
#' 
#' Karl Kumbier 3/8/2023
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(ggsci)
library(ggrepel)
library(hrbrthemes)
library(Matrix)
library(clusterProfiler)
library(patchwork)
library(glmnet)
library(PRROC)
options(stringsAsFactors=FALSE)

theme_set(
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=24, 
    axis_text_size=22,
    base_size=22,
    axis=TRUE,
    axis_col='#000000',
    base_family='sans',
    grid=FALSE
  )
)

################################################################################
# Load data and initialize analysis paramters
################################################################################
# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.dir <- file.path(fig.base.dir, 'data/050322_RNAseq/')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Transcriptomics')
mutation.file <- file.path(fig.base.dir, 'data', 'clinvar_als.txt')
library.path <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')

output.fig.dir <- file.path(figure.dir, 'fig')
dir.create(output.fig.dir, showWarnings=FALSE)

################################################################################
# Load data and initialize analysis parameters
################################################################################
source(file.path(figure.dir, 'utilities.R'))
source(file.path(figure.dir, 'load_fibroblast.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))


# Transformation to gene expression prior to PCA analysis
gene_transform <- function(x) rank(x)
model.type <- 'eigengene'

heat.pal <- pal_material('blue-grey')(10)
disc.col.pal <- c(pal_jama()(7), '#FF827D')
save.fig <- TRUE


# Format metadata feature names
colnames(xmeta.fibro) <- str_replace_all(colnames(xmeta.fibro), ' ', '_')

# Initialize keys for converting between ENSEMBL / Symbol gene IDs
gene.key <- clusterProfiler::bitr(
  colnames(x.fibro),
  'ENSEMBL',
  'SYMBOL',
  'org.Hs.eg.db'
)

gene.key.se <- setNames(gene.key$SYMBOL, gene.key$ENSEMBL)
gene.key.es <- setNames(gene.key$ENSEMBL, gene.key$SYMBOL)

# Load in table of ALS-related mutations
xmut <- fread(mutation.file) %>% 
  filter(str_detect(`Condition(s)`, '(amyo|Amyo)'))

als.genes <- xmut$`Gene(s)` %>% unique

x <- x.fibro[xmeta.fibro$Genetics != 'sporadic',]
xmeta <- xmeta.fibro[xmeta.fibro$Genetics != 'sporadic',]

# Load in cell line library metadata
xcell <- fread(library.path) %>% 
  dplyr::rename(CellLine=`Cell line`) %>%
  dplyr::rename(Fig=`Figure Names`)

# Load pathway enrichment file from spinal cord data
load(file.path(output.fig.dir, 'spine_gsea.Rdata'))

################################################################################
# Gene expression for optimal marker set
################################################################################
deseq <- deseq_design(
  x,
  xmeta,
  design=~ALS + Site,
  fitType='mean'
)

res <- results(deseq, name='ALSTRUE')

# Plot optimal marker set gene expression
genes.select <- gene.key.es[c('FUS', 'EEA1', 'CD63')]

xselect <- data.frame(x[,c(genes.select)], CellLine=xmeta$CellLine) %>%
  reshape2::melt() %>%
  dplyr::rename(Gene=variable, TPM=value) %>%
  left_join(xmeta.fibro, by='CellLine') %>%
  mutate(Gene=gene.key.se[as.character(Gene)])

# Initialize p-values (and adjusted)
pvals <- res[genes.select,]$pvalue %>% round(3)
names(pvals) <- gene.key.se[genes.select]

padj <- res[genes.select,]$padj 
names(padj) <- gene.key.se[genes.select]
pval.lab <- label_pval(padj)

# Initialize table for p-value segments
xseg <- data.frame(p=pval.lab, Gene=names(pvals)) %>%
  mutate(ymax=apply(x[,genes.select], MAR=2, max)) %>%
  mutate(ymax=ymax + 0.05 * ymax) %>%
  mutate(ytext=ymax + 0.05 * ymax)

# Plot marginal distibutions of CD63, EEA1, FUS
p <- filter(xselect, Genetics != 'sporadic') %>%
  mutate(GeneticsL=ifelse(Genetics == 'FUS', 'FUS-ALS', Genetics)) %>%
  mutate(GeneticsL=ifelse(GeneticsL == '3', 'WT', GeneticsL)) %>%
  ggplot(aes(x=GeneticsL, y=TPM)) +
  geom_boxplot(aes(fill=Genetics), coef=NULL, outlier.shape=NA) +
  geom_jitter(height=0, width=0.1, alpha=0.6) +
  geom_segment(data=xseg, x=1, xend=2, col='#0B4668', aes(y=ymax, yend=ymax)) +
  geom_text(data=xseg, x=1.5, aes(y=ytext, label=p), size=5) +
  facet_wrap(~Gene, scales='free_y', ncol=1) +
  scale_fill_manual(values=col.pal.g) +
  theme(legend.position='none') +
  xlab(NULL) +
  ylab('Gene expression (TPM)')

if (save.fig) pdf(file.path(output.fig.dir, 'optimal_markers.pdf'), h=10, w=7)
plot(p)
if(save.fig) dev.off()

################################################################################
#' ## Sequencing models
#' Transcription model of ALS. We fit supervised models trained to predict 
#' disease status. We use a leave-one-out sample split to evaluate predictions 
#' on held-out WT & FUS cell lines.
################################################################################
#+ eigengene
if (model.type == 'als') {
  genes <- intersect(gene.key.es[als.genes], colnames(x))
  x <- x[,genes]
  p.ylab <- 'FUS-ALS phenotype score\nsequencing, ALS genes'
} else if (model.type == 'pathway') {
  genes <- filter(reactome.pathways, Name == "Eukaryotic Translation Elongation")
  genes <- intersect(unlist(genes$ENSEMBL), colnames(x))
  x <- x[,genes]
  p.ylab <- 'FUS-ALS phenotype score\nsequencing, enriched PW genes'
}  else {
  pca <- prcomp(apply(x, MAR=2, rank))
  x <- pca$x
  p.ylab <- 'FUS-ALS phenotype score\nsequencing, eigengenes'
}

# Initialize classification models
fit_ <- fit_lm_classification
predict_ <- predict_lm_classification

# Initialize subset of FUS/WT cell lines for binary classification
xmeta$Y <- as.numeric(xmeta$ALS)
cell.lines <- xmeta$CellLine

fit <- lapply(1:length(cell.lines), function(i) {
  set.seed(i)
  cell <- cell.lines[i]
  out <- fit_holdout(x, xmeta, cell, fit_, predict_, importance_)
  return(out)
})

# Initialize prediction table
ypred.pp <- dplyr::select(xmeta.fibro, Ypred, CellLine)

ypred <- rbindlist(fit) %>%
  group_by(CellLine, Genetics) %>%
  summarize(YpredSeq=mean(YpredSeq), .groups='drop') %>%
  filter(Genetics %in% c('Healthy', 'FUS-ALS')) %>%
  mutate(ALS=as.numeric(Genetics != 'Healthy')) %>%
  left_join(ypred.pp, by='CellLine') %>%
  filter(!is.na(Ypred)) %>%
  mutate(YpredSeqT=YpredSeq > max(YpredSeq[ALS == 0])) %>%
  mutate(YpredT=Ypred > max(Ypred[ALS == 0]))

# Compute recall and AUROC for model predictions
score.imaging <- str_c(sum(ypred$YpredT), ' / ', sum(ypred$ALS))
score.seq <- str_c(sum(ypred$YpredSeqT), ' / ', sum(ypred$ALS))

auc.img <- roc.curve(scores.class0=ypred$Ypred, weights.class0=ypred$ALS)$auc
auc.seq <- roc.curve(scores.class0=ypred$YpredSeq, weights.class0=ypred$ALS)$auc

# Plot imaging vs. sequencing predictions
p <- left_join(ypred, xcell, by='CellLine') %>%
  ggplot(aes(x=Ypred, y=YpredSeq, col=Genetics)) +
  geom_point(size=4) +
  geom_text_repel(aes(label=Fig), size=6) +
  geom_vline(xintercept=max(ypred$Ypred[ypred$ALS == 0]), lty=2, alpha=0.5) +
  geom_hline(yintercept=max(ypred$YpredSeq[ypred$ALS == 0]), lty=2, alpha=0.5) +
  scale_color_manual(values=col.pal) +
  theme(legend.position='none') +
  xlab('FUS-ALS phenotype score\nimaging') +
  ylab(p.ylab) +
  xlim(0:1) +
  ylim(0:1) +
  theme(plot.margin=ggplot2::margin(0, 0, 0, 0))

d1 <- ggplot(ypred, aes(fill=Genetics)) +
  geom_boxplot(aes(x=Ypred), outlier.shape=NA, coef=NULL) + 
  scale_fill_manual(values=col.pal) +
  xlim(0:1) +
  theme_void() + 
  theme(legend.position="none")

d2 <- ggplot(ypred, aes(fill=Genetics)) +
  geom_boxplot(aes(x=YpredSeq), outlier.shape=NA, coef=NULL) + 
  scale_fill_manual(values=col.pal) +
  theme_void() + 
  xlim(0:1) +
  theme(legend.position="none") +
  coord_flip()

pg <- d1 + plot_spacer() + p + d2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 0.5), heights = c(0.5, 4))

if (save.fig) pdf(file.path(output.fig.dir, str_c('models_', model.type, '.pdf')), h=12, w=12)
plot(pg)
if(save.fig) dev.off()


# Correlation tests for prediction scores
ypred.fus <- filter(ypred, Genetics == 'FUS-ALS')
cor.test(ypred.fus$Ypred, ypred.fus$YpredSeq, alternative='greater')
