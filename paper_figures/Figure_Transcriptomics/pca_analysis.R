#' In this script, we run PCA analysis on expression data and report primary
#' drivers of variation in fibroblast / spinal cord expression data.
#' 
#' Requires input datasets
#'   `data/050322_RNAseq`
#'   `data/rna_seq/nygc`
#'   `data/050322_RNAseq`
#'   `data/reactome_pathway.Rdata`
#'   `scripts/Figure1/models/predictions.Rdata`
#' 
#' **Note:** RNA-seq read count datasets for fibroblasts and spinal cord samples
#' can be processed from raw fastq files using the `als_rnaseq_fibro` and 
#' `als_rnaseq_spine` pipelines respectively. Data are loaded from calls to 
#' `load_fibro.R` and `load_spine.R`.
#' 
#' Karl Kumbier 3/3/2023
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(Matrix)
library(gridExtra)
library(ggsci)

options(stringsAsFactors=FALSE)
theme_set(theme_bw(base_size=18))

grid.arrange4 <- function(...) grid.arrange(..., ncol=4)

################################################################################
# Load data and initialize analysis paramters
################################################################################
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.dir <- file.path(fig.base.dir, 'data/050322_RNAseq/')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Transcriptomics')
output.fig.dir <- file.path(figure.dir, 'fig')
dir.create(output.fig.dir, showWarnings=FALSE)

################################################################################
# Load data and initialize analysis parameters
################################################################################
source(file.path(figure.dir, 'utilities.R'))
load(file.path(data.dir, 'fibro_processed.Rdata'))
load(file.path(data.dir, 'spine_processed.Rdata'))

# Initialize analysis parameters
save.fig <- TRUE
disc.pal <- c(pal_jama()(7), '#F57282')
gene_transform <- function(x) log(x + 1)

# Format metadata feature names
colnames(xmeta.fibro) <- str_replace_all(colnames(xmeta.fibro), ' ', '_')
xmeta.fibro <- xmeta.fibro %>% mutate_if(is.character, as.factor)

# Clean names of fibro metadata
xmeta.fibro <- xmeta.fibro %>%
  dplyr::rename(Age=Age_of_biopsy) %>%
  dplyr::rename(Passage=CellLinePassageNumber)

# Clean names of spine sites
xmeta.spine <- mutate(xmeta.spine, Site=str_split(Site, ' ')) %>%
  mutate(Site=lapply(Site, str_subset, pattern='of', negate=TRUE)) %>%
  mutate(Site=lapply(Site, str_sub, start=1, end=1)) %>%
  mutate(Site=sapply(Site, str_c, collapse=''))

################################################################################
# PCA analysis of TPM-normalized fibroblast expression
################################################################################
pca <- prcomp(gene_transform(x.fibro))
xpca <- cbind(pca$x, xmeta.fibro)
pct.var <- round(pca$sdev ^ 2 / sum(pca$sdev ^ 2), 3) * 100

meta.features <- c('Site', 'Sex', 'Genetics', 'Age')

p <- lapply(meta.features, function(f) {
  out1 <- ggplot(xpca, aes_string(col=f)) +
    geom_point(size=2, aes(x=PC1, y=PC2)) +
    scale_size(guide='none') +
    xlab(str_c('PC1 (', pct.var[1], '%)')) +
    ylab(str_c('PC2 (', pct.var[2], '%)'))
  
  out2 <- ggplot(xpca, aes_string(col=f)) +
    geom_point(size=2, aes(x=PC3, y=PC4)) +
    scale_size(guide='none') +
    xlab(str_c('PC3 (', pct.var[3], '%)')) +
    ylab(str_c('PC4 (', pct.var[4], '%)'))
  
  if (is.numeric(xpca[[f]])) {
    out1 <- out1 + scale_color_viridis_c()
    out2 <- out2 + scale_color_viridis_c()
    
  } else if (f == "Genetics") {
    out1 <- out1 + scale_color_manual(values = col.pal)
    out2 <- out2 + scale_color_manual(values = col.pal)
  } else {
    out1 <- out1 + scale_color_jama()
    out2 <- out2 + scale_color_jama()
  }
  
  return(list(out1, out2))
})

if (save.fig) pdf(file=file.path(output.fig.dir, 'pca_fibro.pdf'), h=12, w=24)
p <- unlist(p, recursive=FALSE)
do.call(grid.arrange4, p) 
if (save.fig) dev.off()

################################################################################
# PCA analysis of TPM-normalized spine expression
################################################################################
pca <- prcomp(gene_transform(x.spine))
xpca <- cbind(pca$x, xmeta.spine)
pct.var <- round(pca$sdev ^ 2 / sum(pca$sdev ^ 2), 3) * 100

meta.features <- c('Site', 'Sex', 'ALS', 'Age', 'RIN')

p <- lapply(meta.features, function(f) {
  out1 <- ggplot(xpca, aes_string(col=f)) +
    geom_point(size=2, aes(x=PC1, y=PC2)) +
    scale_size(guide='none') +
    xlab(str_c('PC1 (', pct.var[1], '%)')) +
    ylab(str_c('PC2 (', pct.var[2], '%)'))
  
  out2 <- ggplot(xpca, aes_string(col=f)) +
    geom_point(size=2, aes(x=PC3, y=PC4)) +
    scale_size(guide='none') +
    xlab(str_c('PC3 (', pct.var[3], '%)')) +
    ylab(str_c('PC4 (', pct.var[4], '%)'))
  
  if (is.numeric(xpca[[f]])) {
    out1 <- out1 + scale_color_viridis_c()
    out2 <- out2 + scale_color_viridis_c()
    
  }  else {
    out1 <- out1 + scale_color_manual(values=disc.pal)
    out2 <- out2 + scale_color_manual(values=disc.pal)
  }
  
  return(list(out1, out2))
})

if (save.fig) pdf(file=file.path(output.fig.dir, 'pca_spine.pdf'), h=12, w=24)
p <- unlist(p, recursive=FALSE)
do.call(grid.arrange4, p) 
if (save.fig) dev.off()
