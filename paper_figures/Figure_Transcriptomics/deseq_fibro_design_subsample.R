#' In this script, we investigate the stability of DESeq results relative to 
#' design choice. This analysis is used to select the DESeq design that results
#' in the most consistent p-values for FUS-ALS fibroblast samples across data 
#' subsample replicates.
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
library(hrbrthemes)
library(Matrix)

options(stringsAsFactors=FALSE)

theme_set(
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=24, 
    axis_text_size=22,
    base_size=22,
    axis=TRUE,
    base_family='sans',
    axis_col='#000000',
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
output.fig.dir <- file.path(figure.dir, 'fig')
dir.create(output.fig.dir, showWarnings=FALSE)

################################################################################
# Load data and initialize analysis parameters
################################################################################
source(file.path(figure.dir, 'utilities.R'))
source(file.path(figure.dir, 'load_fibroblast.R'))
source(file.path(figure.dir, 'load_spine.R'))

# Threshold for sporadic groupings
save.fig <- TRUE
disc.pal <- c(pal_jama()(7), '#F57282')
n.subsample <- 25
genes <- intersect(colnames(x.fibro), colnames(x.spine))

# Format metadata feature names
colnames(xmeta.fibro) <- str_replace_all(colnames(xmeta.fibro), ' ', '_')
xmeta.fibro <- xmeta.fibro %>% mutate_if(is.character, as.factor)

xmeta.fibro <- xmeta.fibro %>%
  dplyr::rename(Age=Age_of_biopsy) %>%
  dplyr::rename(Passage=CellLinePassageNumber)

# Fix type in cell line passage number
xmeta.fibro$CellLinePassageNumber[xmeta.fibro$CellLine == 'ND40077'] <- 20

# Initialize IDs and genes for deseq
id.wt <- xmeta.fibro$Genetics == 'Healthy'
id.fus <- xmeta.fibro$Genetics == 'FUS-ALS'
id <- id.wt | id.fus

# Initialize designs as all possible combinations of control factors
# Note: we include site in all designs due to strong effect identified in PCA
factors <- c('Sex', 'Passage', 'Age')

factor.combs <- lapply(1:length(factors), function(k) {
  z <- combn(factors, k, simplify=FALSE)
  z <- lapply(z, function(zz) c('Site', zz))
})

factor.combs <- unlist(factor.combs, recursive=FALSE)
designs <- sapply(factor.combs, str_c, collapse=' + ')
designs <- c('Site', sapply(factor.combs, str_c, collapse=' + '))
designs <- str_c('~Genetics + ', designs)

################################################################################
#' # Stability analysis of DESeq designs
#' To evaluate different designs, we take subsamples of the full data
#' set, compute deseq p-values over each subsample, and compare the correlation 
#' of subsample replicate p-values generated from a fixed design.
################################################################################
#+ stability
set.seed(47)

# run deseq over subsampled data
deseq.ss <- replicate(n.subsample, {
  subsample_deseq(x.fibro[id, genes], xmeta.fibro[id,], designs)
}, simplify=FALSE)

deseq.ss <- lapply(designs, function(d) {
  return(sapply(deseq.ss, function(z) z[,d]))
})

# compute correlation between subsample -log10 p-values
cor.mats <- lapply(deseq.ss, function(z) cor(z, use='complete.obs'))
pw.cors <- lapply(cor.mats, function(z) z[upper.tri(z)])

if (save.fig) pdf(file=file.path(output.fig.dir, 'fibro_stability.pdf'), h=12, w=12)
reshape2::melt(pw.cors) %>%
  mutate(Design=designs[L1]) %>%
  filter(str_detect(Design, 'Site')) %>%
  group_by(Design) %>%
  summarize(Mean=mean(value), SEM=sd(value) / sqrt(length(value))) %>%
  ggplot(aes(x=reorder(Design, Mean), y=Mean, fill=Design)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Mean - SEM, ymax=Mean + SEM), width=0.2) +
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab('Pearson correlation') +
  labs(fill=NULL) +
  scale_fill_manual(values=disc.pal) +
  scale_y_continuous(limits=c(0.5, 0.65), oob=rescale_none) +
  theme(legend.position=c(0.275, 0.9))
if (save.fig) dev.off()
