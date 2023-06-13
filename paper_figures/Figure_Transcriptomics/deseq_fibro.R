#' Runs DESeq and GSEA on patient derived fibroblast RNA-seq data.
#'
#' Requires input datasets
#'   `data/050322_RNAseq`
#'   `data/rna_seq/nygc`
#'   `data/reactome_pathway.Rdata`
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
library(hrbrthemes)
library(Matrix)
options(stringsAsFactors=FALSE)


################################################################################
# Load data and initialize analysis paramters
################################################################################
# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.dir <- file.path(fig.base.dir, 'data/050322_RNAseq/')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure3')
output.fig.dir <- file.path(figure.dir, 'fig3')
dir.create(output.fig.dir, showWarnings=FALSE)

################################################################################
# Load data and initialize analysis parameters
################################################################################
source(file.path(figure.dir, 'utilities.R'))
source(file.path(figure.dir, 'load_fibroblast.R'))
source(file.path(figure.dir, 'load_spine.R'))

# Threshold for sporadic groupings
sp.thresh <- 5 # take top 5 sporadic lines for pw enrichment
genes <- intersect(colnames(x.fibro), colnames(x.spine))
design <- '~ALS + Site'

# Format metadata feature names
colnames(xmeta.fibro) <- str_replace_all(colnames(xmeta.fibro), ' ', '_')
xmeta.fibro <- xmeta.fibro %>% mutate_if(is.character, as.factor)

################################################################################
# FUS-ALS
################################################################################
id.wt <- xmeta.fibro$Genetics == 'Healthy'
id.fus <- xmeta.fibro$Genetics == 'FUS-ALS'
id <- id.wt | id.fus

# compute TPM deseq
deseq.fibro <- deseq_design(
  x.fibro[id, genes], 
  xmeta.fibro[id,],
  design=design,
  fitType='mean'
)

res <- results(deseq.fibro, name='ALSTRUE')

# GSEA comparison
score <- setNames(-log10(res$pvalue), rownames(res)) %>% na.omit
gsea.fus <- gsea(score, pathway.file, minSize=25, maxSize=100, pthr=0.05)

################################################################################
# Sporadic+
################################################################################
sporadic.high <- filter(xmeta.fibro, Genetics == 'sporadic') %>% 
  top_n(sp.thresh, Ypred)

sporadic.high <- sporadic.high$CellLine %>% as.character

id.wt <- xmeta.fibro$Genetics == 'Healthy'
id.sph <- xmeta.fibro$CellLine %in% sporadic.high
id <- id.wt | id.sph

# compute TPM deseq
deseq.sph <- deseq_design(
  x.fibro[id, genes], 
  xmeta.fibro[id,],
  design=design,
  fitType='mean'
)

res <- results(deseq.sph, name='ALSTRUE')

# GSEA comparison
score <- setNames(-log10(res$pvalue), rownames(res)) %>% na.omit
gsea.sph <- gsea(score, pathway.file, minSize=25, maxSize=100, pthr=0.05)

################################################################################
# Sporadic-
################################################################################
id.wt <- xmeta.fibro$Genetics == 'Healthy'
id.spl <- !xmeta.fibro$CellLine %in% sporadic.high 
id.spl <- id.spl & xmeta.fibro$Genetics == 'sporadic'
id <- id.wt | id.spl

# compute TPM deseq
deseq.spl <- deseq_design(
  x.fibro[id, genes], 
  xmeta.fibro[id,],
  design=design,
  fitType='mean'
)

res <- results(deseq.spl, name='ALSTRUE')

# GSEA comparison
score <- setNames(-log10(res$pvalue), rownames(res)) %>% na.omit
gsea.spl <- gsea(score, pathway.file, minSize=25, maxSize=100, pthr=0.05)

save(file=file.path(output.fig.dir, 'fibro_gsea.Rdata'), gsea.fus, gsea.sph, gsea.spl)
