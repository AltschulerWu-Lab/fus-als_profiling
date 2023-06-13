#' Runs DESeq and GSEA on patient derived spinal cord RNA-seq data.
#' 
#' Requires input datasets:
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
library(Matrix)
options(stringsAsFactors=FALSE)

################################################################################
# Load data and initialize analysis parameters
################################################################################
# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.dir <- file.path(fig.base.dir, 'data/050322_RNAseq/')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure3')
output.fig.dir <- file.path(figure.dir, 'fig3')

################################################################################
# Load data and initialize analysis parameters
################################################################################
source(file.path(figure.dir, 'utilities.R'))
source(file.path(figure.dir, 'load_spine.R'))
source(file.path(figure.dir, 'load_fibroblast.R'))

# analysis parameters
genes <- intersect(colnames(x.spine), colnames(x.fibro))
design <- '~ALS + Site + Sex + Age'

# Format metadata feature names
colnames(xmeta.spine) <- str_replace_all(colnames(xmeta.spine), ' ', '_')
xmeta.spine <- mutate_if(xmeta.spine, is.character, as.factor)

################################################################################
# DESeq on spinal cord samples
################################################################################
deseq.spine <- deseq_design(
  x.spine[,genes], 
  xmeta.spine, 
  design=design, 
  fitType='mean'
) 

res <- results(deseq.spine, name='ALSTRUE')

# spine gene set enrichment analysis
score <- setNames(-log10(res$pvalue), rownames(res)) %>% na.omit
gsea.spine <- gsea(score, pathway.file, minSize=25, maxSize=100, pthr=0.05)

save(file=file.path(output.fig.dir, 'spine_gsea.Rdata'), gsea.spine)
