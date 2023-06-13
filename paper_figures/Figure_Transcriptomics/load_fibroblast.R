#' The following loads RNA-seq data for patient derived fibroblast samples.
#' 
#' Requires input datasets
#'   `data/050322_RNAseq`
#'   `data/reactome_pathway.Rdata`
#'   `scripts/Figure_scores/models/predictions.Rdata`
#' 
#' **Note:** RNA-seq read count datasets can be processed from raw fastq files 
#' using the `als_rnaseq_fibro` pipeline.
#' Karl Kumbier 10/19/2022library(data.table)
library(tidyverse)
library(Matrix)
library(clusterProfiler)

options(stringsAsFactors=FALSE)

################################################################################
# Load data and initialize analysis paramters
################################################################################
# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.dir <- file.path(fig.base.dir, 'data/050322_RNAseq/')

pred.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Scores', 'fig_full') 
pathway.file <- str_c(data.dir, 'reactome_pathway.Rdata')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')
qc.file <- file.path(fig.base.dir, 'data', 'data Table 2.1 Sample Sequencing Statistics.csv')

load(file.path(data.dir, 'read_counts.Rdata'))
load(file.path(pred.dir, 'predictions.Rdata'))
xqc <- fread(qc.file) %>% dplyr::rename(CellLine=`Sample ID`)

# Source utility functions
source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
expression.threshold <- 2.5 # minimum log TPM + 1 level for gene filtering

################################################################################
# Initialize cell line metadata
################################################################################
# Read in cell line meta data
cell.lines <- str_remove_all(bam.files, '^.*/') %>% str_remove_all('Aligned.*$')

xcell <- fread(library.file) %>%
  dplyr::rename(CellLine=`Cell line`) %>%
  dplyr::rename(Genetics=`Disease Gene`)

xmeta <- data.frame(CellLine=cell.lines) %>% 
  left_join(xcell, by='CellLine') %>%
  mutate(Genetics=str_replace_all(Genetics, 'ALS with no known mutation', 'sporadic')) %>% 
  mutate(Genetics=str_replace_all(Genetics, 'WT', 'Healthy')) %>%
  mutate(Genetics=str_replace_all(Genetics, 'FUS', 'FUS-ALS')) %>%
  mutate(Sex=ifelse(Sex == 'm', 'M', Sex)) %>%
  mutate(Site=`Origin Institution`) %>%
  mutate(Institution=Site) %>%
  left_join(xqc, by='CellLine')

# Initialize cell line genetics key
gene.tab <- dplyr::select(xmeta, CellLine, Genetics)
gene.tab <- setNames(gene.tab$Genetics, gene.tab$CellLine)

# Merge metadata with imaging predictions
predictions <- group_by(predictions, CellLine) %>% 
  arrange(desc(Ypred)) %>%
  dplyr::select(CellLine, CellLinePassageNumber, Ypred)

xmeta <- left_join(xmeta, predictions, by='CellLine') %>%
  mutate(ALS=Genetics != 'Healthy')

#' # Preprocessing
#' Raw read counts are preprocessed as follows:
#' 
#' 1. Normalize expression within sample as transcript per million (TPM)
#' 2. Drop genes with low expression: average log(TPM + 1) < `r expression.threshold`
#' 3. Filter reactome pathway genes
#+ preprocessing

# Normalize read counts as transcripts per killobase million (TPM)
genes <- filter(genes, ensembl_gene_id %in% colnames(x)) %>%
  arrange(match(ensembl_gene_id, colnames(x)))

exon.length <- setNames(genes$ExonSize, genes$ensembl_gene_id)
genes.drop <- which(is.na(exon.length)) %>% names
x <- as.matrix(x[,!colnames(x) %in% genes.drop])
xcount <- x

x <- sapply(colnames(x), function(g) x[,g] / exon.length[g])
x <- t(apply(x, MAR=1, function(z) 1e6 * z / sum(z)))

# Drop low expression genes
xcount <- xcount[,colMedians(log(x + 1)) >= expression.threshold]
x <- x[,colMedians(log(x + 1)) >= expression.threshold]

# Filter to genes with annotated pathways
load(pathway.file)
xcount <- xcount[,colnames(x) %in% unlist(reactome.pathways$ENSEMBL)]
x <- x[,colnames(x) %in% unlist(reactome.pathways$ENSEMBL)]

# Rename for joint analysis
xcount.fibro <- xcount
x.fibro <- x
xmeta.fibro <- mutate(xmeta, Type='Fibroblast') 
rm(x, xmeta)