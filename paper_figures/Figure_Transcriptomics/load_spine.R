################################################################################
#+ setup
################################################################################
library(Matrix)
library(tidyverse)
library(data.table)
library(DESeq2)

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.dir <- file.path(fig.base.dir, "data/050322_RNAseq/")
nygc.dir <- file.path(data.dir, "nygc/")

pathway.file <- str_c(data.dir, "reactome_pathway.Rdata")
library.file <- file.path(fig.base.dir, "data", "Fibroblasts_Library.csv")

# Source utility functions
source(file.path(fig.base.dir, "paper_figures", "utilities.R"))
expression.threshold <- 2.5 # minimum log read count threshold

# Load exon length table
load(file.path(data.dir, "read_counts.Rdata"))
gene.length <- dplyr::select(genes, ensembl_gene_id, ExonSize)

# Load spinal cord read count data
load(file.path(nygc.dir, "read_counts_aw.Rdata"))
genes <- left_join(genes, gene.length, by = "ensembl_gene_id")

# Initialize patient populations used for analysis
als.groups <- "ALS Spectrum MND"
control.groups <- "Non-Neurological Control"
site <- "Spinal_Cord_Cervical"

# Initialize metadata table
meta.file <- "RNA_Metadata_200212/Delivered Quotes-Table 1.csv"
xmeta <- fread(file.path(nygc.dir, meta.file))

colnames(xmeta) <- colnames(xmeta) %>%
  str_remove_all(" ") %>%
  str_remove_all("\\(.*$")

xmeta <- dplyr::rename(xmeta, Mutation = ReportedGenomicMutations) %>%
  mutate(SubjectGroup = str_replace_all(SubjectGroup, "DIs", "Dis"))

# Filter to metadata with RNA-seq samples
xmeta <- filter(xmeta, ExternalSampleId %in% xbam$ExternalSampleId)

#' # Preprocessing
#'
#' We preprocess data as follows:
#' 1. Filter to samples derived from patient spinal cords
#' 2. Filter to ALS and non-neurological control samples
#' 3. Normalize expression within sample as transcript per million (TPM)
#' 4. Drop genes with low expression: average log(TPM + 1) < `r expression.threshold`
#' 5. Filter reactome pathway genes
#+ preprocessing, fig.height=8, fig.width=12

# Filter to patient spinal cord samples
id.spine <- xmeta$SampleSource %in% site
id.groups <- xmeta$SubjectGroup %in% c(als.groups, control.groups)

x <- x[id.spine & id.groups, ]
xmeta <- xmeta[id.spine & id.groups, ]

# Normalize read counts as transcripts per killobase million (TPM)
genes <- filter(genes, ensembl_gene_id %in% colnames(x)) %>%
  arrange(match(ensembl_gene_id, colnames(x)))

exon.length <- setNames(genes$ExonSize, genes$ensembl_gene_id)

genes.drop <- which(is.na(exon.length)) %>% names()
x <- as.matrix(x[, !colnames(x) %in% genes.drop])
xcount <- x

x <- sapply(colnames(x), function(g) x[, g] / exon.length[g])
x <- t(apply(x, MAR = 1, function(z) 1e6 * z / sum(z)))

# Drop low expression genes
xcount <- xcount[, colMedians(log(x + 1)) >= expression.threshold]
x <- x[, colMedians(log(x + 1)) >= expression.threshold]

# Filter to genes with annotated pathways
load(pathway.file)
xcount <- xcount[, colnames(x) %in% unlist(reactome.pathways$ENSEMBL)]
x <- x[, colnames(x) %in% unlist(reactome.pathways$ENSEMBL)]

# Standardize metadata feature names
xmeta <- mutate(xmeta, ALS = SubjectGroup %in% als.groups) %>%
  dplyr::rename(Site = SiteSpecimenCollected) %>%
  dplyr::rename(Age = AgeatDeath) %>%
  mutate(Age = ifelse(Age == "90 or Older", 90, Age)) %>%
  mutate(Age = as.numeric(Age))

# Rename for joint analysis
x.spine <- x
xcount.spine <- xcount
xmeta.spine <- mutate(xmeta, Type = "Spine")

fout <- file.path(data.dir, "spine_processed.Rdata")
save(file = fout, x.spine, xcount.spine, xmeta.spine)
