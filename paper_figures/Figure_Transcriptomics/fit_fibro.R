################################################################################
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
################################################################################
library(data.table)
library(tidyverse)
library(Matrix)
library(mltools)
library(ggrepel)
library(parallel)

options(stringsAsFactors = FALSE)
theme_set(theme_bw(base_size = 22))

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.dir <- file.path(fig.base.dir, "data/050322_RNAseq/")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Transcriptomics")
library.file <- file.path(fig.base.dir, "data", "Fibroblasts_Library.csv")
pathway.file <- file.path(data.dir, "reactome_pathway.Rdata")

output.fig.dir <- file.path(figure.dir, "fig")
dir.create(output.fig.dir, showWarnings = FALSE)
source(file.path(figure.dir, "utilities.R"))

# Parameters for analysis
n.core <- 12
nrep <- 10
bootstrap <- TRUE

gene_transform <- function(x) rank(x)
model <- fit_lm_classification
model_predict <- predict_lm_classification

################################################################################
#' ## Data preprocessing
#+ load_data, echo=FALSE, warning=FALSE, message=FALSE
################################################################################
load(file.path(data.dir, "fibro_processed.Rdata"))
load(file.path(output.fig.dir, "spine_models.Rdata"))
load(file.path(output.fig.dir, "genes.Rdata"))
load(pathway.file)

colnames(xmeta.fibro) <- str_replace_all(colnames(xmeta.fibro), " ", "_")

x <- x.fibro[, shared.genes]
xmeta <- mutate(xmeta.fibro, Y = as.numeric(ALS))
rm(x.fibro, xmeta.fibro)

# Initialize one-hot encoding of metadata features
xonehot <- dplyr::select(xmeta, Site, Sex) %>%
  mutate_all(as.factor) %>%
  as.data.table() %>%
  mltools::one_hot()

# Normalize features for modeling
x <- apply(x, MAR = 2, gene_transform)

# Load in cell line library metadata
xcell <- fread(library.file) %>%
  dplyr::rename(CellLine = `Cell line`) %>%
  dplyr::rename(Fig = `Figure Names`)

################################################################################
#' ## Modeling
#' Transcription model of ALS. We fit supervised models trained to predict
#' disease status. We use a leave-one-out sample split to evaluate predictions
#' on held-out WT & FUS cell lines. Genes used as features in classification
#' models are taken from pathways enriched in patient-derived spinal cord
#' samples.
#+ models, fig.height=8, fig.width=12, cache=TRUE
################################################################################
pw.genes <- pw.enriched$ENSEMBL
pathways <- pw.enriched$Name

fit <- mclapply(pw.genes, function(g) {
  fit_genes(x, xmeta, xonehot, g, model, model_predict, nrep, bootstrap = bootstrap)
}, mc.cores = n.core)

ypred <- lapply(1:length(fit), function(i) {
  return(mutate(fit[[i]]$Ypred, Name = pathways[i]))
})

ypred <- rbindlist(ypred) %>%
  group_by(CellLine, Genetics, Name, Site, Sex) %>%
  summarize(YpredSeq = mean(YpredSeq), Ypred = mean(Ypred), .groups = "drop") %>%
  mutate(ALS = as.numeric(Genetics != "Healthy"))

coefs <- lapply(1:length(fit), function(i) {
  coef.i <- reshape2::melt(fit[[i]]$Coef) %>%
    mutate(Name = pathways[i]) %>%
    dplyr::rename(CellIdx = L1) %>%
    filter(Var1 != "Intercept") %>%
    mutate(AbsBeta = abs(value)) %>%
    mutate(SignBeta = sign(value)) %>%
    dplyr::select(Name, Var1, CellIdx, AbsBeta, SignBeta)

  return(coef.i)
})

fout <- file.path(output.fig.dir, "fibro_models.Rdata")
save(file = fout, fit, ypred, coefs)
