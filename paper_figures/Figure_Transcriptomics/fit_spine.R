###############################################################################
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
################################################################################
library(data.table)
library(tidyverse)
library(Matrix)
library(caret)
library(iRF)
library(PRROC)
library(parallel)

options(stringsAsFactors = FALSE)
theme_set(theme_bw(base_size = 22))

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.dir <- file.path(fig.base.dir, "data/050322_RNAseq/")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Transcriptomics")
output.fig.dir <- file.path(figure.dir, "fig")

n.core <- 12
train.prop <- 0.5 # training proportion
min.healthy <- 3 # minimum number of healthy samples to include site in analysis
n.null <- 100000 # number of bootstrap subsamples for null distribution
min.size <- 5 # minimum pathway size included in analysis
max.size <- 500 # maximmum pathway size included in analysis

source(file.path(figure.dir, "utilities.R"))
load(file.path(data.dir, "spine_processed.Rdata"))
load(file.path(output.fig.dir, "genes.Rdata"))
load(file.path(data.dir, "reactome_pathway.Rdata"))

auc <- function(x) {
  # Compute AUROC from data frame with columns `Ypred`, `Ytrue`
  roc.curve(scores.class0 = x$Ypred, weights.class0 = x$Ytrue)$auc
}

pw_predict <- function(x, y, id.train, genes) {
  # Fit ALS vs. health classification models on subset of genes
  fselect <- str_c(genes, collapse = "|")
  fselect <- str_c("(", fselect, "|Sex|^Site$|^Age$|pct_)")
  x <- dplyr::select(x, matches(fselect))
  return(fit_sample_split(x, y, id.train))
}

################################################################################
#' ## Load data
#' Below Initialize spinal cord data sets used in classification models. We
#' include sample metadata features: `Sex`, `Site` (of sample collection),
#' `Age`, and `pct_x` (x corresponding to ancestry computed from WGS data).
#+ modeling
################################################################################
xmeta <- dplyr::select(xmeta.spine, matches("(^Sex$|^Site$|^Age$|pct_)"))
x.spine <- x.spine[, shared.genes]

x <- cbind(x.spine, xmeta)
y <- as.factor(xmeta.spine$ALS)

# Filter to sites with both ALS and healthy control samples
sites <- group_by(xmeta.spine, Site) %>% summarize(NHealthy = sum(!ALS))
dplyr::select(xmeta.spine, Site, ALS) %>% table()

genes <- colnames(x.spine)
id.drop <- xmeta$Site %in% sites$Site[sites$NHealthy < min.healthy]

x <- x[!id.drop, ]
y <- y[!id.drop]
xmeta <- xmeta[!id.drop, ]
x.spine <- x.spine[!id.drop, ]

# Take balanced sample within site
set.seed(47)

xtrain <- data.frame(x, Y = y) %>%
  mutate(Idx = 1:n()) %>%
  group_by(Site, Y) %>%
  sample_frac(train.prop)

id.train <- xtrain$Idx
id.test <- setdiff(1:nrow(x), id.train)
ytrue <- as.numeric(y) - 1

################################################################################
#' ## Null model
#' Performance of a null model is estimated
#+ null_model
################################################################################
bs.train <- replicate(n.null, sample(id.train, replace = TRUE), simplify = FALSE)
bs.test <- replicate(n.null, sample(id.test, replace = TRUE), simplify = FALSE)

fit.null <- mcmapply(function(trn, tst) {
  fit <- fit_sample_split(xmeta, y, trn, tst)
  ypred <- data.frame(Ypred = fit$ypred, Ytrue = ytrue[fit$id])
  acc <- data.frame(AUC = auc(ypred))
  return(list(ypred = ypred, acc = acc))
}, bs.train, bs.test, mc.cores = n.core, SIMPLIFY = FALSE)

null.accuracy <- rbindlist(lapply(fit.null, function(z) z$acc))
null.pred <- lapply(fit.null, function(z) z$ypred)

################################################################################
#' ## Pathway models
#' For each reactome pathway, we fit ALS vs. healthy classifiers using the
#' subset of genes associated with that pathway.
#+ pw_models
################################################################################
reactome.pathways <- reactome.pathways %>%
  mutate(ENSEMBL = sapply(ENSEMBL, intersect, y = genes)) %>%
  mutate(Size = sapply(ENSEMBL, length)) %>%
  filter(Size > min.size & Size <= max.size) %>%
  dplyr::select(-Genes, -ID)

# Fit classification models for each reactome pathway
pw.models <- mclapply(reactome.pathways$Name, function(pw) {
  genes <- filter(reactome.pathways, Name == pw)$ENSEMBL
  genes <- na.omit(unlist(genes))
  fit.g <- pw_predict(x, y, id.train, genes)
  ypred <- data.frame(Ypred = fit.g$ypred, Ytrue = ytrue[fit.g$id])

  return(list(
    acc = data.frame(AUC = auc(ypred), Name = pw),
    ypred = ypred,
    importance = fit.g$importance
  ))
}, mc.cores = n.core)

pw.aucs <- rbindlist(lapply(pw.models, function(z) z$acc)) %>%
  arrange(desc(AUC)) %>%
  left_join(reactome.pathways, by = "Name") %>%
  mutate(p = sapply(AUC, function(z) mean(null.accuracy$AUC > z)))

fout <- file.path(output.fig.dir, "spine_models.Rdata")
save(file = fout, pw.aucs, pw.models, null.pred, null.accuracy)
