library(DESeq2)
library(fgsea)

################################################################################
# Differential expression and gene set enrichment analysis
################################################################################
deseq_design <- function(x, xmeta, design, ...) {
  # Wrapper function for running deseq over targeted samples
  deseq <- DESeqDataSetFromMatrix(
    countData=t(round(x)),
    colData=xmeta,
    design=as.formula(design),
  )
  
  deseq <- DESeq(deseq, ...)
  return(deseq)
}


gsea <- function(enrichment, pathway.file, pthr=0.01, ...) {
  # Compute gene set enrichment for selected pathway

  load(pathway.file)
  pathways <- reactome.pathways$Genes
  names(pathways) <- reactome.pathways$Name

  entrez <- clusterProfiler::bitr(
    names(enrichment),
    fromType="ENSEMBL",
    toType="ENTREZID",
    OrgDb='org.Hs.eg.db'
  )
  
  # Join entrez ids with enrichment data
  gene.set <- data.frame(Loading=unname(enrichment)) %>%
    mutate(ENSEMBL=names(enrichment)) %>%
    right_join(entrez, by='ENSEMBL')
  
  # Run gene set enrichment analysis for bag data
  gene.list <- gene.set$Loading
  names(gene.list) <- gene.set$ENTREZID
  
  out <- fgsea(
    pathways=pathways,
    stats=gene.list,
    scoreType='pos',
    eps=1e-30,
    ...
  )
  
  out <- filter(out, padj < pthr) %>% arrange(padj)
  out.collapse <- collapsePathways(out, pathways, gene.list)
  out.f <- filter(out, pathway %in% out.collapse$mainPathways)
  return(list(pw=out.f, filtered=out.collapse, pw.full=out))
}


jaccard <- function(x, y) {
  # Compute jaccard distance between two sets
  return(length(intersect(x, y)) / length(union(x, y)))
}

pathway_distance <- function(gene.set, pw.names=NULL) {
  # Compute jaccard distance matrix over set of pathways
  pairs <- combn(1:length(gene.set), 2, simplify=FALSE)
  dvec <- sapply(pairs, function(p) jaccard(gene.set[[p[1]]], gene.set[[p[2]]]))
  
  # Initialize pathway-level distance matrix
  dmat <- matrix(0, nrow=length(gene.set), ncol=length(gene.set))
  
  dmat[lower.tri(dmat)] <- dvec
  dmat <- dmat + t(dmat)
  diag(dmat) <- 1
  
  if (!is.null(pw.names)) {
    rownames(dmat) <- pw.names
    colnames(dmat) <- pw.names
  }
  
  return(dmat)
}

subsample_deseq <- function(x, xmeta, 
                            designs='~Genetics + Site',
                            name='Genetics_Healthy_vs_FUS.ALS') {
  # Wrapper function for computing deseq over subsamples samples
  
  id <- sample(nrow(x), 0.75 * nrow(x))
  
  out <- sapply(designs, function(d) {
    deseq <- deseq_design(
      x[id,], 
      xmeta[id,],
      design=d,
      fitType='mean'
    )
    
    res <- results(deseq, name=name)
    return(setNames(-log10(res$pvalue), rownames(res)))
  })
  
  return(out)
}

################################################################################
# Supervised analysis
################################################################################
fit_holdout <- function(x, xmeta, lines.predict, fit_, predict_, importance_) {
  # Wrapper function for supervised analysis with holdout cell lines
  train.id <- !xmeta$CellLine %in% lines.predict
  test.id <- xmeta$CellLine %in% lines.predict
  y <- as.numeric(xmeta$Y)
  labels <- as.numeric(xmeta$ALS)
  
  
  x <- mutate(data.frame(x), Y=xmeta$Y)
  xus <- caret::upSample(x[train.id,], as.factor(labels[train.id]))
  xi <- as.matrix(dplyr::select(xus, -Class, -Y))
  yi <- xus$Y
  fit <- fit_(x=xi, y=yi)
  
  x <- dplyr::select(x, -Y) %>% as.matrix
  ypred <- predict_(fit, x[test.id,])
  ypred <- data.frame(CellLine=xmeta$CellLine[test.id]) %>%
    mutate(YpredSeq=ypred[,1]) %>%
    mutate(Genetics=xmeta$Genetics[test.id])
  
  
  return(ypred)
}


fit_model <- function(x, y, id.train, model, model_predict) {
  # Wrapper function for fitting model and predicting on held-out samples
  fit <- model(x[id.train,], y[id.train])
  ypred <- model_predict(fit, x[!id.train,])
  return(list(fit=fit, ypred=ypred, ytest=y[!id.train]))
}

################################################################################
# Wrapper functions for running supervised analysis
################################################################################
fit_lm_classification <- function(x, y, ...) {
  # Standardized function for sparse logistic regression
  return(cv.glmnet(x=x, y=y, nfolds=3, family='binomial', intercept=FALSE, ...))
}

predict_lm_classification <- function(fit, x) {
  # Standardized function for sparse logistic regression predictions
  lambda.min <- fit$lambda.min
  
  if (which(fit$lambda == lambda.min) == 1) {
    ypred <- predict(fit$glmnet.fit, x, type='response')
    null.model <- colMeans(ypred == 0.5) == 1
    ypred <- ypred[,min(which(!null.model))] %>% as.matrix
  } else {
    ypred <- predict(fit, x, s='lambda.min', type='response')
  }
    
  
  return(ypred)
}
