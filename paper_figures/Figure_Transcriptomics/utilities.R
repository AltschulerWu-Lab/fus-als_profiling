library(glmnet)
library(iRF)

################################################################################
# Supervised analysis
################################################################################
fit_holdout <- function(x, xmeta, lines.predict, fit_, predict_, bootstrap=FALSE) {
  # Wrapper function for supervised analysis with holdout cell lines
  train.id <- !xmeta$CellLine %in% lines.predict
  test.id <- xmeta$CellLine %in% lines.predict
  
  x <- mutate(data.frame(x), Y=xmeta$Y)
  xus <- caret::upSample(x[train.id,], as.factor(xmeta$ALS[train.id]))
  
  if (bootstrap) {
    xus <- group_by(xus, Class) %>% 
      sample_frac(1, replace=TRUE) %>%
      ungroup()
  }
  
  xi <- as.matrix(dplyr::select(xus, -Class, -Y))
  fit <- fit_(x=as.matrix(xi), y=xus$Y)
  
  # Get model coefficients
  betas <- coef(fit, s='lambda.min')
  if (sum(betas) == 0) {
    betas <- coef(fit$glmnet.fit)
    betas.sum <- Matrix::colSums(betas)
    betas <- betas[,min(which(betas.sum != 0))]
  }
  
  x <- dplyr::select(x, -Y) %>% as.matrix
  ypred <- predict_(fit, x[test.id,])
  ypred <- mutate(xmeta[test.id,], YpredSeq=ypred[,1])
  return(list(Ypred=ypred, Coef=betas))
}

fit_genes <- function(x, xmeta, xonehot, genes, fit_, predict_, 
                      nrep=10, bootstrap=FALSE) {
  
  # Runs leave-one out prediction analysis, fitting models on selected subset of 
  # full gene set.
  xg <- x[,intersect(genes, colnames(x))]
  xg <- cbind(xg, xonehot)
  
  # Initialize subset of FUS/WT cell lines for binary classification
  cell.lines <- filter(xmeta, Genetics != 'sporadic')$CellLine
  s.cell.lines <- filter(xmeta, Genetics == 'sporadic')$CellLine

  out <- lapply(1:length(cell.lines), function(i) {
    set.seed(i)
    lines.predict <- c(s.cell.lines, cell.lines[i])
    
    out <- replicate(nrep, {
      fit_holdout(xg, xmeta, lines.predict, fit_, predict_, bootstrap)
    }, simplify=FALSE)
    
    ypred <- lapply(out, function(z) z$Ypred) %>% rbindlist
    coefs <- sapply(out, function(z) as.matrix(z$Coef))
    rownames(coefs) <- c('Intercept', colnames(xg))
    return(list(Ypred=ypred, Coef=coefs))
  })
  
  ypred <- lapply(out, function(z) z$Ypred) %>% rbindlist
  coefs <- lapply(out, function(z) z$Coef)
  return(list(Ypred=ypred, Coef=coefs))
}


fit_sample_split <- function(x, y, id.train, id.test=NULL, n.iter=1, n.core=1) {
  # Fits iRF and predicts RF models on sample split
  if (is.null(id.test)) id.test <- setdiff(1:nrow(x), id.train)
  
  fit.trn <- iRF(
    x=x[id.train,],
    y=y[id.train],
    type='ranger',
    n.iter=n.iter,
    n.core=n.core,
    respect.unordered.factors=TRUE
  )
  
  fit.tst <-  iRF(
    x=x[id.test,],
    y=y[id.test],
    type='ranger',
    n.iter=n.iter,
    n.core=n.core,
    respect.unordered.factors=TRUE
  )
  
  ypred.train <- predict(fit.tst$rf.list, x[id.train], predict.all=TRUE)
  ypred.train <- rowMeans(ypred.train$predictions)
  
  ypred.test <- predict(fit.trn$rf.list, x[id.test], predict.all=TRUE)
  ypred.test <- rowMeans(ypred.test$predictions)
  
  ypred <- c(ypred.train, ypred.test)
  id.split <- c(id.train, id.test)
  
  importance.trn <- fit.trn$rf.list$variable.importance
  importance.tst <- fit.tst$rf.list$variable.importance
  importance <- (importance.trn + importance.tst) / 2
  
  return(list(id=id.split, ypred=ypred, importance=importance))
}

################################################################################
# Wrapper functions for running supervised analysis
################################################################################
fit_lm <- function(x, y, ...) {
  # Standardized function for sparse logistic regression
  return(cv.glmnet(x=x, y=y, nfolds=3, intercept=FALSE, ...))
}

predict_lm <- function(fit, x) {
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

################################################################################
# General utility functions
################################################################################
reweight_importance <- function(importance, xcor) {
  # Computes correlation-weighted feature importance. 
  importance <- importance[rownames(xcor)]
  out <- abs(xcor) %*% as.matrix(importance)
  return(out[,1])
}

get_children <- function(nodes, edges, ids, children=list()) {
  # Traverse graph to select children of selected node
  if (!any(ids %in% edges$from)) return(children)
  
  new.children <- edges$to[edges$from %in% ids]
  children <- c(children, list(new.children))
  return(get_children(nodes, edges, new.children, children))
}

set_node_gradient <- function(x, col.pal=viridis(100), col.range=NULL) {
  if (all(is.na(x))) return(rep('#707579', length(x)))

  if (is.null(col.range)) {
    col.range <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  }

  if (length(unique(col.range)) == 1) stop('Color range must contain unique levels')
  cut.seq <- seq(col.range[1], col.range[2], length.out=length(col.pal))
  out <- rev(col.pal)[cut(x, cut.seq, include.lowest=TRUE)]
  out[is.na(out)] <- '#707579'
  return(out)
}

