# " This script contains helper functions for running the predictive analyses.
# " The core function `fit_model` can be used with any supervised model by
# " passing in functions with standardized inputs/outputs for model fitting and
# " prediction. Currently, models for iterative random forests and logistic
# " regression with L1 penalty are supported.
# "
# " Note: irf is set to use all available cores.
# "
# " Karl Kumbier 10/4/2023
library(iRF)
library(glmnet)
library(keras)
library(caret)
library(PRROC)

cell_fit <- function(
    x, model, model_predict,
    genetics.train = "FUS-ALS", # nolint
    verbose = TRUE,
    ncell = 250) {
  # Wrapper function to call cell_fit_single predict over all cell lines.

  # Normalize features within plate
  ecdf <- function(x) rank(x) / length(x)
  xfeat <- ungroup(x) %>%
    dplyr::select(matches("(^X|PlateID)")) %>%
    group_by(PlateID) %>%
    mutate_if(is.numeric, ecdf) %>%
    ungroup() %>%
    dplyr::select(-PlateID)

  x <- dplyr::select(x, -matches("^X"))
  x <- cbind(x, xfeat)

  # Initialize grid for hold-out cell lines, training genetics
  cell.lines <- unique(x$CellLine)

  # Drop sporadic lines from model training — evaluated in all models
  sporadics <- filter(x, Genetics == "sporadic")$CellLine
  if (!"sporadic" %in% genetics.train) {
    cell.lines <- setdiff(cell.lines, sporadics)
  }

  # Fit model
  fits <- lapply(cell.lines, function(cl) {
    if (verbose) print(cl)
    cell_fit_single(
      x = x,
      cell.test = cl,
      genetics.train = genetics.train,
      model = model,
      model_predict = model_predict,
      ncell = ncell
    )
  })

  id.drop <- sapply(fits, is.null)
  fits <- fits[!id.drop]

  # Initialize importance
  importance <- sapply(fits, function(f) f$fit$fit$importance)
  ypred <- lapply(fits, function(f) f$ypred) %>% rbindlist()
  return(list(ypred = ypred, importance = importance))
}

cell_fit_single <- function(x, cell.test, model, model_predict,
                            genetics.train = "FUS-ALS", ncell = 250) {
  # Runs leave one cell line out prediction analysis.
  #
  # x: data frame containing predictive features (named as ^X.*), response
  #   (Genetics), CellLine, WellID, PlateID, and Markers.
  # cell.test: held-out cell line(s) to be evaluated.
  # model: standardized model fitting function. Takes arguments x, y only.
  # model_predict: standardized model prediction function. Takes arguments
  #   fit (fitted model) and x only
  # genetics.train: ALS genetic backgrounds to be used for model training.
  # ncell: number of cells per well to use in training

  # Check for valid input
  model.line <- cell.test
  print(str_c("Fitting ", cell.test, "..."))
  if (nrow(x) == 0) {
    return(NULL)
  }

  # Set training/test indices based on hold-out cell line
  sporadics <- filter(x, Genetics == "sporadic")$CellLine
  if (!"sporadic" %in% genetics.train) cell.test <- c(cell.test, unique(sporadics))
  cells.train <- setdiff(unique(x$CellLine), cell.test)

  # Set training/test indices based on hold-out cell line
  id.cell.train <- which(x$CellLine %in% cells.train)
  id.genetics <- which(x$Genetics %in% c(genetics.train, "Healthy"))

  id.train <- intersect(id.cell.train, id.genetics)
  id.test <- which(x$CellLine %in% cell.test)
  if (length(id.test) == 0) {
    return(NULL)
  }

  # Set feature matrix and response vector
  y <- as.numeric(x$Genetics != "Healthy")
  xx <- select(x, matches("^X"))

  # Initialize cell index table for down sampling
  xidx <- dplyr::select(x, ID) %>%
    mutate(Y = y, Idx = 1:n()) %>%
    filter(Idx %in% id.train)

  # Down sample for equal cell counts from each well
  xcount <- group_by(xidx, ID) %>% dplyr::count()
  sample.size <- min(c(xcount$n, ncell))
  xidx <- group_by(xidx, ID) %>% sample_n(sample.size)

  # Down sample for class balance
  xidx <- downSample(xidx, as.factor(xidx$Y))

  # Fit model
  id.train <- intersect(id.train, xidx$Idx)
  fit <- fit_model(xx, y, id.train, model, model_predict, test.id = id.test)

  # Initialize prediction table output
  f.select <- c("CellLine", "Genetics", "WellID", "PlateID", "Markers")
  out <- dplyr::select(x[id.test, ], one_of(f.select)) %>%
    mutate(Ypred = fit$ypred) %>%
    mutate(Ytest = as.numeric(Genetics != "Healthy")) %>%
    mutate(ModelLine = model.line) %>%
    mutate(TestIdx = id.test)

  return(list(fit = fit, ypred = out))
}

fit_model <- function(x, y, train.id, model, model_predict, test.id = NULL) {
  # Core function for model fitting. Fits the specified model, generates
  # predictions, and evaluates AUROC.
  #
  # x: input matrix of predictor features
  # y: input binary vector of response values
  # train.id: rows of x to be used for model training
  # model: standardized model fitting function. Takes arguments x, y only.
  # model_predict: standardized model prediction function. Takes arguments
  #   fit (fitted model) and x only
  # test.id: rows of x to generate predictions for. If NULL, all non training
  #   samples will be used.

  # Set test indices if not specified
  if (is.null(test.id)) test.id <- setdiff(seq_len(nrow(x)), train.id)

  # Fit model
  fit <- model(x = x[train.id, ], y = y[train.id])

  # Fit model and generate predictions
  ypred <- model_predict(fit, x[test.id, ])
  ytest <- y[test.id]
  roc <- roc.curve(scores.class0 = ypred, weights.class0 = ytest)

  return(list(
    fit = fit,
    auc = roc$auc,
    y = y,
    ypred = ypred,
    id.test = test.id
  ))
}


logistic <- function(x, y) {
  # Standardized function for L1 logistic regression.

  print("Fitting LR")
  fit <- cv.glmnet(
    x = as.matrix(x),
    y = y,
    family = "binomial",
    type.logistic = "modified.Newton",
    thresh = 1e-5
  )
  print("LR fit")

  fit$importance <- coef(fit, s = "lambda.min")[-1, 1]
  return(fit)
}


predict_logistic <- function(fit, x) {
  # Standardized function for logistic regression prediction.
  ypred <- predict(fit, as.matrix(x), s = "lambda.min", type = "response")
  return(ypred[, 1])
}

irf <- function(x, y, n.core = 1) {
  # Standardized function for irf.
  if (all(y %in% c(0, 1))) y <- as.factor(y)

  fit <- iRF(
    x = x,
    y = y,
    type = "ranger",
    verbose = FALSE,
    n.iter = 2,
    n.core = n.core
  )

  fit$importance <- fit$rf.list$variable.importance
  return(fit)
}

predict_irf <- function(fit, x) {
  # Standardized function for irf prediction.
  ypred <- predict(fit$rf.list, x, predict.all = TRUE)
  ypred <- rowMeans(ypred$prediction)
  return(ypred)
}

fcnn <- function(x, y) {
  # Standardized function for FC network

  if (all(y %in% 0:1)) {
    y <- to_categorical(y, 2)
    x <- as.matrix(x)
  } else {
    stop("Multiclass and regression not supported")
  }

  p <- c(ncol(x))

  model <- keras_model_sequential() %>%
    layer_dense(units = 16, activation = "relu", input_shape = p) %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 2, activation = "softmax")

  model %>% compile(
    loss = "categorical_crossentropy",
    optimizer = optimizer_rmsprop(),
    metrics = c("accuracy")
  )

  history <- model %>% fit(
    x, y,
    epochs = 30,
    batch_size = 128,
    validation_split = 0.2
  )


  # Evaluate feature importance
  ximportance <- diag(ncol(x))
  importance <- predict(model, ximportance)[, 2]
  importance <- abs(importance - mean(importance))
  names(importance) <- colnames(x)

  fit <- list()
  fit$model <- model
  fit$importance <- importance
  return(fit)
}

predict_fcnn <- function(fit, x) {
  return(predict(fit$model, as.matrix(x))[, 2])
}


###############################################################################
# Additional utility functions used in analysis
###############################################################################
label_pval <- function(x) {
  # Generate significance *s for p-values
  pval.lab.s <- " ns"
  pval.lab.s <- ifelse(x < 0.05, "*", pval.lab.s)
  pval.lab.s <- ifelse(x < 0.01, "**", pval.lab.s)
  pval.lab.s <- ifelse(x < 0.001, "***", pval.lab.s)
  return(pval.lab.s)
}

select_cytosolic <- function(x) {
  regions <- sapply(str_split(colnames(x), "\\.\\."), "[", 2)
  regions <- lapply(regions, str_split, pattern = "_")
  regions <- sapply(regions, function(z) tail(z[[1]], 1))
  regions[is.na(regions)] <- "X"

  # Filter to selected region / marker
  return(dplyr::select_if(x, regions %in% c("X", "N", "NN")))
}
