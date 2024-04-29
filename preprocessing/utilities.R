library(tidyverse)
library(data.table)

median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

subset_features <- function(x, marker, channels, features, pct.var = 0) {
  # Filter dataset to selected markers/channels/features

  # Filter to selected marker set
  xm <- filter(x, str_detect(Markers, marker)) %>%
    select_if(function(col) !all(is.na(col)))

  # Filter to select channels
  channels.select <- unlist(str_split(channels, ""))
  c.drop <- str_c(setdiff(1:4, channels.select), collapse = "|")
  f.drop <- str_c("X\\.\\..*\\.\\..*(", c.drop, ").*$")
  id.drop <- str_detect(colnames(xm), f.drop)
  xm <- select_if(xm, !id.drop)

  # Filter to select features
  features <- str_replace_all(features, "\\*", "\\.\\*")
  features.select <- str_c("X..", features, "..", channels)
  xf <- lapply(features.select, function(f) select(xm, matches(f)))

  # Reduce dimensions of filtered data
  if (pct.var != 0) {
    xf <- lapply(xf, reduce_dim, pct.var = pct.var)
    p <- sapply(xf, ncol)
    features.select <- rep(features.select, times = p)
    xf <- do.call(cbind, xf)
    colnames(xf) <- str_c(features.select, ":", colnames(xf))
  } else if (pct.var == "average") {
    xf <- sapply(xf, rowMeans)
    colnames(xf) <- features.select
  } else {
    xf <- do.call(cbind, xf)
    xf <- mutate_all(xf, median_impute)
  }

  # Aggregatedata and metadata
  meta.cols <- str_subset(colnames(xm), "^X", negate = TRUE)
  xm <- select(xm, one_of(meta.cols))
  out <- cbind(xm, xf)
  return(out)
}

reduce_dim <- function(x, pct.var = 0.9) {
  # Reduce dataset dimensionality with PCA
  x <- mutate_all(x, median_impute)
  pca <- prcomp(x)

  # Evaluate # PCs required to obtain specified level of explained variability
  var.exp <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n.pc <- min(which(var.exp > pct.var))
  xout <- pca$x[, 1:n.pc, drop = FALSE]
  return(xout)
}
