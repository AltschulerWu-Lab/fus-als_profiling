#' This script is used to generate eigenfeatures from raw imaging data.
#' Eigenfeatures are computed for a single screen / marker set for either
#' single cells or KS normalized well replicates.
#'
#' Requires input data
#'
#'  `<screen>`: raw imaging data from a screen with both
#'     `cbfeature_profiles.csv` and `metadata.csv` files. These can be generated
#'     using the `rprofiler` package in awlabtools (see `GENERATE_PROFILE` shell
#'     script in each screen subdirectory).
#'  `Fibroblasts_Library.csv`: cell line metadata
#'  `features_select_expanded.csv`: table of features to be included and their
#'     associated feature groups.
#'
#' Karl Kumbier 10/20/2022
library(tidyverse)
library(data.table)

als.dir <- Sys.getenv("ALS_PAPER")
source(file.path(als.dir, "preprocessing/utilities.R"))
args <- R.utils::commandArgs(asValues = TRUE)

# Variable parameters for eigenfeature generation
pct.var <- 0.95 # % variance explained by eigenfeatures within category
type <- "cell" # type pf sample, being 'cell' or 'well'
screen <- args$SCREEN #' 032422_WTFUSporadics' # screen to be analyzed
raw <- as.logical(args$RAW) # FALSE

# Directory containing cell line and feature metadata files
library.path <- file.path(als.dir, "data", "Fibroblasts_Library.csv")
feature.path <- file.path(als.dir, "data", "features_select_expanded.csv")

# Directories with raw imaging data, metadata, and eigenfeature output
base.dir <- "/awlab/projects/2017_01_ProjectALS/"
data.dir <- file.path(base.dir, "plates")
metadata.dir <- file.path(base.dir, "data", screen, "Cell_Filtered/")
output.dir <- file.path(base.dir, "data", screen, "Cell_Filtered/")

################################################################################
# Load data
################################################################################
# Load in metadata
xmeta <- fread(str_c(metadata.dir, "metadata.csv")) %>%
  mutate(Markers = str_remove_all(Markers, " ")) %>%
  mutate(Markers = str_remove_all(Markers, "Hoechst_")) %>%
  mutate(Markers = str_remove_all(Markers, "_Mitotracker"))

# Load in cell line data
xcell <- fread(library.path) %>%
  dplyr::rename(CellLine = `Cell line`) %>%
  dplyr::rename(Genetics = `Disease Gene`) %>%
  dplyr::select(CellLine, Genetics) %>%
  mutate(CellLine = str_remove_all(CellLine, " ")) %>%
  mutate(Genetics = ifelse(str_detect(Genetics, "mutation"), "sporadic", Genetics)) %>%
  mutate(Genetics = ifelse(Genetics == "WT", "Healthy", Genetics)) %>%
  mutate(Genetics = ifelse(Genetics == "FUS", "FUS-ALS", Genetics))

xmeta <- left_join(xmeta, xcell, by = "CellLine")

# Clean markerset names
markers <- unique(xmeta$Markers)

# Initialize feature grouping key
xmarker <- fread(feature.path) %>%
  rename(Feature = feature_name) %>%
  rename(FeatureAbbr = feature_short) %>%
  rename(Channel = channel) %>%
  rename(Group = group) %>%
  distinct()

# Generate key for mapping feature full names to feature short names
feature.key <- select(xmarker, Feature, FeatureAbbr) %>% distinct()

# Drop channels for select markersets
lapply(markers, function(m) {
  print(m)
  xmeta.m <- dplyr::filter(xmeta, Markers == m)
  plate.ids <- unique(xmeta.m$PlateID)

  xm <- lapply(plate.ids, function(p) {
    data.file <- str_c("cbfeature_", type, "_profiles.csv")
    data.path <- file.path(data.dir, p, "ks_profiles", data.file)
    fread(data.path)
  })

  # Merge data / metadata
  xm <- rbindlist(xm) %>%
    mutate(ID = str_c(PlateID, "_", ID)) %>%
    dplyr::select(-PlateID) %>%
    left_join(xmeta.m, by = "ID")

  if (raw) {
    features <- str_replace_all(xmarker$Feature, "\\*", "\\.\\*")
    features <- str_c("X\\.\\.", features) %>% str_c(collapse = "|")
    xf <- dplyr::select(xm, matches(str_c("(", features, ")")))

    channels <- str_remove_all(colnames(xf), "^.*\\.\\.") %>% str_split("")
    id.channel <- sapply(channels, function(ch) all(ch %in% 2:3))
    xf <- dplyr::select_if(xf, id.channel)
    xf.meta <- dplyr::select(xm, -matches("^X"))
    x <- cbind(xf.meta, xf) %>% mutate_if(is.numeric, median_impute)
    save(file = file.path(output.dir, "profiles_raw.Rdata"), x)
  } else {
    # Reduce dataset based on select features
    xm <- subset_features(xm, m, xmarker$Channel, xmarker$Feature, pct.var)

    # Update column names with feature abbrevaiation
    col.names <- str_split(colnames(xm), "\\.\\.")
    col.names <- sapply(col.names, function(z) {
      if (length(z) == 1) {
        return(z)
      }
      z.feat <- str_remove_all(z[2], "\\.")
      z.feat <- feature.key$FeatureAbbr[feature.key$Feature == z.feat]
      out <- str_c(z[1], "..", z.feat, "..", z[3])
      return(str_replace_all(out, ":PC", ".."))
    })

    colnames(xm) <- col.names
    x <- xm

    output.path <- str_c(output.dir, "profiles_", m, ".Rdata")
    save(file = output.path, x)
  }
})
