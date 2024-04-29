library(tidyverse)
library(data.table)
library(readxl)

setwd("~/github/projectALS/")
source("121820_SporadicScreen/utilities.R")
output.dir <- "data/profiles/121820_SporadicScreen_MergedWell_Condensed/"
dir.create(output.dir, showWarnings = FALSE)

################################################################################
# Data loading
################################################################################
# Set paths to datasets used for integrated analysis
data.path.1 <- "data/profiles/110420_BiomarkerScreen/cbfeature_profiles.csv"
data.path.2 <- "data/profiles/121820_SporadicScreen/cbfeature_profiles.csv"

meta.path.1 <- "data/profiles/110420_BiomarkerScreen/metadata.csv"
meta.path.2 <- "data/profiles/121820_SporadicScreen/metadata.csv"

desc.path.1 <- "data/profiles/110420_BiomarkerScreen/feature_descriptors.csv"
desc.path.2 <- "data/profiles/121820_SporadicScreen/feature_descriptors.csv"
feature.path <- "data/features_select_condensed.csv"

# Load in data from 11/20 and 12/20 screens
x1 <- fread(data.path.1) %>%
  dplyr::mutate(ID = str_c(PlateID, ID, sep = "_")) %>%
  dplyr::select(-PlateID) %>%
  dplyr::filter(!Bootstrap)

x2 <- fread(data.path.2) %>%
  dplyr::mutate(ID = str_c(PlateID, ID, sep = "_")) %>%
  dplyr::select(-PlateID) %>%
  dplyr::filter(!Bootstrap)

# Load in metadata from 11/20 and 12/20 screens
xmeta1 <- fread(meta.path.1)
x1 <- left_join(x1, xmeta1, by = "ID")
x1 <- as.data.frame(x1) %>% mutate(Screen = "1120")

# Match marker set names
x1$Markers <- str_replace_all(x1$Markers, "_EEA1 _CD63", "_CD63 _EEA1")
x1$Markers <- str_replace_all(x1$Markers, "p-p62", "phospho-p62")

xmeta2 <- fread(meta.path.2)
x2 <- left_join(x2, xmeta2, by = "ID")
x2 <- as.data.frame(x2) %>% mutate(Screen = "1220")

# Merge data and metadata
x <- rbindlist(list(x1, x2), fill = TRUE)

# Update marker set names
x <- mutate(x, Markers = str_remove_all(Markers, " "))
markers <- unique(x$Markers)

# Update compound names
x <- mutate(x, Compounds = str_replace_all(Compounds, "None", "Control")) %>%
  mutate(Compounds = str_replace_all(Compounds, "0.1%_DMSO", "Control")) %>%
  mutate(CellLine = str_remove_all(CellLine, " ")) %>%
  mutate(CellLine = str_replace_all(CellLine, "308", "ALS308")) %>%
  mutate(CellLine = str_replace_all(CellLine, "309", "ALS309"))

# Load in cell line data
xcl <- fread("data/2020_Fibroblasts_Library.csv") %>%
  dplyr::rename(CellLine = `Cell line`) %>%
  dplyr::rename(CellLineGenetics = `Disease Gene`) %>%
  dplyr::select(CellLine, CellLineGenetics) %>%
  mutate(CellLine = str_remove_all(CellLine, " "))

x <- left_join(x, xcl, by = "CellLine")

# Load in feature desc file
xdesc1 <- fread(desc.path.1) %>%
  mutate(Full = str_c("X..", Feature, "..", Channel)) %>%
  mutate(Full = str_replace_all(Full, "-", "."))

xdesc2 <- fread(desc.path.2) %>%
  mutate(Full = str_c("X..", Feature, "..", Channel)) %>%
  mutate(Full = str_replace_all(Full, "-", "."))

xdesc <- rbindlist(list(xdesc1, xdesc2)) %>% distinct()

# Load in select features file
xmarker <- fread(feature.path) %>%
  rename(Feature = feature_name) %>%
  rename(FeatureAbbr = feature_short) %>%
  rename(Channel = channel) %>%
  rename(Group = group) %>%
  distinct() %>%
  select(-V5, -V6)

################################################################################
# Data preprocessing
################################################################################
# Markers for which to generate datasets
markers <- c(
  "CD63_EEA1",
  "EEA1",
  "CD63",
  "p62_phospho-p62S403",
  "p62",
  "phospho-p62S403",
  "LAMP1_TFEB",
  "LAMP1",
  "TFEB",
  "FUS",
  "TIAR"
)

# Filter to WT cell lines shared across screens
x.lines <- filter(x, CellLineGenetics == "WT") %>%
  group_by(Screen) %>%
  summarize(CellLine = list(unique(CellLine)), .groups = "drop")

# Normalize data by median shift over shared WT cell lines
shared.lines <- reduce(x.lines$CellLine, intersect)

# Generate key for mapping feature full names to feature short names
feature.key <- select(xmarker, Feature, FeatureAbbr) %>% distinct()

# Process each markerset independently and return as list
lapply(markers, function(m) {
  # Drop channels for single marker evaluation
  if (m %in% c("FUS", "EEA1", "phospho-p62S403", "TFEB")) {
    xmarker <- filter(xmarker, Channel == 3)
  } else if (m %in% c("TIAR", "p62", "CD63", "LAMP1")) {
    xmarker <- filter(xmarker, Channel == 2)
  }

  print(str_c("Processing markerset: ", m, "..."))
  xm <- subset_features(x, m, xmarker$Channel, xmarker$Feature, 0.95)
  xm <- batch_normalize(as.data.frame(xm), shared.lines)

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

  output.path <- str_c(output.dir, "cbfeature_profiles_", m, ".Rdata")
  save(file = output.path, x)
})
