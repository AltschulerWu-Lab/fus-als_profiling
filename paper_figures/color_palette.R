library(tidyverse)
library(RColorBrewer)
library(data.table)

genetics.include <- c('C9orf72', 'FUS-ALS', 'sporadic', 'TBK1', 'Healthy')

library.file <- 'Fibroblasts_Library.csv'
library.path <- file.path(Sys.getenv('ALS_PAPER'), 'data', library.file)

# Load in cell line metadata
xcell <- fread(library.path) %>% 
  dplyr::rename(Genetics=`Disease Gene`, CellLine=`Cell line`) %>%
  mutate(Genetics=ifelse(Genetics == 'WT', 'Healthy', Genetics)) %>%
  mutate(Genetics=ifelse(Genetics == 'FUS', 'FUS-ALS', Genetics)) %>%
  mutate(Genetics=str_replace_all(Genetics, 'ALS with no known mutation', 'sporadic'))

# Initialize cell line/age key
age.key <- setNames(xcell$`Age of biopsy`, xcell$CellLine)

# Initialize cell line/genetics key
xg <- dplyr::select(xcell, CellLine, Genetics) %>% 
  distinct() %>%
  group_by(Genetics) %>% 
  summarize(CellLine=list(unique(CellLine))) %>%
  filter(Genetics %in% genetics.include)

# Initialize base color palette
col.pals <- list(
  brewer.pal(8, 'Purples')[-(1:3)],
  brewer.pal(8, 'Reds')[-(1:2)],
  brewer.pal(8, 'Greys')[-(1:3)],
  brewer.pal(8, 'YlOrRd')[4:5],
  brewer.pal(8, 'Greens')[-(1:3)]
)

# Initialize color palette for cell lines
col.pal.c <- lapply(1:nrow(xg), function(i) {
  ncell <- length(xg$CellLine[[i]])
  cols <- colorRampPalette(col.pals[[i]])(ncell)
  
  cells.i <- xg$CellLine[[i]]
  age.i <- age.key[cells.i] %>% sort %>% rev
  
  names(cols) <- names(age.i)
  return(cols)
})

col.pal.c <- unlist(col.pal.c)

# Initialize color palette for genetics
col.pal.g <- sapply(1:nrow(xg), function(i) {
  ncell <- length(xg$CellLine[[i]]) + 1
  cols <- colorRampPalette(col.pals[[i]])(ncell)
  return(cols[floor(ncell / 2)])
})

names(col.pal.g) <- xg$Genetics
col.pal <- c(col.pal.c, col.pal.g)
