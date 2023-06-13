#' # Figure S2
#' The following script evaluates primary sources of variability in raw imaging
#' feature data through PCA analysis.
#' 
#' Requires input datasets
#'   `data/032422_WTFUSSporadics/Cell_Filtered/profiles_raw.Rdata`
#' 
#' Karl Kumbier 6/9/2023
#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(hrbrthemes)
library(ggsci)
library(gridExtra)
library(caret)
library(viridis)

theme_set(
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=24, 
    axis_text_size=22,
    base_size=22,
    base_family='sans',
    axis=TRUE,
    axis_col='#000000',
    grid=FALSE
  )
)

size <- 16

# Paths to data directories
als.dir <- Sys.getenv('ALS_PAPER')
figure.dir <- file.path(als.dir, 'paper_figures', 'Figure_Scores', 'fig_qc')
dir.create(figure.dir, showWarnings=FALSE)

# Source utility functions
source(file.path(als.dir, 'paper_figures', 'utilities.R'))
source(file.path(als.dir, 'paper_figures', 'color_palette.R'))
library.file <- file.path(als.dir, 'data', 'Fibroblasts_Library.csv')

# Parameters for analysis
save.fig <- TRUE # save figures or generate in notebook

# Initialize cell line metadata table
xcell <- fread(library.file) %>% 
  dplyr::rename(CellLine=`Cell line`) %>%
  mutate(Sex=toupper(Sex)) %>%
  dplyr::rename(Genetics=`Disease Gene`) %>%
  dplyr::rename(Figure=`Figure Names`) %>%
  mutate(Genetics=str_replace_all(Genetics, 'ALS with no known mutation', 'sporadic'))

data.dir <- file.path(als.dir, 'data_profiles', '032422_WTFUSSporadics', 'Cell_Filtered')
load(file.path(data.dir, 'profiles_raw.Rdata'))

x <- left_join(x, xcell, by='CellLine') %>%
  mutate(Plate=str_remove_all(PlateID, '^.*plate_')) %>%
  mutate(DiseaseStatus=ifelse(Genetics.x != 'Healthy', 'ALS', 'Healthy')) %>%
  dplyr::rename(Age=`Age of biopsy`) %>%
  dplyr::rename(Passage=CellLinePassageNumber) %>%
  dplyr::rename(Genetics=Genetics.x) %>%
  dplyr::rename(Site=`Origin Institution`) %>%
  mutate(Passage=ifelse(CellLine == 'F7', 20, Passage))


################################################################################
# PCA plots
################################################################################
xfeat <- dplyr::select(x, matches('^X'))
xpca <- prcomp(apply(xfeat, MAR=2, rank))
pct.var <- xpca$sdev ^ 2 / sum(xpca$sdev ^ 2)

plot_pc <- function(x, var, pal, pct.var, facet=FALSE) {
  # Wrapper function to plot PCs 1-4
  
  cols <- c('PlateID', 'WellID', var)
  if (facet) cols <- c(cols, 'ALS')
  pct.var <- round(pct.var, 3) * 100
  
  p1 <- ggplot(x, aes(x=PC1, y=PC2)) +
    geom_point(aes_string(col=var)) +
    scale_y_continuous(labels = scales::scientific, breaks=c(-2e5, 0, 2e5)) +
    scale_x_continuous(labels = scales::scientific, breaks=c(-4e5, 0, 4e5)) +
    xlab(str_c('PC1, (', pct.var[1], '%)')) +
    ylab(str_c('PC2, (', pct.var[2], '%)'))
  
  p2 <- ggplot(x, aes(x=PC3, y=PC4)) +
    geom_point(aes_string(col=var)) +
    scale_y_continuous(labels = scales::scientific, breaks=c(-2e5, 0, 2e5)) +
    scale_x_continuous(labels = scales::scientific, breaks=c(-5e5, 0, 5e5)) +
    xlab(str_c('PC3, (', pct.var[3], '%)')) +
    ylab(str_c('PC4, (', pct.var[4], '%)')) 

  if (is.numeric(x[[var]])) {
    p1 <- p1 + scale_color_gradientn(colors=pal, limits=force)
    p2 <- p2 + scale_color_gradientn(colors=pal, limits=force)
  } else {
    p1 <- p1 + scale_color_manual(values=pal, limits=force)
    p2 <- p2 + scale_color_manual(values=pal, limits=force)
  }
  
  if (facet) {
    p1 <- p1 + facet_wrap(~DiseaseStatus)
    p2 <- p2 + facet_wrap(~DiseaseStatus)
  }
  
  if (!facet) return(grid.arrange(p1, p2, nrow=1))
  return(grid.arrange(p1, p2, nrow=2))
}

# Compute well-level averages of PCA-projections
cols <- c('PlateID', 'WellID', 'DiseaseStatus', 'Genetics', 'Sex', 'Age', 'Passage', 'Plate')

xplot <- cbind(xpca$x, x) %>%
  group_by_at(vars(one_of(cols))) %>%
  summarize_if(is.numeric, mean)

pdf(file.path(figure.dir, 'pca_als.pdf'), height=size / 2, width=size)
plot_pc(xplot, 'DiseaseStatus', pal_jama()(7), pct.var)
dev.off()

pdf(file.path(figure.dir, 'pca_gene.pdf'), height=size / 2, width=size)
plot_pc(xplot, 'Genetics', col.pal, pct.var)
dev.off()

pdf(file.path(figure.dir, 'pca_sex.pdf'), height=size, width=size)
plot_pc(xplot, 'Sex', pal_jama()(7), pct.var, TRUE)
dev.off()

pdf(file.path(figure.dir, 'pca_age.pdf'), height=size, width=size)
plot_pc(xplot, 'Age', viridis(10), pct.var, TRUE)
dev.off()

pdf(file.path(figure.dir, 'pca_passage.pdf'), height=size, width=size)
plot_pc(xplot, 'Passage', viridis(10), pct.var, TRUE)
dev.off()

pdf(file.path(figure.dir, 'pca_plate.pdf'), height=size, width=size)
plot_pc(xplot, 'Plate', pal_jama()(7), pct.var, TRUE)
dev.off()
