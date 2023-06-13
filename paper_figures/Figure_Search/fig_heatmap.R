#' # Figure 1a
#' The following script generates heatmaps of prediction accuracy by genetics 
#' and marker set. Predictions are averaged across all cells from a given line 
#' or well replicate. For each subtype, we report (i) recall for ALS cell lines
#' and (ii) AUROC evaluated on single cell predictions. Class predictions 
#' at the cell line level are defined by thresholding prediction scores at the 
#' maximum Healthy prediction score. Under this definition, all Healthy cell 
#' lines are accurately classified.
#' 
#' Requires input datasets
#'   `data_profiles/110420_BiomarkerScreen/Cell_Filtered` 
#'   `data_profiles/121820_SporadicScreen/Cell_Filtered`
#'   `data_profiles/030321_FUSWT/Cell_Filtered`
#' 
#' Requires models trained using scripts:
#'    `train_marker_screens.R`
#'    `train_fus_screen.R`
#' 
#' Karl Kumbier, 2/27/2023
#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(scales)
library(superheat)
library(ggsci)

# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA') 
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Search')

#screen <- '110420_BiomarkerScreen'
#screen <- '121820_SporadicScreen'
screen <- '030321_FUSWT'

# Source utility functions
source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Parameters for analysis
save.fig <- TRUE # save figures (TRUE) or generate in notebook (FALSE)
heat.pal <- RColorBrewer::brewer.pal(8, 'Blues')
heat.thresh <- 0.8

# Initialize plot ordering for genetics
genetic.order <- c('Healthy', 'TBK1', 'sporadic', 'C9orf72', 'FUS-ALS')

# Load in fitted model for current screen
output.fig.dir <- file.path(figure.dir, 'fig')
load(file.path(output.fig.dir, str_c('model_', screen, '.Rdata')))


################################################################################
# Generate plots
################################################################################
#+ cell_prediction_plots, fig.width=18, fig.height=12, warning=FALSE
# Predictions by genetics
predictions <- rbindlist(fits) %>%
  mutate(Genetics=factor(Genetics, levels=genetic.order)) %>%
  filter(!str_detect(Markers, ':'))

markers <- unique(predictions$Markers)
genetics <- unique(predictions$Genetics) %>% as.character
genetics <- intersect(genetic.order, genetics)

# Compute recall for cell line aggregate predictions
cell.line.acc <- sapply(genetics, function(g) {
  sapply(markers, function(m) {
    
    if (g != 'Healthy') {
      pm <- filter(predictions, Markers == m, Genetics %in% c(g, 'Healthy')) %>%
        group_by(CellLine, Ytest) %>%
        summarize(Ypred=mean(Ypred), .groups='drop') %>%
        mutate(YpredT=as.numeric(Ypred > max(Ypred[Ytest == 0])))
    
      nals <- sum((pm$Ytest == pm$YpredT)[pm$Ytest == 1])
      ntot <- sum(pm$Ytest == 1)
      return(str_c(nals, ' / ', ntot))
    } else {
      pm <- filter(predictions, Markers == m) %>%
        group_by(CellLine, Ytest) %>%
        summarize(Ypred=mean(Ypred), .groups='drop') %>%
        mutate(YpredT=as.numeric(Ypred > max(Ypred[Ytest == 0])))
      
      nwt <- sum((pm$Ytest == pm$YpredT)[pm$Ytest == 0])
      ntot <- sum(pm$Ytest == 0)
      return(str_c(nwt, ' / ', ntot))
    }
  })
})

# Compute AUROC for single cell predictions
single.cell.acc <- sapply(genetics, function(g) {
  sapply(markers, function(m) {
    if (g != 'Healthy') {
      
      pm <- filter(predictions, Markers == m, Genetics %in% c(g, 'Healthy'))
      roc <- roc.curve(scores.class0=pm$Ypred, weights.class0=pm$Ytest)
      return(pmax(roc$auc, 0.5))
    } else {
      pm <- filter(predictions, Markers == m)
      roc <- roc.curve(scores.class0=pm$Ypred, weights.class0=pm$Ytest)
      return(pmax(roc$auc, 0.5))
    }
  })
})

# Threshold for visualization
single.cell.acc[single.cell.acc > heat.thresh] <- heat.thresh

# Initialize row/column ordering for visualization
row.order <- order(single.cell.acc[,'FUS-ALS'])

# Initialize color palette for visualization
bottom.label.col <- col.pal[colnames(single.cell.acc)]

# Format table names for visualization
rownames(single.cell.acc) <- rownames(single.cell.acc) %>%
  str_replace_all('_', ', ') %>%
  str_remove_all('\\.t$')

gene.table <- dplyr::select(predictions, Genetics, CellLine) %>% 
  distinct() %>%
  group_by(Genetics) %>%
  count() %>%
  arrange(match(Genetics, colnames(single.cell.acc)))

colnames(single.cell.acc) <- str_c(colnames(single.cell.acc), ' (', gene.table$n, ')')
colnames(cell.line.acc) <- colnames(single.cell.acc)

# Format figure size
height <- 8
if (screen == '110420_BiomarkerScreen') height <- 11

width <- 10
if (screen == '030321_FUSWT') width <- 8

bls <- 0.7
if (screen == '110420_BiomarkerScreen') bls <-  0.3

fout <- file.path(output.fig.dir, str_c('fig1a_', screen, '.png'))
if (save.fig) png(fout, height=height, width=width, units='in', res=300)
superheat(
  X=single.cell.acc[row.order,],
  X.text=cell.line.acc[row.order,],
  X.text.size=8,
  left.label='variable',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  left.label.text.size=9.5,
  left.label.size=1.5,
  grid.hline.size=1,
  grid.vline.size=1,
  bottom.label.size=bls,
  bottom.label.text.size=9.5,
  bottom.label.col=bottom.label.col,
  bottom.label.text.angle=90,
  heat.lim=c(0.5, heat.thresh)
)
if (save.fig) dev.off()
