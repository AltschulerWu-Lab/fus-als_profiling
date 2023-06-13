#' # Figure S1
#' The following script fits models by marker set for the FUS, CD63, EEA1 
#' comparison screen (Fig. S1). We consider model performance with respect to
#' both marker pairs and individual markers. 
#' 
#' Requires input datasets
#'   `data_profiles/030321_FUSWT_Cell`
#'   
#' Requires models trained using scripts:
#'    `train_fus_screen.R`
#'   
#'   Karl Kumbier 2/23/2022`
#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(tidytext)
library(PRROC)
library(scales)
library(ggsci)
library(hrbrthemes)
library(gridExtra)

theme_set(
  theme_ipsum(
    plot_title_size=24,
    axis_title_size=24,
    strip_text_size=24, 
    axis_text_size=22,
    base_size=24,
    base_family='sans',
    axis=TRUE,
    axis_col='#000000',
    grid=FALSE
  )
)

# Set paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Search')

source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Set paths to data directories
data.dir <- file.path(data.base.dir, '030321_FUSWT', 'Cell_Filtered')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')

# Parameters for analysis
save.fig <- TRUE

# Fit model to select data
output.fig.dir <- file.path(figure.dir, 'fig/')
load(file.path(output.fig.dir, 'model_030321_FUSWT.Rdata'))


#' ## Accuracy by marker set
#' To evaluate performance of marker sets, we compute AUROC on held-out cell 
#' lines at three data resolutions: (i) cell line (ii) well replicate 
#' (iii) single cell. Since NS003 cannot be accurately classified using any
#' marker set, but is an outlier wrt some, we drop this line when computing
#' AUROC.
#+ accuracy
accuracy <- lapply(fits, function(z) {
  
  # Drop cell line not classified by any marker set
  z <- filter(z, CellLine != 'F4')
  
  # Compute aggregate predicitons by well, cell line
  xwell <- group_by(z, PlateID, WellID, CellLine) %>%
    summarize(Ypred=mean(Ypred), Ytest=mean(Ytest), .groups='drop')
  
  xline <- group_by(z, CellLine) %>%
    summarize(Ypred=mean(Ypred), Ytest=mean(Ytest), .groups='drop')
  
  # Compute AUROC
  rs <- roc.curve(scores.class0=z$Ypred, weights.class0=z$Ytest)
  rw <- roc.curve(scores.class0=xwell$Ypred, weights.class0=xwell$Ytest)
  rl <- roc.curve(scores.class0=xline$Ypred, weights.class0=xline$Ytest)
  
  # aggregate for output
  out <- data.frame(Marker=unique(z$Marker)) %>%
    mutate(AUROC_SC=rs$auc) %>%
    mutate(AUROC_WL=rw$auc) %>%
    mutate(AUROC_CL=rl$auc)
  
  return(out)
})

accuracy <- rbindlist(accuracy)

#' The figure below compares model prediction scores across well replicates
#' for each pari of EEA1, CD63, and FUS.
#+ prediction_scores
xcell <- fread(library.file) %>% 
  dplyr::rename(CellLine=`Cell line`) %>%
  dplyr::rename(Figure=`Figure Names`)

predictions <- rbindlist(fits) %>%
  group_by(CellLine, WellID, Markers, Ytest, PlateID, Genetics) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine') %>%
  filter(Genetics %in% c('Healthy', 'FUS-ALS'))


xplot <- lapply(unique(predictions$Markers), function(m) {
  
  predictions.m <- filter(predictions, Markers == m)
  
  m <- str_remove_all(m, '^.*:')
  if (str_detect(m, '_')) {
    m <- str_c(str_replace_all(m, '_', ', '), ' combined features')
  } else {
    m <- str_c(m, ' features only')
  }
  
  p <- predictions.m  %>%
    ggplot(aes(x=reorder(Figure, Ypred, median), y=Ypred, col=Genetics)) +
    geom_boxplot(aes(fill=Genetics), outlier.shape=NA, coef=NULL, alpha=0.75) +
    geom_jitter(width=0.15, alpha=0.9) +
    theme(legend.position='none') +
    theme(axis.text.x=element_text(angle=90)) +
    xlab(NULL) +
    ylab('FUS-ALS phenotype score') +
    scale_color_manual(values=col.pal) +
    scale_fill_manual(values=col.pal) +
    ylim(0:1) +
    ggtitle(m) +
    theme(plot.subtitle=element_text(size=14))
})


if (save.fig) pdf(file=file.path(output.fig.dir, 'fig_supp_tuning_grid.pdf'), h=22, w=22)
do.call(grid.arrange, xplot)
if (save.fig) dev.off()

################################################################################
# Summarize accuracy - population-level
################################################################################
accuracy <- mutate(accuracy, MarkerSet=str_remove_all(Marker, ':.*$')) %>%
  mutate(Marker=str_remove_all(Marker, '^.*:')) %>%
  mutate(Marker=str_replace_all(Marker, '_', ', ')) %>%
  mutate(MarkerSet=str_replace_all(MarkerSet, '_', ', '))

accuracy.g <- group_by(accuracy, Marker) %>% 
  summarize_if(is.numeric, mean) %>%
  arrange(desc(AUROC_SC))

marker.levels <- accuracy.g$Marker

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig_supp_tuning_auc.pdf'), h=8, w=22)
reshape2::melt(accuracy.g) %>%
  mutate(variable=str_replace_all(variable, 'SC', 'single cell')) %>%
  mutate(variable=str_replace_all(variable, 'WL', 'well replicate')) %>%
  mutate(variable=str_replace_all(variable, 'CL', 'cell line')) %>%
  mutate(variable=str_remove_all(variable, 'AUROC_')) %>%
  mutate(Marker=factor(Marker, levels=marker.levels)) %>%
  mutate(variable=factor(variable, levels=c('single cell', 'well replicate', 'cell line'))) %>%
  ggplot(aes(x=variable, y=value, fill=Marker)) +
  geom_bar(stat='identity', position='dodge') +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_jama() +
  xlab(NULL) +
  ylab('AUROC') +    
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(fill=NULL)
if (save.fig) dev.off()

################################################################################
# Summarize accuracy - cell line
################################################################################
if (save.fig) pdf(file=file.path(output.fig.dir, 'fig_supp_tuning_cell.pdf'), h=8, w=22)
mutate(predictions, MarkerSet=str_remove_all(Markers, ':.*$')) %>%
  mutate(Marker=str_remove_all(Markers, '^.*:')) %>%
  mutate(Marker=str_replace_all(Markers, '_', ', ')) %>%
  mutate(Marker=str_remove_all(Marker, '^.*:')) %>%
  mutate(Marker=factor(Marker, levels=marker.levels)) %>%
  mutate(MarkerSet=str_replace_all(MarkerSet, '_', ', ')) %>%
  mutate(Count=str_detect(Marker, ',')) %>%
  group_by(CellLine, WellID, PlateID, Marker, Count, Figure) %>%
  summarize(Ypred=mean(Ypred)) %>%
  ggplot(aes(x=Count, y=Ypred, fill=Marker, color=Marker)) +
  geom_boxplot(alpha=0.8) +
  facet_wrap(~Figure, nrow=3) +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.text.x=element_blank()) +
  theme(legend.position='none') +
  ylab(NULL) +
  xlab(NULL) +
  ylim(0:1)
if (save.fig) dev.off()
