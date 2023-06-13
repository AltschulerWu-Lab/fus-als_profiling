#' # Figure 5b
#' The following script generates plots of feature importance for cytosolic and
#' full models (Fig. 5b). 
#' 
#' Requires input datasets
#'   `data/080122_ASOScreenV2_Cell/`
#'   
#' **Note:** full and cytosolic models must be fit prior to running this script.
#' Models can be fit with calls to `fig3_cd.R`.
#' 
#' Karl Kumbier 10/19/2022
library(data.table)
library(tidyverse)
library(tidytext)
library(ggsci)
library(ggpubr)
library(hrbrthemes)

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

# Set parameters for analysis
save.fig <- TRUE

# Set paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure5')

source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Set paths to data directories
data.dir <- file.path(data.base.dir, '080122_ASOScreenV2', 'Cell_Filtered')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')

output.fig.dir <- file.path(figure.dir, 'fig5')
dir.create(output.fig.dir, recursive=TRUE, showWarnings=FALSE)

# Load in fitted models
load(file.path(output.fig.dir, 'model.Rdata'))
cell.models.full <- cell.models

load(file.path(str_c(output.fig.dir, '_cytosolic'), 'model.Rdata'))
cell.models.cyto <- cell.models
  
#' ### Fig 3E
#' Model feature importance
#+ feature_importance
importance.full <- data.frame(Importance=cell.models.full$fit$fit$importance) %>%
  mutate(Feature=names(cell.models.full$fit$fit$importance)) %>%
  mutate(Type='all features')

importance.cyto <- data.frame(Importance=cell.models.cyto$fit$fit$importance) %>%
  mutate(Feature=names(cell.models.cyto$fit$fit$importance)) %>%
  mutate(Type='excluding nuclear FUS')

importance <- rbind(importance.full, importance.cyto) %>%
  mutate(Feature=str_remove_all(Feature, '^X\\.\\.')) %>%
  mutate(Feature=str_remove_all(Feature, '\\.\\.[0-9]*$')) %>%
  group_by(Feature, Type) %>%
  summarize(Importance=max(Importance), .groups='drop') %>%
  group_by(Type) %>%
  mutate(Importance=Importance / max(Importance)) %>%
  ungroup()

top.features <- group_by(importance, Type) %>% top_n(5, Importance)


clean_feature <- function(x) {
  # Reformatting for feature name
  if (!str_detect(x, 'Intensity Ratio')) return(x)
  
  xsplit <- str_split(x, '\n ')[[1]]
  feature <- xsplit[1]
  xloc.split <- str_split(xsplit[2], ', ')[[1]] %>% str_split('/')
  
  region <- str_c(
    str_c(xloc.split[[1]][1], ' ', xloc.split[[2]][1]), ' / ',
    str_c(xloc.split[[1]][2], ' ', xloc.split[[2]][2])
  )
  
  return(str_c(feature, '\n', region))
}

p <- filter(importance, Feature %in% top.features$Feature) %>%
  mutate(Feature=str_replace_all(Feature, 'Number', 'number')) %>%
  mutate(Feature=str_replace_all(Feature, '32', 'EEA1/FUS')) %>%
  mutate(Feature=str_replace_all(Feature, '23', 'FUS/EEA1')) %>%
  mutate(Feature=str_replace_all(Feature, '3', 'EEA1')) %>%
  mutate(Feature=str_replace_all(Feature, '2', 'FUS')) %>%
  mutate(Feature=str_replace_all(Feature, 'DN', 'nucleus/cytoplasm,')) %>%
  mutate(Feature=str_replace_all(Feature, 'DC', 'nucleus/cell,')) %>%
  mutate(Feature=str_replace_all(Feature, 'NC', 'cytoplasm/cell,')) %>%
  mutate(Feature=str_replace_all(Feature, 'CC', 'cell/cell,')) %>%
  mutate(Feature=str_replace_all(Feature, 'NN', 'cytoplasm/cytoplasm,')) %>%
  mutate(Feature=str_replace_all(Feature, 'N', 'cytoplasm,')) %>%
  mutate(Feature=str_replace_all(Feature, 'C', 'cell ')) %>%
  mutate(Feature=str_replace_all(Feature, 'D', 'nucleus ')) %>%
  mutate(Feature=str_replace_all(Feature, 'cell,o', 'Co')) %>%
  mutate(Feature=str_replace_all(Feature, '_', ' ')) %>%  
  mutate(Feature=str_replace_all(Feature, 'Area', 'Area\n')) %>%
  mutate(Feature=str_replace_all(Feature, 'Intensity Ratio', 'Intensity Ratio\n')) %>%
  mutate(Feature=str_replace_all(Feature, 'Intensity Region', 'Intensity Region\n')) %>%
  mutate(Feature=str_replace_all(Feature, '  ', ' ')) %>%
  mutate(Feature=str_replace_all(Feature, '\\.\\.', ' ')) %>%
  mutate(Feature=sapply(Feature, clean_feature)) %>%
  ggplot(aes(x=reorder(Feature, Importance), y=Importance)) +
  geom_bar(aes(fill=Type), position=position_dodge(preserve="single"), stat='identity', color='black') +
  theme(axis.text.x=element_text(angle=90)) +
  coord_flip() +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values=pal_jama()(7)[c(5, 7)]) +
  theme(legend.position=c(0.7, 0.1)) +
  labs(fill=NULL)

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig5b.pdf'), height=8, width=10)
plot(p)
if (save.fig) dev.off()
