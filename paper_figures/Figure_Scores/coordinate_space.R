#' # Figure 2e
#' The following script projects cell lines from FUS vs. WT age matched screen 
#' into feature space based on important features as determined by FUS vs. WT 
#' classification models (paper Fig 2e). We plot raw image features to assess 
#' whether eigen features & model feature importance recapitulate raw signal.
#' 
#' Requires input datasets
#'   `data_profiles/032422_WTFUSSporadics_Well`
#'   
#' Karl Kumbier 10/18/2022
#+ setup, message=FALSElibrary(tidyverse)
library(data.table)
library(gridExtra)
library(ggsci)
library(hrbrthemes)
library(ggpubr)
library(patchwork)
library(Cairo)

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

fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Scores')

source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Set paths to data directories
data.dir <- file.path(data.base.dir, '032422_WTFUSSporadics', 'Well')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')

# Fit model to select data
output.fig.dir <- file.path(figure.dir, 'fig_full/')
load(file.path(data.dir, 'profiles_raw.Rdata'))

# Analysis parameters
marker <- 'FUS_EEA1'
f1 <- 'X..Iav_in_non_dna_region..2' # x-axis feature
f2 <- 'X..Iav_in_cell_region..3' # y-axis feature
save.fig <- FALSE

################################################################################
# Load in raw feature and merge with cell line metadata
################################################################################
# Read in cell line genetics
xgen <- fread(library.file) %>%
  dplyr::rename(CellLine=`Cell line`, Genetics=`Disease Gene`) %>%
  mutate(Genetics=ifelse(Genetics == 'FUS', 'FUS-ALS', Genetics)) %>%
  mutate(Genetics=ifelse(Genetics == 'WT', 'Healthy', Genetics))

x <- left_join(x, xgen, by='CellLine')

# Filter to FUS/WT cell lines
x <- filter(x, Genetics %in% c('FUS-ALS', 'Healthy'))

################################################################################
# Plot raw feature distributions
################################################################################
pxlim <- range(x[[f1]])
pylim <- range(x[[f2]])

# Generate raw feature plot
g <- 'FUS-ALS'
cols <- unlist(col.pal.c[c('Healthy', g)])
names(cols) <- str_remove_all(names(cols), '^.*\\.')

# Compute centroids & standard errors
xmean <- dplyr::select(x, one_of(f1, f2, 'CellLine', 'Genetics', 'Mutation')) %>%
  group_by(CellLine, Genetics, Mutation) %>% 
  summarize_if(is.numeric, mean) %>%
  mutate(Mutation=str_remove_all(Mutation, 'FUS-')) %>%
  mutate(Mutation=str_remove_all(Mutation, 'FUS '))

xsd <- dplyr::select(x, one_of(f1, f2, 'CellLine')) %>%
  group_by(CellLine) %>% 
  summarize_if(is.numeric, sd)

xrange <- data.frame(CellLine=xmean$CellLine, Genetics=xmean$Genetics) %>%
  mutate(MinX=xmean[[f1]] - xsd[[f1]]) %>%
  mutate(MaxX=xmean[[f1]] + xsd[[f1]]) %>%
  mutate(AvgX=xmean[[f1]]) %>%
  mutate(MinY=xmean[[f2]] - xsd[[f2]]) %>%
  mutate(MaxY=xmean[[f2]] + xsd[[f2]]) %>%
  mutate(AvgY=xmean[[f2]])

pxlim <- c(min(xrange$MinX), max(xrange$MaxX))
pylim <- c(min(xrange$MinY), max(xrange$MaxY))

# Generate plots
s1 <- ggplot(xmean, aes_string(col='Genetics', x=f1, y=f2)) +
  geom_point(size=4) +
  geom_segment(data=xrange, aes(x=MinX, xend=MaxX, y=AvgY, yend=AvgY), alpha=0.5, size=1.5) +
  geom_segment(data=xrange, aes(x=AvgX, xend=AvgX, y=MinY, yend=MaxY), alpha=0.5, size=1.5) +
  xlim(pxlim) +
  ylim(pylim) +
  xlab('Cytosolic FUS intensity (KS)') +
  ylab('Cellular EEA1 intensity (KS)') +
  scale_color_manual(values=col.pal) +
  theme(legend.position='none') +
  theme(plot.margin=ggplot2::margin(0, 0, 0, 0))

d1 <- filter(xmean, Genetics %in% c('Healthy', g)) %>%
  ggplot(aes(fill=Genetics)) +
  geom_boxplot(aes_string(x=f1), coef=NULL) + 
  scale_fill_manual(values=col.pal) +
  theme_void() + 
  theme(legend.position="none") +
  theme(plot.title=element_text(face="bold", size=18)) +
  xlim(pxlim)

d2 <- filter(xmean, Genetics %in% c('Healthy', g)) %>%
  ggplot(aes(fill=Genetics)) +
  geom_boxplot(aes_string(x=f2), coef=NULL) + 
  scale_fill_manual(values=col.pal) +
  theme_void() + 
  theme(legend.position="none") +
  coord_flip() +
  xlim(pylim)

pg <- d1 + plot_spacer() + s1 + d2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 0.5), heights = c(0.5, 4))

if (save.fig) pdf(str_c(output.fig.dir, 'fig2e.pdf'), height=8, width=8)
plot(pg)
if (save.fig) dev.off()