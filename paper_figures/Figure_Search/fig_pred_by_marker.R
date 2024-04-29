################################################################################
#+ setup, message=FALSE
################################################################################
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
source(file.path(fig.base.dir, "paper_figures", "theme.R"))
theme_set(theme_als())

# Parameters for visualization
save.fig <- TRUE

# Load in fitted model for current screen
output.fig.dir <- file.path(figure.dir, 'fig')
load(file.path(output.fig.dir, str_c('model_', screen, '.Rdata')))


################################################################################
#+ cell_prediction_plots, fig.width=18, fig.height=12, warning=FALSE
################################################################################
# Predictions by genetics
predictions <- rbindlist(fits) %>%
  mutate(Genetics=factor(Genetics, levels=genetic.order)) %>%
  filter(!str_detect(Markers, ':')) %>%
  group_by(PlateID, WellID, CellLine, Genetics, Markers) %>%
  summarize(Ypred = mean(Ypred))

nrow <- ifelse(screen == "110420_BiomarkerScreen", 2, 1)
p <- ggplot(predictions, aes(x = reorder(CellLine, Ypred), y = Ypred)) +
  geom_boxplot(aes(fill = Genetics, col = Genetics), alpha = 0.5) + 
  geom_jitter(width = 0.2, alpha = 0.8, aes(col = Genetics)) +
  facet_wrap(~Markers, nrow = nrow) + 
  scale_fill_manual(values = col.pal) +
  scale_color_manual(values = col.pal) +
  theme(legend.position = "none") +
  ylab("i-MAP score") + 
  xlab(NULL) +
  ylim(c(0, 1)) +
  theme(text = element_text(size = 13))

fout <- file.path(output.fig.dir, str_c("fig_", screen, ".pdf"))
if (save.fig) pdf(file = fout, height = nrow * 5, width=16)
plot(p)
if (save.fig) dev.off()
