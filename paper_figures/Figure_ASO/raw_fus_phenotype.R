#' # Figure 5
#' The following script generates plots for qPCR levels following ASO knockdown
#' of FUS (Fig. S5a) and FUS intensity levels in microscopy screens (Fig. 5a).
#' 
#' Requires input datasets
#'   `data/2022_ASOscreen_v3_qPCR results.xls`
#'   `data/080122_ASOScreenV2_Cell`
#'   
#' Karl Kumbier, 3/6/2023
library(data.table)
library(readxl)
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
treatment.levels <- c('untreated', 'H2O', 'NTC 10uM', 'ASO 1uM', 'ASO 10uM')
heat.pal <- c('#FFFFFF', pal_material('blue-grey')(10))

# Set paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_ASO')

source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Set paths to data directories
data.dir <- file.path(data.base.dir, '080122_ASOScreenV2', 'Cell_Filtered')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')
qpcr.file <- file.path(fig.base.dir, 'data', '2022_ASOscreen_qPCR.xls')

output.fig.dir <- file.path(figure.dir, 'fig/')
dir.create(output.fig.dir, recursive=TRUE, showWarnings=FALSE)


# Load mutation table
xcell <- fread(library.file) %>%
  dplyr::rename(Figure=`Figure Names`) %>%
  dplyr::rename(CellLine=`Cell line`)

mutation.key <- setNames(xcell$Figure, xcell$CellLine)

################################################################################
# FUS qPCR by cell line.
################################################################################
qpcr <- read_excel(qpcr.file, sheet=5)

colnames(qpcr)[1] <- 'CellLine'
colnames(qpcr) <- str_remove_all(colnames(qpcr), ' ')

p <- reshape2::melt(qpcr) %>%
  filter(variable != 'Untreated') %>%
  mutate(variable=str_replace_all(variable, 'ASO', ' ASO')) %>%
  mutate(variable=str_replace_all(variable, 'NTC', 'NTC 10uM')) %>%
  mutate(variable=str_replace_all(variable, 'low ASO', 'ASO 1uM')) %>%
  mutate(variable=str_replace_all(variable, 'high ASO', 'ASO 10uM')) %>%
  mutate(CellLineL=mutation.key[CellLine]) %>%
  mutate(ASO=factor(variable, levels=treatment.levels)) %>%
  mutate(CellLineL=str_remove_all(CellLineL, '-.*$')) %>%
  ggplot(aes(x=CellLine, y=value, fill=ASO)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  scale_fill_manual(values=heat.pal[c(3, 6, 9)]) +
  xlab(NULL) +
  ylab('qPCR-FUS') +
  facet_wrap(~CellLineL, nrow=2, scales='free_x') +
  theme(axis.text.x=element_blank())

pdf(file.path(output.fig.dir, 'fig_supp_qpcr.pdf'), height=8, width=16)
plot(p)
dev.off()


################################################################################
# Raw nuclear and cytoplasmic FUS intensities
################################################################################
load(file.path(data.dir, 'profiles_raw.Rdata'))
x <- mutate(x, ASO=str_remove_all(ASOtreatment, 'EP_')) %>%
  mutate(ASO=str_replace_all(ASO, 'ASO', ' ASO')) %>%
  mutate(ASO=str_replace_all(ASO, 'NTC', 'NTC 10uM')) %>%
  mutate(ASO=str_replace_all(ASO, 'low ASO', 'ASO 1uM')) %>%
  mutate(ASO=str_replace_all(ASO, 'high ASO', 'ASO 10uM')) %>%
  mutate(ASO=factor(ASO, treatment.levels)) %>%
  mutate(Genetics=str_replace_all(Genetics, 'FUS', 'FUS-ALS'))

# Plot change in nuclear/cytoplasmic FUS
fselect <- c('X..Iav_in_dna_region..2', 'X..Iav_in_non_dna_region..2')
mselect <- c('Genetics', 'ASO')

# Compute well-level median intensities
xgroup <- group_by(x, ID, ASO, Genetics) %>% 
  filter(ASO != 'H2O') %>%
  filter(ASO != 'untreated') %>%
  summarize_if(is.numeric, median) %>%
  ungroup() %>%
  dplyr::select(-ID)

comparisons <- list(c("NTC 10uM", "ASO 1uM"), c("ASO 1uM", "ASO 10uM"))

p1 <- dplyr::rename(xgroup, NucFUS=`X..Iav_in_dna_region..2`) %>% 
  dplyr::select(one_of(c('NucFUS', mselect))) %>%
  ggbarplot(x="ASO", y="NucFUS", fill="ASO", facet.by="Genetics", add="mean_sd") +
  ylab('Intensity') +
  scale_fill_manual(values=heat.pal[c(3, 6, 9)], labels=NULL) +
  xlab(NULL) +
  ylab('Nuclear FUS pixel intensity') +
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=20, 
    axis_text_size=20,
    base_size=20,
    base_family='sans',
    axis=TRUE,
    axis_col='#000000',
    grid=FALSE
  ) +
  theme(legend.position='none') +
  stat_compare_means(comparisons=comparisons, label="p.signif", size=7, p.adjust.method='bonferroni') +
  theme(axis.text.x=element_blank()) 

p2 <- dplyr::rename(xgroup, CytFUS=`X..Iav_in_non_dna_region..2`) %>% 
  dplyr::select(one_of(c('CytFUS', mselect))) %>%
  ggbarplot(x="ASO", y="CytFUS", fill="ASO", facet.by="Genetics", add="mean_sd")+
  ylab('Intensity') +
  scale_fill_manual(values=heat.pal[c(3, 6, 9)], labels=NULL) +
  xlab(NULL) +
  ylab('Cytosolic FUS pixel intensity') +
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=20, 
    axis_text_size=20,
    base_size=20,
    base_family='sans',
    axis=TRUE,
    axis_col='#000000',
    grid=FALSE
  ) +
  theme(legend.position='none') +
  stat_compare_means(comparisons=comparisons, label="p.signif", size=7, p.adjust.method='bonferroni') +
  theme(axis.text.x=element_blank())

p3 <- mutate(xgroup, ratio=(`X..Iav_in_dna_region..2` / `X..Iav_in_non_dna_region..2`)) %>% 
  dplyr::select(one_of(c('ratio', mselect))) %>%
  ggbarplot(x="ASO", y="ratio", fill="ASO", facet.by="Genetics",add="mean_sd")+
  ylab('Intensity') +
  scale_fill_manual(values=heat.pal[c(3, 6, 9)]) +
  labs(fill=NULL) +
  xlab(NULL) +
  ylab('Nuclea') +
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=20, 
    axis_text_size=20,
    base_size=20,
    base_family='sans',
    axis=TRUE,
    axis_col='#000000',
    grid=FALSE
  ) +
  stat_compare_means(comparisons=comparisons, label="p.signif", size=7, p.adjust.method='bonferroni') +
  ylab('Nuclear / cytosolic FUS pixel intensity ratio') +
  theme(axis.text.x=element_blank())



if (save.fig) pdf(file=str_c(output.fig.dir, 'fig5a_nuc.pdf'), height=8, width=10)
plot(p1)
if (save.fig) dev.off()

if (save.fig) pdf(file=str_c(output.fig.dir, 'fig5a_cyt.pdf'), height=8, width=10)
plot(p2)
if (save.fig) dev.off()

if (save.fig) pdf(file=str_c(output.fig.dir, 'fig5a_rat.pdf'), height=8, width=10)
plot(p3)
if (save.fig) dev.off()
