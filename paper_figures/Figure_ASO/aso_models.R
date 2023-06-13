#' # Figure 5c-d
#' The following script fits and evaluates FUS vs. WT classifiers on ASO-treated
#' cell lines. Models can be fit for the full feature set or cytosolic features
#' only by setting the variable `cytosolic = FALSE` or `cytosolic = TRUE`
#' respectively. Models can be fit using EEA1/FUS features or only FUS features 
#' by setting the variable `fus.only = FALSE` or `fus.only = TRUE` respectively.
#' 
#' Models are trained using data from EP cells and evaluated on EP cells treated
#' with NTC, low ASO and high ASO.
library(data.table)
library(readxl)
library(tidyverse)
library(tidytext)
library(ggsci)
library(ggpubr)
library(hrbrthemes)
library(caret)
library(gridExtra)
library(see)

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

# Paths to data directories
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.base.dir <- Sys.getenv('ALS_DATA')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_ASO')
platemap.dir <- file.path(fig.base.dir, 'data', 'platemaps')
dir.create(platemap.dir)

source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Set paths to data directories
data.dir <- file.path(data.base.dir, '080122_ASOScreenV2', 'Cell_Filtered')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')
plate.file <- file.path(fig.base.dir, 'data', 'plate_list.txt')

# Set parameters for analysis
n.core <- 16
compound <- 'DMSO'
marker <- 'FUS_EEA1'
cytosolic <- FALSE
fus.only <- FALSE
save.fig <- TRUE
save.plates <- TRUE
train.condition <- 'H2O'
treatment.levels <- c('untreated', 'H2O', 'NTC 10uM', 'ASO 1uM', 'ASO 10uM')

# Set paths to data directories
output.fig.dir <- file.path(figure.dir, 'fig')
if (cytosolic) output.fig.dir <- str_c(output.fig.dir, '_cytosolic')
if (fus.only) output.fig.dir <- str_c(output.fig.dir, '_fus')
dir.create(output.fig.dir, recursive=TRUE, showWarnings=FALSE)

# Load mutation table
xcell <- fread(library.file) %>% 
    dplyr::rename(CellLine=`Cell line`) %>%
    mutate(Sex=toupper(Sex)) %>%
    dplyr::rename(Genetics=`Disease Gene`) %>%
    mutate(Genetics=ifelse(Genetics == 'WT', 'Healthy', Genetics)) %>%
    mutate(Genetics=ifelse(Genetics == 'FUS', 'FUS-ALS', Genetics)) %>%
    mutate(Genetics=str_replace_all(Genetics, 'ALS with no known mutation', 'sporadic'))

# Initialize models to be used for supervised learning
model <- function(x, y) irf(x, y, n.core=n.core)
model_predict <- predict_irf
importance_norm <- function(x) as.matrix(x / max(x))

#' ## Loading data
#' Our dataset consists of phenotypic profiles derived from `r marker`. Each 
#' observation in the dataset corresponds to a single `r type`. We subset to
#' features as indicated by `cytosolic` and `fus.only` variables.
#+load_data
# Load in dataset for select marker set
load(file.path(data.dir, str_c('profiles_', marker, '.Rdata')))

x <- mutate(x, PlateID=str_remove_all(PlateID, '.*plate_')) %>%
  mutate(ID=str_c(ID, 1:n())) %>%
  dplyr::rename(ASO=ASOtreatment) %>%
  mutate(ASO=str_remove_all(ASO, 'EP_')) %>%
  mutate(ASO=str_replace_all(ASO, 'ASO', ' ASO')) %>%
  mutate(ASO=str_replace_all(ASO, 'NTC', 'NTC 10uM')) %>%
  mutate(ASO=str_replace_all(ASO, 'low ASO', 'ASO 1uM')) %>%
  mutate(ASO=str_replace_all(ASO, 'high ASO', 'ASO 10uM')) %>%
  mutate(ASO=factor(ASO, treatment.levels)) %>%
  dplyr::select(-Genetics) %>%
  left_join(xcell, by='CellLine')

# Initialize gene/cell line key
genetic.key <- dplyr::select(x, CellLine, Genetics) %>% distinct()
mutation.key <- setNames(xcell$`Figure Names`, xcell$CellLine)

# Initialize channel / region IDs
channels <- sapply(str_split(colnames(x), '\\.\\.'), '[', 3)

regions <- sapply(str_split(colnames(x), '\\.\\.'), '[', 2)
regions <- lapply(regions, str_split, pattern='_') 
regions <- sapply(regions, function(z) tail(z[[1]], 1))

# Initialize and (optionally) drop FUS nuclear features
channels <- str_split(channels, '')
regions <- str_split(regions, '')

nuc.fus <- mapply(function(ch, rg) {
  any(ch == 2 & rg == 'D') | any(ch == 2 & rg == 'C')
}, channels, regions)

nuc.fus[is.na(nuc.fus)] <- FALSE
if (cytosolic) x <- dplyr::select_if(x, !nuc.fus)

# Initialize and (optionally) drop EEA1 features
channels <- sapply(str_split(colnames(x), '\\.\\.'), '[', 3)
channels[is.na(channels)] <- 0
if (fus.only) x <- dplyr::select_if(x, !str_detect(channels, '3'))


#' ## ASO models
#' We evaluate the degree to which ASO treatment "pushes" cell lines (disease 
#' to health, health to disease) in the context of supervised models trained
#' to discriminate between ALS/FUS. Models are trained on DMSO treated cell
#' lines and evaluated on lines treated with (i) H2O (ii) NTC (iii) low ASO
#' and (iv) high ASO.
#+ leave_cell_out
aso_fit <- function(x,
                    model=irf,
                    model_predict=predict_irf) {
  
  # Check for valid input
  if (nrow(x) == 0) return(NULL)
  
  # Set feature matrix and reponse vector
  y <- as.numeric(x$Genetics != 'Healthy')
  
  # Downsample for class balance
  xds <- downSample(x, as.factor(y))
  xx <- dplyr::select(xds, matches('^X'))  
  y <- as.numeric(xds$Class) - 1
  
  # Set training/test indices
  id.train <- which(xds$ASO == train.condition)
  id.test <- which(xds$ASO != train.condition)
  fit <- fit_model(xx, y, id.train, model, model_predict)
  
  predicted <- data.frame(
    Ypred=fit$ypred,
    Ytest=y[id.test],
    Genetics=xds$Genetics[id.test],
    CellLine=xds$CellLine[id.test],
    Treatment=xds$ASO[id.test],
    WellID=xds$WellID[id.test],
    ID=xds$ID[id.test],
    Marker=marker,
    Plate=xds$PlateID[id.test]
  )
  
  return(list(fit=fit, predicted=predicted))
}


# Filter data to selected screen/compound
set.seed(47)
cell.models <- aso_fit(x, model, model_predict)

#' ### Fig 3c
#' Model prediction scores by genetics, treatment. 
#+ raw_predictions, fig.height=8, fig.width=12
if (cytosolic & fus.only) {
  main <- 'FUS only,\nexcluding nuclear FUS features'
} else if (cytosolic) {
  main <- 'FUS & EEA1,\nexcluding nuclear FUS features'
} else if (fus.only) {
  main <- 'FUS only,\nall features'
} else {
  main <- 'FUS & EEA1,\nall features'
}

# Comparisons of raw predictions
xplot <- filter(cell.models$predicted, !Treatment %in% train.condition) %>%
  group_by(Genetics, CellLine, Treatment) %>%
  summarize(Ypred=mean(Ypred), .groups='drop') %>%
  mutate(ChannelS='FUS, EEA1') %>% 
  filter(Treatment != 'untreated') %>%
  mutate(Treatment=factor(Treatment, levels=levels(Treatment)[3:5]))

comparison <- list(c("Healthy", "FUS-ALS"))
                   
p <- ggboxplot(xplot, x="Genetics", fill="Genetics", y="Ypred", facet.by='Treatment', 
          outlier.shape=NA, coef=NULL) +
  stat_boxplot(coef=NULL, aes(fill=Genetics)) +
  geom_jitter(size=3, width=0.1) +
  stat_compare_means(comparisons=comparison, label="p.signif", size=7, p.adjust.method='bonferroni') +
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_manual(values=col.pal) +
  ylim(0:1) +
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
  theme(axis.text.x=element_blank())

pdf(file.path(output.fig.dir, 'fig5c.pdf'), height=8, width=8)
plot(p)
dev.off()

#' ### Fig 3d
#' Model prediction scores by cell line, treatment. 
#+ raw_predictions_cell_line, fig.height=8, fig.width=12
# Raw predictions
p <- filter(cell.models$predicted, !Treatment %in% train.condition) %>%
  mutate(Treatment=factor(Treatment, levels=treatment.levels)) %>%
  group_by(Genetics, CellLine, Treatment, WellID, Plate) %>%
  summarize(Ypred=mean(Ypred), Count=n(), .groups='drop') %>%
  filter(Treatment != 'untreated') %>%
  mutate(CellLineL=mutation.key[CellLine]) %>%
  ggplot(aes(x=Treatment, y=Ypred, fill=Genetics)) +
  geom_boxplot(outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.25) +
  facet_wrap(~CellLineL, nrow=2) +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=col.pal) +
  scale_color_manual(values=col.pal) +
  ggtitle(main) +

pdf(file.path(output.fig.dir, 'fig5d.pdf'), height=9, width=20)
plot(p)
dev.off()

save(file=file.path(output.fig.dir, 'model.Rdata'), cell.models)

# Save table of plates used
if (save.plates) {
  write.table(
    file=plate.file, 
    unique(x$PlateID), 
    append=TRUE, 
    row.names=FALSE, 
    col.names=FALSE, 
    quote=FALSE
  )
  
  # Save plate maps
  plate.layout <- cell.models$predicted %>%
    mutate(PlateID=str_remove_all(ID, '_[0-9]*-[0-9]*$')) %>%
    dplyr::select(PlateID, CellLine, WellID) %>%
    distinct() %>%
    mutate(Row=as.numeric(str_remove_all(WellID, '-.*$'))) %>%
    mutate(Col=as.numeric(str_remove_all(WellID, '^.*-')))
  
  for (p in unique(plate.layout$PlateID)) {
    platemap <- matrix('', nrow=16, ncol=24)
    pm <- filter(plate.layout, PlateID == p)
    for (i in 1:nrow(pm)) platemap[pm$Row[i], pm$Col[i]] <- pm$CellLine[i]
    colnames(platemap) <- 1:24
    
    print(table(platemap))
    fout <- file.path(platemap.dir, str_c(p, '.csv'))
    write.csv(file=fout, platemap)
  }
}
