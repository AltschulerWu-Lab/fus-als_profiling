#' # Figure 2
#' The following script fits FUS vs. Healthy classifiers and generates plots of:
#' predicted probability FUS by genetics (Fig. 2a),and cell line (Fig. 2b), ROC
#' curves for model predictions (Fig. 2c), feature importance (Fig. 2d), and
#' against patient age (Fig 2g). 
#' 
#' The script also generates plots for age of onset vs. model predictions 
#' (Fig. s3), and Healthy vs. sporadic predictions (Fig. 4a-b).
#' 
#' Requires input datasets
#'   `data_profiles/032422_HealthyFUSSporadics_Cell`
#'   
#' Karl Kumbier 10/18/2022
#+ setup, message=FALSE
library(tidyverse)
library(data.table)
library(tidytext)
library(scales)
library(parallel)
library(lme4)
library(ggsci)
library(hrbrthemes)
library(ggrepel)
library(PRROC)

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
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Scores')
platemap.dir <- file.path(fig.base.dir, 'data', 'platemaps')
dir.create(platemap.dir)

source(file.path(fig.base.dir, 'paper_figures', 'utilities.R'))
source(file.path(fig.base.dir, 'paper_figures', 'color_palette.R'))

# Set paths to data directories
data.dir <- file.path(data.base.dir, '032422_WTFUSSporadics', 'Cell_Filtered')
library.file <- file.path(fig.base.dir, 'data', 'Fibroblasts_Library.csv')
plate.file <- file.path(fig.base.dir, 'data', 'plate_list.txt')

# Parameters for analysis
n.core <- 16
type <- 'Cell'
cytosolic <- FALSE
marker <- 'FUS_EEA1'
save.fig <- TRUE
save.plate <- FALSE

# Initialize models to be used for supervised learning
model <- function(x, y) irf(x, y, n.core=n.core)
model_predict <- predict_irf
importance_norm <- function(x) as.matrix(x / max(x))

# Initialize output directories
ext <- ifelse(cytosolic, '_cytosolic', '_full')
output.fig.dir <- file.path(figure.dir, str_c('fig', ext))
dir.create(output.fig.dir, showWarnings=FALSE)

#' # Loading data
#' Our dataset consists of phenotypic profiles derived from `r marker`. Each 
#' observation in the dataset corresponds to a single `r type`.
#+load_data
# Load in datasets for selected marker
input.file <- str_c('profiles_', marker, '.Rdata')
load(file.path(data.dir, input.file))

# Initialize cell line metadata table
xcell <- fread(library.file) %>% dplyr::rename(CellLine=`Cell line`)

xcell <- dplyr::select(x, CellLine, CellLinePassageNumber) %>%
  distinct() %>%
  left_join(xcell, by='CellLine') %>%
  mutate(Sex=toupper(Sex)) %>%
  dplyr::rename(Genetics=`Disease Gene`) %>%
  mutate(Genetics=ifelse(Genetics == 'WT', 'Healthy', Genetics)) %>%
  mutate(Genetics=ifelse(Genetics == 'FUS', 'FUS-ALS', Genetics)) %>%
  mutate(Genetics=str_replace_all(Genetics, 'ALS with no known mutation', 'sporadic'))


# Filter to cytosolic features
if (cytosolic) {
  # Initialize channel / region IDs
  channels <- sapply(str_split(colnames(x), '\\.\\.'), '[', 3)
  
  regions <- sapply(str_split(colnames(x), '\\.\\.'), '[', 2)
  regions <- lapply(regions, str_split, pattern='_') 
  regions <- sapply(regions, function(z) tail(z[[1]], 1))
  
  # Initialize FUS nuclear features
  channels <- str_split(channels, '')
  regions <- str_split(regions, '')
  
  nuc.fus <- mapply(function(ch, rg) {
    any(ch == 2 & rg == 'D') |   any(ch == 2 & rg == 'C')
  }, channels, regions)
  
  nuc.fus[is.na(nuc.fus)] <- FALSE
  
  # Drop nuclear features
  x <- dplyr::select_if(x, !nuc.fus)
}

#' ## Predictive analysis
#' Here we address the question of whether FUS-ALS cell lines can be separated 
#' from Healthy cell lines in a single screen. Toward this end, we train binary 
#' classifiers to discriminate between FUS-ALS and Healthy `r type`s. 
#' 
#' To assess the generalizability of our predictive models to new cell lines, we 
#' consider a leave one cell line out strategy. Specifically, we remove a single 
#' cell line, train models on the remaining cell lines, and evaluate predictions 
#' on the held-out line.
#+ leave_cell_out

# Fit model to select data
set.seed(47)

if (!file.exists(file.path(output.fig.dir, 'model.Rdata'))) {
  cell.models <- cell_fit(x, model, model_predict)
  save(file=file.path(output.fig.dir, 'model.Rdata'), cell.models)
} else {
  load(file.path(output.fig.dir, 'model.Rdata'))
}

#' ## Fig 2A
#' Model predictions between Healthy and FUS-ALS cell lines. Predictions scores 
#' for a given cell line are computed as the average model prediction over all 
#' cells from a given cell line. Supplementary figures report prediction scores
#' for WT vs. sporadic cell lines and by metadata features.
#+ fig2a

# Aggregate FUS v. Healthy predictions by cell line
predictions <- cell.models$ypred %>% 
  group_by(CellLine, Ytest) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine')

################################################################################
# FUS-ALS vs. WT, Fig. 2a
################################################################################
# Compute test statistic for difference in prediction scores
x.wt.fus <- filter(predictions, Genetics != 'sporadic')
fit <- lm(Ypred ~ Genetics + Sex + `Origin Institution` + CellLinePassageNumber, data=x.wt.fus)
pval.lm <- summary(fit)$coefficients['GeneticsHealthy', 'Pr(>|t|)']
print(pval.lm)

yfus <- predictions$Ypred[predictions$Genetics == 'FUS-ALS']
ywt <- predictions$Ypred[predictions$Genetics == 'Healthy']
pval.sr <- wilcox.test(yfus, ywt, alternative='greater')$p.value
print(pval.sr)

pval.lab <- label_pval(pval.sr)
ymax <- max(predictions$Ypred) + 1e-2

p <- predictions %>%
  filter(Genetics != 'sporadic') %>%
  ggplot(aes(x=reorder(Genetics, Ypred, median), y=Ypred)) +
  geom_boxplot(aes(fill=Genetics), outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, size=3, alpha=0.6) +
  geom_text(x=1.5, y=ymax + 5e-2, label=pval.lab, size=8) +
  geom_segment(x=1, xend=2, y=ymax, yend=ymax, col='#0B4668', linewidth=1) +
  ylab('FUS-ALS phenotype score') +
  xlab(NULL) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=col.pal.g) +
  scale_color_manual(values=col.pal.g) +
  ylim(c(0, 1))

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig2a.pdf'), h=8, w=4)
plot(p)
if (save.fig) dev.off()

################################################################################
# FUS-ALS vs. WT, Fig. 4a
################################################################################
# Compute test statistic for difference in prediction scores
x.wt.sp <- filter(predictions, Genetics != 'FUS-ALS')
fit <- lm(Ypred ~ Genetics + Sex + `Origin Institution` + CellLinePassageNumber, data=x.wt.sp)
pval.lm <- summary(fit)$coefficients['Geneticssporadic', 'Pr(>|t|)']

ysporadic <- predictions$Ypred[predictions$Genetics == 'sporadic']
pval.sr <- wilcox.test(ysporadic, ywt, alternative='greater')$p.value

pval.lab.s <- label_pval(pval.lm)
ymax.s <- max(predictions$Ypred[predictions$Genetics == 'sporadic']) + 5e-2

p <- predictions %>%
  ggplot(aes(x=reorder(Genetics, Ypred, median), y=Ypred)) +
  geom_boxplot(aes(fill=Genetics), outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, size=3, alpha=0.6) +
  geom_text(x=1.5, y=ymax + 5e-2, label=pval.lab, size=8) +
  geom_segment(x=1, xend=3, y=ymax, yend=ymax, col='#0B4668', size=1) +
  geom_text(x=1.5, y=ymax.s + 5e-2, label=pval.lab.s, size=8) +
  geom_segment(x=1, xend=2, y=ymax.s, yend=ymax.s, col='#0B4668', size=1) +
  ylab('FUS-ALS phenotype score') +
  xlab(NULL) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=col.pal.g) +
  scale_color_manual(values=col.pal.g) +
  ylim(0:1)

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig4a.pdf'), h=8, w=4)
plot(p)
if (save.fig) dev.off()

################################################################################
# Predictions vs. metadata features
################################################################################
p1 <- predictions %>%
  filter(Genetics == 'Healthy') %>%
  ggplot(aes(x=Sex, y=Ypred, fill=Sex)) +
  geom_boxplot(outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, size=3, alpha=0.6) +
  xlab('Sex') +
  theme(legend.position='none') +
  ylab('FUS-ALS phenotype score') +
  scale_fill_hue(l=55) +
  labs(fill=NULL) +
  ylim(0:1)

p2 <- predictions %>%
  filter(Genetics == 'Healthy') %>%
  ggplot(aes(x=`Origin Institution`, y=Ypred, fill=`Origin Institution`)) +
  geom_boxplot(outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, size=3, alpha=0.6) +
  xlab('Origin institution') +
  theme(legend.position='none') +
  ylab(NULL) +
  scale_fill_jama() +
  labs(fill=NULL) +
  ylim(0:1)

p3 <- predictions %>%
  filter(Genetics == 'Healthy') %>%
  ggplot(aes(x=`CellLinePassageNumber`, y=Ypred)) +
  geom_smooth(method='lm', se=FALSE) +
  geom_point(size=3) +
  theme(legend.position='none') +
  ylab(NULL) +
  xlab('Cell Line Passage Number')

if (save.fig) pdf(file=file.path(output.fig.dir, 'supp_metadata.pdf'), h=8, w=24)
gridExtra::grid.arrange(p1, p2, p3, nrow=1)
if (save.fig) dev.off()

# Aggregate by cell line / plate
predictions <- cell.models$ypred %>% 
  group_by(CellLine, Ytest, PlateID) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine') %>%
  dplyr::rename(Figure=`Figure Names`) %>%
  mutate(Plate=as.numeric(as.factor(PlateID))) %>%
  mutate(Plate=str_c('Plate ', Plate))

p4 <- predictions %>%
  filter(Genetics == 'Healthy') %>%
  ggplot(aes(x=Plate, y=Ypred, fill=PlateID)) +
  geom_boxplot(outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, size=3, alpha=0.6) +
  theme(legend.position='none') +
  ylab('FUS-ALS phenotype score') +
  xlab(NULL) +
  scale_fill_hue(l=55) +
  labs(fill=NULL) +
  ylim(0:1)

# Aggregate by cell line / plate
predictions <- cell.models$ypred %>% 
  group_by(CellLine, Ytest, PlateID, WellID) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine') %>%
  dplyr::rename(Figure=`Figure Names`) %>%
  mutate(Row=str_remove_all(WellID, '-.*$')) %>%
  mutate(Row=as.numeric(Row)) %>%
  mutate(Col=str_remove_all(WellID, '^.*-')) %>%
  mutate(Col=as.numeric(Col)) %>%
  mutate(Plate=as.numeric(as.factor(PlateID))) %>%
  mutate(Plate=str_c('Plate ', Plate))

p5 <- predictions %>%
  mutate(Ypred=ifelse(Genetics == 'Healthy', Ypred, NA)) %>%
  ggplot(aes(x=Col, y=Row, fill=Ypred)) +
  geom_tile() +
  ylab('Plate Column') + 
  xlab('Plate Row') +
  scale_fill_viridis_c() +
  facet_wrap(~Plate) +
  labs(fill='FUS-ALS phenotype score')


if (save.fig) pdf(file=file.path(output.fig.dir, 'supp_plate.pdf'), h=8, w=8)
plot(p4)
if (save.fig) dev.off()

if (save.fig) pdf(file=file.path(output.fig.dir, 'supp_well.pdf'), h=12, w=24)
plot(p5)
if (save.fig) dev.off()

#' ## Fig 2B
#' Variability of model predictions across well replicates. Prediction scores 
#' for a given replicate are computed as the average model prediction across all 
#' cells from that well.
#+ fig2b

# Aggregate predictions by well replicate
predictions <- cell.models$ypred %>% 
  group_by(CellLine, WellID, Ytest, PlateID) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine') %>%
  dplyr::rename(Figure=`Figure Names`)

################################################################################
# FUS-ALS vs. WT, Fig. 2b
################################################################################
predictions.agg <- filter(predictions, Genetics %in% c('FUS-ALS', 'Healthy')) %>%
  group_by(CellLine, Genetics) %>%
  summarize(Ypred=median(Ypred), .groups='drop') %>%
  arrange(Ypred) %>%
  mutate(Face=ifelse(Genetics == 'Healthy', 'bold', 'plain')) %>%
  mutate(Face=ifelse(CellLine %in% c('NS020', 'NS044'), 'bold', Face))

# Plot predictions, FUS v. Healthy
p <- filter(predictions, Genetics %in% c('FUS-ALS', 'Healthy')) %>%
  ggplot(aes(x=reorder(Figure, Ypred, median), y=Ypred)) +
  geom_boxplot(aes(fill=Genetics), outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, alpha=0.6) +
  theme(axis.text.x=element_text(angle=90, face=predictions.agg$Face)) +
  theme(legend.position='none') +
  xlab(NULL) +
  ylab('FUS-ALS phenotype score') +
  scale_color_manual(values=col.pal) +
  scale_fill_manual(values=col.pal) +
  ylim(c(0, 1))

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig2b.pdf'), h=8, w=16)
plot(p)
if (save.fig) dev.off()

################################################################################
# FUS-ALS vs. WT, Fig. 4b
################################################################################
# plot predictions, sporadic v. Healthy
predictions.agg <- filter(predictions, Genetics %in% c('sporadic', 'Healthy')) %>%
  group_by(CellLine, Genetics) %>%
  summarize(Ypred=median(Ypred), .groups='drop') %>%
  arrange(Ypred) %>%
  mutate(Face=ifelse(Genetics == 'Healthy', 'bold', 'plain'))

p <- filter(predictions, Genetics %in% c('sporadic', 'Healthy')) %>%
  ggplot(aes(x=reorder(Figure, Ypred, median), y=Ypred)) +
  geom_boxplot(aes(fill=Genetics), outlier.shape=NA, coef=NULL) +
  geom_jitter(width=0.2, alpha=0.6) +
  theme(axis.text.x=element_text(angle=90, face=predictions.agg$Face)) +
  theme(legend.position='none') +
  xlab(NULL) +
  ylab('FUS-ALS phenotype') +
  scale_color_manual(values=col.pal) +
  scale_fill_manual(values=col.pal) +
  ylim(c(0, 1)) +
  geom_vline(xintercept=18.5, lty=2) +
  geom_text(x=9.25, y=0.9, label='sporadic-', size=10, col=col.pal['sporadic']) +
  geom_text(x=21.25, y=0.9, label='sporadic+', size=10, col=col.pal['sporadic']) 

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig4b.pdf'), h=8, w=16)
plot(p)
if (save.fig) dev.off()


#' ## Fig 2C
#' ROC curve for Healthy / FUS classification performance
#' roc
# Aggregate predictions by single cell
predictions.cell <- filter(cell.models$ypred, Genetics %in% c('FUS-ALS', 'Healthy'))

# Aggregate predictions by single well
predictions.well <- filter(cell.models$ypred, Genetics %in% c('FUS-ALS', 'Healthy')) %>%
  group_by(CellLine, WellID, Ytest, PlateID) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine')

# Aggregate predictions by cell line
predictions.line <- filter(cell.models$ypred, Genetics %in% c('FUS-ALS', 'Healthy')) %>%
  group_by(CellLine, Ytest) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine')

# Compute ROC curves
roc.cell <- roc.curve(
  scores.class0=predictions.cell$Ypred, 
  weights.class0=predictions.cell$Ytest,
  curve=TRUE
)

roc.well <- roc.curve(
  scores.class0=predictions.well$Ypred, 
  weights.class0=predictions.well$Ytest,
  curve=TRUE
)

roc.line <- roc.curve(
  scores.class0=predictions.line$Ypred, 
  weights.class0=predictions.line$Ytest,
  curve=TRUE
)

# Initialize type labels
type.levels <- c('Single cell', 'Well replicate', 'Cell line')
aucs <- c(roc.cell$auc, roc.well$auc, roc.line$auc) %>% round(2)
type.levels <- str_c(type.levels, ', AUROC = ', aucs)

# Plot ROC curves
xplot = rbind(
  data.frame(roc.cell$curve[,1:2], Type=type.levels[1]),
  data.frame(roc.well$curve[,1:2], Type=type.levels[2]),
  data.frame(roc.line$curve[,1:2], Type=type.levels[3])
) 

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig2c.pdf'), h=8, w=8)
mutate(xplot, Type=factor(Type, levels=type.levels)) %>%
  ggplot(aes(x=X1, y=X2, col=Type)) +
  geom_line(size=1.5, alpha=0.8) +
  geom_abline(col='grey', lty=2) +
  scale_color_jama() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  theme(legend.position=c(0.75, 0.15)) +
  theme(legend.text=element_text(size=16)) +
  labs(color=NULL)
if (save.fig) dev.off()


#' ## Fig 2D
#' Feature importance in classification models.
#+ feature_importance
# Max-normalize feature importance for each hold-out cell line 
importance <- apply(cell.models$importance, MAR=2, importance_norm) 
features <- rownames(cell.models$importance) %>% str_remove_all('^X\\.\\.')

# Compute maximum importance within feature category, average across models
importance.avg <- reshape2::melt(importance) %>%
  mutate(Feature=features[Var1]) %>%
  dplyr::rename(Model=Var2) %>%
  mutate(Feature=str_remove_all(Feature, '\\.\\.[0-9]*$')) %>%
  group_by(Model, Feature) %>%
  summarize(Importance=max(value), .groups='drop') %>%
  group_by(Feature) %>%
  summarize(SD=sd(Importance), Importance=mean(Importance))

# Plot ranked importance for grouped features
clean_feature <- function(x) {
  # Reformatting for feature name
  xsplit <- str_split(x, '\n ')[[1]]
  feature <- xsplit[1]
  xloc.split <- str_split(xsplit[2], ', ')[[1]] %>% str_split('/')
  
  if (length(xloc.split[[1]]) == 2) {
      region <- str_c(
        str_c(xloc.split[[1]][1], ' ', xloc.split[[2]][1]), ' / ',
        str_c(xloc.split[[1]][2], ' ', xloc.split[[2]][2])
      )
  } else {
      region <- str_c(xloc.split[[1]][1], ' ', xloc.split[[2]][1])
  }
  
  return(str_c(feature, '\n', region))
}

if (save.fig) pdf(file=file.path(output.fig.dir, 'fig2d.pdf'), h=7, w=8)
top_n(importance.avg, 5, Importance) %>%
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
  mutate(Feature=str_replace_all(Feature, 'C', 'cell,')) %>%
  mutate(Feature=str_replace_all(Feature, '_', ' ')) %>%
  mutate(Feature=str_replace_all(Feature, 'Intensity Ratio', 'Intensity Ratio\n')) %>%
  mutate(Feature=str_replace_all(Feature, 'Intensity Region', 'Intensity Region\n')) %>%
  mutate(Feature=str_replace_all(Feature, '  ', ' ')) %>%
  mutate(Feature=str_replace_all(Feature, '\\.\\.', ' ')) %>%
  mutate(Feature=sapply(Feature, clean_feature)) %>%
  ggplot(aes(x=reorder(Feature, Importance), y=Importance)) +
  geom_bar(stat='identity', fill='#0B4668', color='black') +
  geom_errorbar(aes(ymin=Importance, ymax=Importance + SD), width=0.25) +
  theme(axis.text.y=element_text(size=16)) +
  coord_flip() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle(NULL) +
  theme(text=element_text(size=22))
if (save.fig) dev.off()

# Save table of feature importances
clean_region <- function(x) {
  x <- str_replace_all(x, 'C', 'Cell') %>%
    str_replace_all('N', 'Cytoplasm') %>%
    str_replace_all('D', 'DNA')
  
  return(x)
}

clean_marker <- function(x) {
  x <- str_replace_all(x, '2', 'FUS') %>%
    str_replace_all('3', 'EEA1')
  
  return(x)
}

importance.avg <- importance.avg %>%
  mutate(Region=str_remove_all(Feature, '^.*_')) %>%
  mutate(Marker=str_remove_all(Region, '^.*\\.\\.')) %>%
  mutate(Region=str_remove_all(Region, '\\.\\.[0-9]*')) %>%
  mutate(Feature=str_remove_all(Feature, '_(C|D|N)*\\.\\.[0-9]*$')) %>%
  mutate(Region=str_split(Region, '')) %>%
  mutate(Region=sapply(Region, str_c, collapse=' ')) %>%
  mutate(Marker=str_split(Marker, '')) %>%
  mutate(Marker=sapply(Marker, str_c, collapse=' ')) %>%
  mutate(Region=clean_region(Region)) %>%
  mutate(Marker=clean_marker(Marker)) %>%
  arrange(desc(Importance))

fout <- file.path(output.fig.dir, 'importance.csv')
write.csv(file=fout, importance.avg)

#' ## Fig 2F
#' Correlation between model predicted probability FUS and patient age.
#' predictions_v_age
# Aggregate model predictions by cell line
predictions <- cell.models$ypred %>%
  dplyr::rename(Plate=PlateID) %>%
  group_by(CellLine, Ytest, Genetics, WellID, Plate) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')

# Plot age v. predictions by genetics
xcell.s <- dplyr::select(xcell, -Genetics) %>%
  dplyr::rename(Age=`Age of biopsy`) %>%
  dplyr::rename(AgeOfOnset=`Age of onset`) %>%
  dplyr::rename(Figure=`Figure Names`)

xplot <- predictions %>% left_join(xcell.s, by='CellLine')

xplot.avg <- group_by(predictions, CellLine, Genetics) %>%
  summarize(Ypred=mean(Ypred), .groups='drop') %>%
  left_join(xcell.s, by='CellLine') %>%
  mutate(AgeOfOnset=ifelse(AgeOfOnset == 'Control', Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset=ifelse(AgeOfOnset == 'control', Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset=ifelse(AgeOfOnset == '', Age, AgeOfOnset)) %>%
  mutate(AgeOfOnset=as.numeric(AgeOfOnset)) %>%
  mutate(Figure=str_remove_all(Figure, '-.*$')) %>%
  dplyr::rename(Site=`Origin Institution`)

# Compute age/prediction correlations within genetic background
xfus <- filter(xplot.avg, Genetics == 'FUS-ALS')
cor.fus <- cor.test(xfus$Ypred, xfus$Age, alternative='less')

xwt <- filter(xplot.avg, Genetics == 'Healthy')
cor.wt <- cor.test(xwt$Ypred, xwt$Age, alternative='less')

xsporadic <- filter(xplot.avg, Genetics == 'sporadic')
cor.spor <- cor.test(xsporadic$Ypred, xsporadic$Age, alternative='less')

rhos <- round(c(cor.fus$estimate, cor.wt$estimate, cor.spor$estimate), 3)
pvals <- c(cor.fus$p.value, cor.wt$p.value, cor.spor$p.value)
pvals <- pvals * length(pvals)
pval.lab <- sapply(pvals, label_pval)

xplot.text <- data.frame(Genetics=c('FUS-ALS', 'Healthy', 'sporadic')) %>%
  mutate(text=str_c('Cor = ', rhos, pval.lab)) %>%
  mutate(Ypred=c(0.9, 0.85, 0.8)) %>%
  mutate(Age=55)

# Plot predictions v. age, FUS v. Healthy
if (save.fig) pdf(file=file.path(output.fig.dir, 'fig_supp_aob.pdf'), h=10, w=10)
filter(xplot.avg, Genetics %in% c('Healthy', 'FUS-ALS')) %>%
  ggplot(aes(x=Age, y=Ypred, col=Genetics)) +
  geom_point(size=3) +
  geom_line(stat="smooth", method="lm", formula=y~x, alpha=0.7, size=1.5) +
  geom_text_repel(aes(label=Figure), size=7, min.segment.length=0.01) +
  geom_text(data=filter(xplot.text, Genetics != 'sporadic'), aes(label=text), size=10) +
  scale_color_manual(values=col.pal.g) +
  ylim(0:1) +
  ylab('FUS-ALS phenotype score') +
  theme(legend.position='none')
if (save.fig) dev.off()

if (save.fig) pdf(file=file.path(output.fig.dir, 'supp_sporadic_age.pdf'), h=10, w=10)
filter(xplot.avg, Genetics %in% c('Healthy', 'sporadic')) %>%
  ggplot(aes(x=Age, y=Ypred, col=Genetics)) +
  geom_point(size=3) +
  geom_line(stat="smooth", method="lm", formula=y~x, alpha=0.7, size=1.5) +
  geom_text_repel(aes(label=Figure), size=7, min.segment.length=0.01) +
  geom_text(data=filter(xplot.text, Genetics != 'FUS-ALS'), aes(label=text), size=10) +
  scale_color_manual(values=col.pal.g) +
  ylim(0:1) +
  ylab('FUS-ALS phenotype score') +
  theme(legend.position='none')
if (save.fig) dev.off()

#' ## Fig S3 age of onset
#' Correlation between model predicted probability FUS and patient age of onset.
#' predictions_v_age_of_onset
# Aggregate model predictions by cell line
# Compute age/prediction correlations within genetic background
cor.fus.aoo <- cor.test(xfus$Ypred, xfus$AgeOfOnset, alternative='less')

rhos <- round(c(cor.fus.aoo$estimate, cor.wt$estimate, cor.spor$estimate), 3)
pvals <- c(cor.fus.aoo$p.value, cor.wt$p.value, cor.spor$p.value)
pvals <- pvals * length(pvals)
pval.lab <- sapply(pvals, label_pval)

xplot.text <- data.frame(Genetics=c('FUS-ALS', 'Healthy', 'sporadic')) %>%
  mutate(text=str_c('Cor = ', rhos, pval.lab)) %>%
  mutate(Ypred=c(0.9, 0.85, 0.8)) %>%
  mutate(AgeOfOnset=55)

# Plot predictions v. age of onset, FUS v. Healthy
if (save.fig) pdf(file=file.path(output.fig.dir, 'fig2g.pdf'), h=12, w=12)
filter(xplot.avg, Genetics %in% c('Healthy', 'FUS-ALS')) %>%
  ggplot(aes(x=AgeOfOnset, y=Ypred, col=Genetics)) +
  geom_point(size=3) +
  geom_line(stat="smooth", method="lm", formula=y~x, alpha=0.7, size=1.5) +
  geom_text_repel(aes(label=Figure), size=7, min.segment.length=0.01) +
  geom_text(data=filter(xplot.text, Genetics != 'sporadic'), aes(label=text), size=10) +
  scale_color_manual(values=col.pal.g) +
  ylim(0:1) +
  ylab('FUS-ALS phenotype score') +
  xlab('Age of onset') +
  theme(legend.position='none')
if (save.fig) dev.off()


################################################################################
# Save table of predictions by cell line
################################################################################
# Predictions by genetics
predictions <- cell.models$ypred %>% 
  group_by(CellLine, Ytest) %>%
  summarize(Ypred=mean(Ypred), .groups='drop')  %>%
  left_join(xcell, by='CellLine')

save(file=file.path(output.fig.dir, 'predictions.Rdata'), predictions)

################################################################################
# Save table of plates used
################################################################################
if (save.plate) {
  plate.list <- unique(predictions$Plate)
  
  write.table(
    file=plate.file, 
    plate.list, 
    append=TRUE, 
    row.names=FALSE, 
    col.names=FALSE, 
    quote=FALSE
  )
  
  # Save plate maps
  plate.layout <- cell.models$ypred %>%
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
