#' The following generates pathway enrichment figures comparing FUS-ALS, 
#' sporadics, and spinal cord samples. To run pathway enrichment analyses, see
#' `deseq_spine.R` and `deseq_fibro.R`.
#' 
#' Karl Kumbier 3/8/2023
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(ggsci)
library(hrbrthemes)
library(patchwork)
library(dendextend)
options(stringsAsFactors=FALSE)

theme_set(
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24,
    strip_text_size=24, 
    axis_text_size=22,
    base_size=22,
    axis=TRUE,
    axis_col='#000000',
    base_family='sans',
    grid=FALSE
  )
)

################################################################################
# Load data and initialize analysis paramters
################################################################################
fig.base.dir <- Sys.getenv('ALS_PAPER')
data.dir <- file.path(fig.base.dir, 'data/050322_RNAseq/')
figure.dir <- file.path(fig.base.dir, 'paper_figures', 'Figure_Transcriptomics')
pathway.file <- file.path(data.dir, 'reactome_pathway.Rdata')

output.fig.dir <- file.path(figure.dir, 'fig')
dir.create(output.fig.dir, showWarnings=FALSE)

source(file.path(fig.base.dir, 'paper_figures/color_palette.R'))
source(file.path(figure.dir, 'utilities.R'))

load(file.path(output.fig.dir, 'fibro_gsea.Rdata'))
load(file.path(output.fig.dir, 'spine_gsea.Rdata'))

heat.pal <- pal_material('blue-grey')(7)
save.fig <- TRUE

################################################################################
# Pathway enrichment for FUS/sporadic fibroblasts, spinal cord
################################################################################
# Load gene sets associated with each pathway annotation
load(pathway.file)
reactome.genes <- reactome.pathways$Genes
names(reactome.genes) <- reactome.pathways$Name

# Initialize set of pathways to visualize
types <- c('FUS-ALS', 'sporadic+', 'sporadic-')

gsea.fibro <- rbind(
  mutate(gsea.fus$pw.full, Type='FUS-ALS'),
  mutate(gsea.sph$pw.full, Type='sporadic+'),
  mutate(gsea.spl$pw.full, Type='sporadic-')
)
  
gsea.spine <- mutate(gsea.spine$pw.full, Type='Spinal cord')

pathways <- c(gsea.fibro$pathway, gsea.spine$pathway)
pathways <- unique(pathways)

# Compute jaccard similarity scores for pathway gene sets
pathway.sim <- sapply(pathways, function(p1) {
  sapply(pathways, function(p2) {
    jaccard(reactome.genes[[p1]], reactome.genes[[p2]])
  })
})

# Plot pathway clusters
clusters <- hclust(dist(1 - pathway.sim))
dend <- as.dendrogram(clusters) %>% set("branches_lwd", 2) %>% as.ggdend

if (save.fig) pdf(file.path(output.fig.dir, 'pathway_dend.pdf'), h=8, w=24)
ggplot(dend)
if (save.fig) dev.off()


# Initialize pathway enrichment table
types <- c(types, 'Spinal cord')
types <- types[c(4, 1, 2, 3)]

xplot.enrich <- rbind(gsea.fibro, gsea.spine) %>%
  dplyr::rename(Pathway=pathway) %>%
  mutate(logp=-log10(padj)) %>%
  mutate(logp=ifelse(logp > 4, 4, logp)) %>%
  mutate(Pathway=str_remove_all(Pathway, 'Nonsense Mediated Decay ')) %>%
  mutate(Pathway=str_remove_all(Pathway, 'Nonsense-Mediated Decay ')) %>%
  mutate(Pathway=str_replace_all(Pathway, '\\(NMD\\)', 'NMD')) %>%
  mutate(Pathway=str_remove_all(Pathway, 'Exon Junction Complex ')) %>%
  mutate(Pathway=str_replace_all(Pathway, '\\(EJC\\)', 'EJC')) %>%
  mutate(Pathway=str_remove_all(Pathway, '\\ *$')) %>%
  mutate(Pathway=str_remove_all(Pathway, ', and subsequent binding to 43S')) %>%
  mutate(Type=factor(Type, levels=types))

# Add low sporadics to table for figure
xplot.enrich.s <- xplot.enrich[1,]
xplot.enrich.s$Type <- "sporadic-"
xplot.enrich.s$logp <- NA

pathways <- pathways %>%
  str_remove_all('Nonsense Mediated Decay ') %>%
  str_remove_all('Nonsense-Mediated Decay ') %>%
  str_replace_all('\\(NMD\\)', 'NMD') %>%
  str_remove_all('Exon Junction Complex ') %>%
  str_replace_all('\\(EJC\\)', 'EJC') %>%
  str_remove_all('\\ *$') %>%
  str_remove_all(', and subsequent binding to 43S')
  
pw.order <- pathways[clusters$order]

if (save.fig) pdf(file.path(output.fig.dir, 'pathway_enrich.pdf'), h=22, w=24)
cols <- c('#C79652', col.pal['FUS-ALS'], col.pal['sporadic'], col.pal['sporadic'])
rbind(xplot.enrich, xplot.enrich.s) %>%
  mutate(Pathway=factor(Pathway, levels=pw.order)) %>%
  mutate(logp=ifelse(logp > 3, 3, logp)) %>%
  ggplot(aes(x=Type, y=Pathway, size=logp, col=logp)) + 
  geom_point() +
  scale_color_gradientn(colours=heat.pal) +
  scale_size(range=c(5, 15)) +
  guides(color=guide_legend(), size=guide_legend()) +
  theme(legend.text=element_text(size=28), legend.title=element_text(size=28)) +
  ylab(NULL) +
  xlab(NULL) +
  labs(colour="-log10 p-value", size="-log10 p-value") +
  theme(axis.text.x=element_text(size=36, angle=70, hjust=-0.025, color=cols)) +
  theme(axis.text.y=element_text(size=32)) +
  theme(plot.margin = unit(c(1,1,1,4), "cm")) +
  scale_x_discrete(position = "top")
if (save.fig) dev.off()
