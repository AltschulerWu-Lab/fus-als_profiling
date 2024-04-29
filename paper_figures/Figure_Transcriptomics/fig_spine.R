################################################################################
#+ setup, echo=FALSE, warning=FALSE, message=FALSE
################################################################################
library(data.table)
library(tidyverse)
library(Matrix)
library(ggsci)
library(PRROC)
library(ggh4x)
library(parallel)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
theme_set(theme_bw(base_size = 22))

# Paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.dir <- file.path(fig.base.dir, "data/050322_RNAseq/")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_Transcriptomics")
pathway.file <- file.path(data.dir, "reactome_pathway.Rdata")
hierarchy.file <- file.path(fig.base.dir, "data", "reactome_hierarchy.Rdata")
output.fig.dir <- file.path(figure.dir, "fig")

disc.col.pal <- pal_npg()(10)[-8]
save.fig <- FALSE
n.null <- 100000
p.thresh <- 0

load(file.path(output.fig.dir, "spine_models.Rdata"))
load(hierarchy.file)

source(file.path(fig.base.dir, "paper_figures", "color_palette.R"))

clean_name <- function(x, n = 10) {
  # Prettify pathway names for figure axes
  x <- str_split(x, " ")[[1]]
  #if (length(x) > 12) n <- 6
  nbreak <- floor(length(x) / n)
  
  # adjust n for roughly equal words per line
  n <- ceiling(length(x) / (nbreak + 1))
  
  if (nbreak > 0) {
    for (i in 1:nbreak) {
      if ((i * n) > length(x)) break
      x[i * n] <- str_c(x[i * n], "\n")
    }
  }
  
  return(str_c(x, collapse = " "))
}

################################################################################
#+ pathways, fig.height=12, fig.width=16
# Initialize top-level pathways in hierarchy
################################################################################
root.nodes <- filter(nodes, Top)$id
root.children <- lapply(root.nodes, get_children, nodes = nodes, edges = edges)
root.children <- mapply(function(r, ch) c(r, ch), root.nodes, root.children, SIMPLIFY = FALSE)

root.children <- reshape2::melt(root.children) %>%
  dplyr::rename(Root = L1) %>%
  dplyr::rename(Depth = L2) %>%
  dplyr::rename(id = value) %>%
  left_join(nodes, by = "id") %>%
  dplyr::rename(Name = label)

root.key <- filter(root.children, Top) %>%
  dplyr::select(Root, Name) %>%
  dplyr::rename(RootName = Name)

pw.enriched <- filter(pw.aucs, p <= p.thresh) %>%
  left_join(root.children, by = "Name") %>%
  left_join(root.key, by = "Root") %>%
  group_by(RootName) %>%
  mutate(Count = n()) %>%
  mutate(RootName = ifelse(Count < 3, "Other", RootName)) %>%
  ungroup() %>%
  filter(!(Name == "Retinoid metabolism and transport" & RootName == "Other"))

# Plot pathways with enrichment
limits <- c(0.25, 1)
breaks <- seq(0.25, 1, by=0.25)

p1 <- filter(pw.enriched, p == 0) %>%
  mutate(Name=str_sub(Name, 1, 50)) %>%
  ggplot(aes(x=AUC, y=RootName)) +
  #geom_boxplot() +
  geom_jitter(height=0.25, width=0) +
  ylab(NULL) +
  scale_x_continuous(limits=limits, breaks=breaks) +
  xlab(NULL) +
  theme(axis.text.x = element_blank()) +
  labs(fill=NULL)

p2 <- ggplot(null.accuracy) +
  geom_histogram(aes(x=AUC)) +
  scale_x_continuous(limits=limits, breaks=breaks) +
  ylab(NULL) +
  xlab('AUROC')

# build the plots 
p1.common <- ggplot_gtable(ggplot_build(p1))
p2.common <- ggplot_gtable(ggplot_build(p2))
p2.common$widths <- p1.common$widths

fout <- file.path(output.fig.dir, 'spine_auroc.pdf')
pdf(file=fout, height=12, width=8)
grid.arrange(p1.common, p2.common, nrow=2, heights=c(5, 1))
dev.off()


# Save Rdata object
fout <- file.path(output.fig.dir, "spine_models.Rdata")
save(file = fout, pw.aucs, pw.models, null.pred, null.accuracy, pw.enriched)

# Save table of pathway enrichment
fout <- file.path(output.fig.dir, "spine_enrichment_table.csv")
xtab <- dplyr::select(pw.enriched, Name, RootName, Size, AUC, p)
write.csv(file = fout, xtab, row.names=FALSE)
