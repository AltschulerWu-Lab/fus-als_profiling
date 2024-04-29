#' The following script generates plots for qPCR levels following ASO knockdown
#' of FUS and FUS intensity levels in microscopy screens.
#'
#' Requires input datasets
#'   `data/2022_ASOscreen_v3_qPCR results.xls`
#'   `data/080122_ASOScreenV2_Cell`
#'
#' Karl Kumbier, 12/19/2023
library(data.table)
library(tidyverse)
library(tidytext)
library(ggsci)
library(ggpubr)
library(lmerTest)

theme_set(theme_bw(base_size = 22))

# Set parameters for analysis
save.fig <- TRUE
treatment.levels <- c("untreated", "NTC 10uM", "ASO 1uM", "ASO 10uM")
heat.pal <- c("#FFFFFF", pal_material("blue-grey")(10))
correction <- "BH"

# Set paths to data directories
fig.base.dir <- Sys.getenv("ALS_PAPER")
data.base.dir <- Sys.getenv("ALS_DATA")
figure.dir <- file.path(fig.base.dir, "paper_figures", "Figure_ASO")

source(file.path(fig.base.dir, "paper_figures", "utilities.R"))
source(file.path(fig.base.dir, "paper_figures", "color_palette.R"))

# Set paths to data directories
data.dir <- file.path(data.base.dir, "120423_ASO14", "Cell_Filtered")
library.file <- file.path(fig.base.dir, "data", "Fibroblasts_Library.csv")
qpcr.file <- file.path(fig.base.dir, "data", "20231007_ASO14_qPCRresults.csv")

output.fig.dir <- file.path(figure.dir, "fig/")
dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)

# Load mutation table
xcell <- fread(library.file) %>%
  dplyr::rename(Figure = `Figure Names`) %>%
  dplyr::rename(CellLine = `Cell line`)

mutation.key <- setNames(xcell$Figure, xcell$CellLine)

################################################################################
# FUS qPCR by cell line.
################################################################################
qpcr <- fread(qpcr.file)

p <- qpcr %>%
  filter(ASO != "Un") %>%
  mutate(ASO = str_replace_all(ASO, "1uM", "ASO 1uM")) %>%
  mutate(ASO = str_replace_all(ASO, "10uM", "ASO 10uM")) %>%
  mutate(ASO = str_replace_all(ASO, "NTC", "NTC 10uM")) %>%
  mutate(CellLineL = mutation.key[CellLine]) %>%
  mutate(ASO = factor(ASO, levels = treatment.levels)) %>%
  mutate(CellLineL = str_remove_all(CellLineL, "-.*$")) %>%
  ggplot(aes(x = CellLine, y = `FUS (normalized to NTC)`, fill = ASO)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_fill_manual(values = heat.pal[c(3, 6, 9)]) +
  xlab(NULL) +
  facet_wrap(~CellLineL, nrow = 2, scales = "free_x") +
  theme(axis.text.x = element_blank())

pdf(file.path(output.fig.dir, "fig_supp_qpcr.pdf"), height = 8, width = 16)
plot(p)
dev.off()


################################################################################
# Raw nuclear and cytoplasmic FUS intensities
################################################################################
load(file.path(data.dir, "profiles_raw.Rdata"))
x <- mutate(x, ASO = str_remove_all(ASOtreatment, "EP_")) %>%
  mutate(ASO = str_replace_all(ASO, "ASO", " ASO")) %>%
  mutate(ASO = str_replace_all(ASO, "NTC", "NTC 10uM")) %>%
  mutate(ASO = str_replace_all(ASO, "low ASO", "ASO 1uM")) %>%
  mutate(ASO = str_replace_all(ASO, "high ASO", "ASO 10uM")) %>%
  mutate(ASO = factor(ASO, treatment.levels)) %>%
  filter(Genetics != "sporadic")

params.plot <- list(
  list(
    fselect = as.symbol("X..Iav_in_dna_region..2"),
    main = "Nuclear FUS, average pixel intensity",
    fout = "nuclear_fus.pdf"
  ),
  list(
    fselect = as.symbol("X..Iav_in_non_dna_region..2"),
    main = "Cytosolic FUS, average pixel intensity",
    fout = "cytosolic_fus.pdf"
  ),
  list(
    fselect = as.symbol("X..Imax_in_non_dna_region..3"),
    main = "Cytosolic EEA1, max pixel intensity",
    fout = "cytosolic_eea1.pdf"
  )
)

mselect <- c("Genetics", "ASO")

# Compute well-level median intensities
for (prm in params.plot) {
  fselect <- prm$fselect
  main <- prm$main
  fout <- file.path(output.fig.dir, prm$fout)

  xplot <- group_by(x, ID, ASO, Genetics, CellLine, WellID, PlateID) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup() %>%
    dplyr::select(-ID) %>%
    filter(ASO != "untreated") %>%
    dplyr::rename(Treatment = ASO)

  # Initialize plot range
  xrange <- group_by(xplot, Genetics, Treatment) %>%
    summarize(
      Mean = mean(!!fselect),
      SD = sd(!!fselect)
    ) %>%
    mutate(MaxRange = Mean + SD) %>%
    mutate(MinRange = Mean - SD)


  xrange <- c(min(xrange$MinRange), max(xrange$MaxRange))

  eps <- diff(xrange) / 10

  # Compute significance levels for treatment effect
  lme.test <- lapply(unique(xplot$Genetics), function(g) {
    xg <- filter(xplot, Genetics == g)
    f <- str_c(as.character(fselect), "~ Treatment + (1 | CellLine)")
    mixed <- lmer(as.formula(f), data = xg)
    out <- summary(mixed)$coefficients[3, 5]
    out <- data.frame(Treatment = "ASO 10uM", p = out) %>%
      mutate(Treatment = str_remove_all(Treatment, "Treatment")) %>%
      mutate(Genetics = g)
    return(out)
  })

  test.plot <- rbindlist(lme.test) %>%
    mutate(padj = p.adjust(p, method = correction)) %>%
    mutate(padj = signif(padj, 3)) %>%
    mutate(x = ifelse(Treatment == "ASO 1uM", 1, 2)) %>%
    mutate(x = 2) %>%
    mutate(y = max(xrange) - eps) %>%
    mutate(y = ifelse(Treatment == "ASO 1uM", y + eps / 2, y)) %>%
    mutate(xend = x + 1) %>%
    mutate(yend = y) %>%
    mutate(xtext = 2) %>%
    mutate(ytext = y + eps / 4) %>%
    mutate(x = 1) %>%
    mutate(plab = label_pval(padj))

  p <- xplot %>%
    ggbarplot(
      x = "Treatment",
      fill = "Genetics",
      y = as.character(fselect),
      facet.by = "Genetics",
      merge = TRUE,
      add = "mean_se",
      error.plot = "upper_errorbar"
    ) +
    geom_text(
      data = filter(test.plot, str_detect(Treatment, "1uM")),
      aes(x = xtext, y = ytext, label = plab), size = 8,
    ) +
    geom_text(
      data = filter(test.plot, str_detect(Treatment, "10uM")),
      aes(x = xtext, y = ytext, label = plab), size = 8,
    ) +
    geom_segment(
      data = test.plot,
      aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    ylab(main) +
    xlab(NULL) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = col.pal) +
    coord_cartesian(ylim = xrange) +
    geom_hline(yintercept = 1, col = "grey", lty = 2)  +
    theme(text = element_text(size = 22))


  if (save.fig) pdf(file = fout, height = 8, width = 10)
  plot(p)
  if (save.fig) dev.off()
}
