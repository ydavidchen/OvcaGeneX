# Exploration / Visualization: Co-expression w/ Candidate Genes
# Notes:
# - This script is for convenient interactive EDA (for discovery purposes) only. Not for final pub

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggplot2)

CANDIDATES <- "YOUR GENES HERE in format c(<gene1>,...,<genek>)"

eda_wrapper <- function(mat, clin, prim_key, g, testByGroup=FALSE) {
  if(g %in% rownames(mat)) {
    df <- as.data.frame(t(mat[g, , drop=FALSE]))
    colnames(df) <- "Gene"
    df <- merge(df, clin, by.x="row.names", by.y=prim_key)
    
    if(testByGroup) print(t.test(Gene ~ Group0, data=df))
    
    plt_box <- ggplot(df, aes(Group0, Gene)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.1) +
      ylab(g) +
      THEME_BOX
    
    plt_scatter <- ggplot(df, aes(StandardScale(log2(mRNA+1)), Gene)) +
      geom_point() +
      geom_smooth(method="lm") +
      labs(x=GENE, y=g) + 
      THEME_SCATTER
    
    gridExtra::grid.arrange(grobs=list(plt_box, plt_scatter), ncol=2)
  }
  else
    paste(g, "NOT available!")
}

# ------------------------------ TCGA ------------------------------
load(paste0(DIR_OUT,"240504a_tcgaov.RData"))
expr_tcga <- winsorize(expr_tcga, -5, 5)

for(g in CANDIDATES) eda_wrapper(expr_tcga, clin_tcga, "bcr_sample_barcode", g, TRUE)

# ------------------------------ AOCS ------------------------------
load(paste0(DIR_OUT, "240504b_aocs_mocog.RData"))
expr_aus <- winsorize(expr_aus, -5, 5)

for(g in CANDIDATES) eda_wrapper(expr_aus, clin_aus, "Sample_title", g, TRUE)
