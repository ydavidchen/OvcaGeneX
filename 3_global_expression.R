# Global Expression by PCA
# Notes:
# - Objective: Investigate relationship btwn Groups vs. global expression summaries

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggplot2)

PROP_SELE <- 0.75 #most variable genes to keep in each cohort

load(paste0(DIR_OUT,"240504a_tcgaov.RData"))
load(paste0(DIR_OUT, "240504b_aocs_mocog.RData"))

## Keep shared most variable genes:
expr_tcga <- select_most_var(expr_tcga, PROP_SELE)
expr_aus <- select_most_var(expr_aus, PROP_SELE)

gShared <- intersect(rownames(expr_tcga), rownames(expr_aus))
length(gShared) #11163
length(gShared) / nrow(expr_tcga) #0.8527884

expr_tcga <- expr_tcga[gShared, ]
expr_aus <- expr_aus[gShared, ]

## PCA:
wrapper_pca <- function(mat, clin=NULL, prim_key=NULL) {
  #'@description Wrapper to systematically run 2-dimensional PCA & downstream stats
  #'@param mat Expression matrix w/ standard format, i.e. row=genes
  #'@describeIn Adapted from Chen 2023 JOMES
  
  pcaObj <- prcomp(t(mat), center=TRUE, scale.=TRUE) #requires columns=features
  
  resPca <- scale(pcaObj$x[ , 1:2]) #columnwise
  resPca <- as.data.frame(resPca)
  
  pcaPerf <- summary(pcaObj)$importance[c(2,3,1), 1:2] #2D performance stats
  explVars <- round(100*c(pcaPerf[1,1], pcaPerf[1,2]), 2)
  print(pcaPerf)
  
  ## Optional Visualizations:
  if(! is.null(clin) & ! is.null(prim_key)) {
    resPca[prim_key] <- rownames(resPca)
    resPca <- merge(resPca, clin[,c(prim_key,"Group0")], by=prim_key)
    
    pltS <- ggplot(resPca, aes(PC1, PC2, color=Group0)) + 
      geom_point(size=3, alpha=0.75) + 
      scale_x_continuous(limits=c(-3,3)) +
      scale_y_continuous(limits=c(-3,3)) +
      xlab(paste0("PC1 (", explVars[1], "% var.)")) +
      ylab(paste0("PC2 (", explVars[2], "% var.)")) +
      geom_vline(xintercept=0, linetype=2) +
      THEME_SCATTER
    
    mResPca <- melt(resPca)
    pltB <- ggplot(mResPca, aes(variable, value, color=Group0)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(position=position_jitterdodge(jitter.width=0.2)) +
      THEME_BOX
    
    print( gridExtra::grid.arrange(grobs=list(pltS,pltB), widths=c(0.6,0.4)) ) 
    print(t.test(PC1 ~ Group0, data=resPca))
    print(t.test(PC2 ~ Group0, data=resPca))
  }
  
  ## Arbitrary grouping by only 1st component: 
  resPca$pcGroup <- ifelse(resPca$PC1 > 0, "A", "B")
  cTab <- table(Expr=resPca$Group0=="High", PC=resPca$pcGroup=="A")
  print(cTab)
  print(fisher.test(cTab))
  
  return(resPca)
}

pca_tcga <- wrapper_pca(expr_tcga, clin_tcga, "bcr_sample_barcode")
pca_aus <- wrapper_pca(expr_aus, clin_aus, "Sample_title")
