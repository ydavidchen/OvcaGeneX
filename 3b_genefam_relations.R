# Relationship w/ Gene Family

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)

load(paste0(DIR_OUT,"240504a_tcgaov.RData"))
load(paste0(DIR_OUT,"240504b_aocs_mocog.RData"))

## Linear transformation for heatmap visualization:
expr_tcga <- t(apply(expr_tcga[GFAM, ], 1, MinMaxScale))
expr_aus <- t(apply(expr_aus[GFAM, ], 1, MinMaxScale))

## Heatmap tracking bar annotation: 
HM_COLS <- list(Group=c(High="black", Low="lightgray"))
HM_COLS[["Age60Up"]] <- HM_COLS[["StageIV"]] <- BINARY_COLORS

hm_samps_tcga <- data.frame(
  row.names = clin_tcga$bcr_sample_barcode,
  Group = clin_tcga$Group,
  Age60Up = ifelse(clin_tcga$AgeAtDx >= 60, "Yes", "No"),
  StageIV = ifelse(clin_tcga$simple_stage=="IV", "Yes", "No")
)

hm_samps_aocs <- data.frame(
  row.names = clin_aus$Sample_title,
  Group = clin_aus$Group,
  Age60Up = ifelse(clin_aus$AgeAtDx >= 60, "Yes", "No"),
  StageIV = ifelse(clin_aus$simple_stage=="IV", "Yes", "No")
)

## Heatmaps:
ph_tcga <- pheatmap(
  expr_tcga,
  cutree_cols = 2, 
  annotation_col = hm_samps_tcga,
  annotation_colors = HM_COLS,
  color = HEAT_COLS,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  clustering_distance_cols = CL_PARAMS[2],
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "TCGA"
)

ph_aocs <- pheatmap(
  expr_aus,
  cutree_cols = 2,
  annotation_col = hm_samps_aocs,
  annotation_colors = HM_COLS,
  color = HEAT_COLS,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  clustering_distance_cols = CL_PARAMS[2],
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "AOCS"
)

# --------------------------- Extract Clusters ---------------------------
getClusters <- function(hcObj, num_cl=2, rowname=NULL) {
  #'@param hcObj hclust object e.g. <pheatmap>$tree_col
  #'@references Chen Y 2024 JMCCP; Chen Y 2022 Med Res Arch
  res_df <- data.frame(cutree(hcObj, k=num_cl))
  colnames(res_df) <- "Cluster"
  if(! is.null(rowname)) {
    res_df[rowname] <- rownames(res_df)
    rownames(res_df) <- NULL
    res_df <- res_df[ , c(2,1)]
  }
  return(res_df)
}

cl_tcga <- getClusters(ph_tcga$tree_col, rowname="bcr_sample_barcode")
cl_aocs <- getClusters(ph_aocs$tree_col, rowname="Sample_title")

clin_tcga <- merge(cl_tcga, clin_tcga, by="bcr_sample_barcode")
table(clin_tcga$Group, clin_tcga$Cluster)

clin_aus <- merge(cl_aocs, clin_aus, by="Sample_title")
table(clin_aus$Group, clin_aus$Cluster)
