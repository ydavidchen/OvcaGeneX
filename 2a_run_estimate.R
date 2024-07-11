# Compute Stromal & Immune Components by Gene Strata
# Notes:
# - The estimate pacakge is not built with the intention to use in Live R sessions.
# - Wrapper written to format data for & run ESTIMATE.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(estimate)

run_estimate <- function(path_in, dir_out=paste0(DIR_OUT,"240506_estimate/")) {
  #'@description End-to-end I/O wrapper to load & format data & execute ESTIMATE
  #'@param path_in Input path to log2-transformed RNAseq data
  #'@param dir_out Writes intermediary files which will be deleted AND final output
  
  path_estinp <- paste0(dir_out, "estInput.tsv")
  path_estout <- paste0(dir_out, "estOutput.tsv")
  
  ## Format RNA-seq data:
  expr <- read.table(path_in, check.names=FALSE)
  expr <- cbind(NAME=rownames(expr), Description=rownames(expr), expr)
  rownames(expr) <- NULL
  
  write.table(expr, path_estinp, sep="\t", row.names=FALSE, quote=FALSE) #intermediate file
  
  ## Core algorithm w/ intermediate files as I/O:
  estimateScore(input.ds=path_estinp, output.ds=path_estout, platform="illumina")
  
  file.remove(path_estinp)
  
  ## Format output:
  mEst <- read.table(path_estout, skip=3)
  rownames(mEst) <- mEst$V1
  mEst$V1 <- mEst$V2 <- NULL
  mEst <- t(mEst)
  rownames(mEst) <- colnames(expr)[-c(1:2)]
  
  file.remove(path_estout)
  
  return(mEst)
}

## Main:
est_tcga <- run_estimate(paste0(DIR_OUT, "240505a_hugo_log2_expr_tcga.tsv"))
est_aus <- run_estimate(paste0(DIR_OUT, "240505b_hugo_expr_log2count_aocs.tsv")) 

save(
  list = c("est_tcga", "est_aus"), 
  file = paste0(DIR_OUT, "240506_estimate.RData"),
  compress = TRUE
)
