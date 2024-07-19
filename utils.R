# Utility Module: Constants & Functions
# Notes:

library(data.table)
library(ggplot2)

## Paths ***** MASKED *****:
DIR <- "***** MASKED *****"
DIR_OUT <- "***** MASKED *****"
PATH_AUS <- "***** MASKED *****"
PATH_CSC <- "***** MASKED *****"

GENE <- "***** MASKED *****"
ENTREZ <- "***** MASKED *****" #ncbi.nlm.nih.gov/gene/
HG19_COORDS <- "***** MASKED *****"
GFAM <- c("***** MASKED *****")

STRNAS <- c("", "NA", "n/a", "NaN")
CHROMS <- 1:22
CL_PARAMS <- c("ward.D", "euclidean")
ADMIN_CENSOR <- 36 #mo
DAYSPERMO <- 30.437 #conversion factor

## Helper methods:
StandardScale <- function(vec) {
  #'@usage Use directly or as callback: apply(<df>, <margin>, FUN=StandardScale)
  sig <- sd(vec, na.rm=TRUE)
  if(sig == 0) return(vec)
  mu <- mean(vec, na.rm=TRUE)
  return((vec-mu) / sig)
}

MinMaxScale <- function(vec) {
  xmin <- min(vec, na.rm=TRUE)
  xmax <- max(vec, na.rm=TRUE)
  if(xmax-xmin == 0) return(rep(0, length(vec)))
  return((vec- xmin)/(xmax-xmin))
}

SelectMostVar <- function(mat, size) {
  #'@param size Number or proportion to select
  if(size < 1) size <- as.integer(nrow(mat) * size)
  sele <- order(matrixStats::rowVars(mat), decreasing=TRUE)[1:size]
  return(mat[sele, ])
}

drop_flat_vals <- function(mat, val=0) return(mat[rowSums(mat==val) < ncol(mat), ])

drop_uniform_var <- function(mat) mat[matrixStats::rowVars(mat) > 0, ]

winsorize <- function(mat, lower=NA, upper=NA) {
  #'@description Constrains extreme values in matrix `mat`
  #'@describeIn Chen 2023 JOMES & Chen 2024 JMCCPL
  if(! is.na(lower)) mat[mat < lower] <- lower
  if(! is.na(upper)) mat[mat > upper] <- upper
  return(mat)
}

build_gene_gr <- function(chr, start, end) {
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start=start, end=end))
  return(gr)
}

## Data Visualization:
BINARY_COLORS <- c(Yes="black", No="lightgray")
HEAT_COLS <- colorRampPalette(c("blue","lightgray","red"))(1024)
HEAT_COLS_YB <- colorRampPalette(c("yellow","black","blue"))(1024)
HEAT_COLS_BW <- colorRampPalette(c("gray95","gray9"))(1024)
SURV_COLS <- c("blue","red")

THEME_SCATTER <- theme_gray() + 
  theme(axis.text.x=element_text(size=10,color="black"), axis.title.x=element_text(size=15,color="black"),
        axis.text.y=element_text(size=10,color="black"), axis.title.y=element_text(size=15,color="black"),
        strip.text.x=element_text(size=15,color="black"),
        legend.position="top", legend.text=element_text(size=10,color="black"))

THEME_BOX <- theme_bw() +
  theme(axis.text.x=element_text(size=12,color="black"), axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=15,color="black"),
        strip.text.x=element_text(size=15,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black"))

THEME_BAR <- theme_classic() +
  theme(axis.text.x=element_text(size=10,color="black",angle=90), axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=15,color="black"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=15), legend.text=element_text(size=12,color="black"))

THEME_SURV <- theme_classic() +
  theme(axis.text=element_text(size=20,color="black"), 
        axis.title=element_text(size=20,color="black"),
        title=element_text(size=20,color="black"),
        legend.title=element_text(size=16,color="black",face="bold"), legend.text=element_text(size=16,color="black",face="bold"), 
        legend.box.background=element_rect(colour="black",size=2),
        text=element_text(size=20),
        strip.text.x=element_text(size=20,colour="black",face="bold"))

