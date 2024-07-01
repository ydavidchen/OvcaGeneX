# Identify Shared Genes between Cohorts
# Notes:

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

CHROMS <- c(1:22,"X")

## TCGA-OV
eTCGA <- read.csv(paste0(DIR,"SynapseTCGAlive/OV/unc.edu_OV_IlluminaHiSeq_RNASeqV2.geneExp.tsv"), sep="\t")
anyDuplicated(eTCGA$gene_id)

rownames(eTCGA) <- eTCGA$gene_id
eTCGA$gene_id <- NULL
eTCGA <- data.matrix(eTCGA)

nrow(eTCGA) #20531 
eTCGA <- drop_flat_vals(eTCGA)
nrow(eTCGA) #20189

eTCGA <- as.data.frame(eTCGA)
eTCGA$gene_id <- gsub("\\|.*", "", rownames(eTCGA))
length(unique(eTCGA$gene_id)) #20163

## Combined AOCS+MOCOG expression:
eAus <- read.table(paste0(DIR,PATH_AUS,"GSE211669_20_04_2020.HTSeq.all.raw.txt"), header=TRUE)
anyDuplicated(eAus$GeneID)
rownames(eAus) <- eAus$GeneID
eAus$GeneID <- NULL
eAus <- data.matrix(eAus)

nrow(eAus) #57773
eAus <- drop_flat_vals(eAus)
nrow(eAus) #56090

eAus <- as.data.frame(eAus)
eAus$GeneID <- gsub("\\|.*", "", rownames(eAus))
length(unique(eAus$GeneID)) #56090

ENSEMBL <- read.csv(paste0(DIR,PATH_AUS,"ensembl_dict_GSE211669_column1.csv"))
ENSEMBL <- subset(ENSEMBL, chr %in% CHROMS)
ENSEMBL <- subset(ENSEMBL, ensembl %in% eAus$GeneID)
eAus <- subset(eAus, GeneID %in% ENSEMBL$ensembl)
stopifnot(identical(eAus$GeneID, ENSEMBL$ensembl)) #checkpoint
eAus$GeneID <- ENSEMBL$HUGO

## Find common & export:
gShared <- intersect(eTCGA$gene_id, eAus$GeneID)
length(gShared) #18034
# write.table(gShared, paste0(DIR_OUT,"240501_shared_genes_nozero.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
