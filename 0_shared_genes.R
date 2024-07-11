# Identify Shared Genes between Cohorts
# Notes:

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## ENSEMBL-indexed combined AOCS+MOCOG expression data:
ENSEMBL <- read.csv(paste0(DIR,PATH_AUS,"ensembl_dict_GSE211669_column1.csv"))
ENSEMBL <- subset(ENSEMBL, chr %in% CHROMS)
ENSEMBL <- subset(ENSEMBL, category == "protein_coding")

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

eAus <- subset(eAus, GeneID %in% ENSEMBL$ensembl)
ENSEMBL <- subset(ENSEMBL, ensembl %in% eAus$GeneID)
stopifnot(identical(eAus$GeneID, ENSEMBL$ensembl)) #checkpoint
eAus$GeneID <- ENSEMBL$HUGO

## HUGO indexed TCGA-OV:
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

## Find common HUGO & export:
gShared <- intersect(unique(eTCGA$gene_id), unique(eAus$GeneID))
length(gShared) #17455
# write.table(gShared, paste0(DIR_OUT,"240501_shared_genes_nozero.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
