# Compute SBS 3 for TCGA from Trinucleotide Matrix
# Notes:
# - Only need to compute for TCGA. AOCS alaready available as Suppl Table 9 of Garsed et al.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(maftools)
library(NMF)

## Trinucleotide Matrix from Alexandrov 2020 Nature:
tcga96 <- read.csv(paste0(DIR, "TCGA_MISC/Alexandrove2020_TCGA/WES_TCGA.96.csv"), check.names=FALSE)
rownames(tcga96) <- paste0(substr(tcga96$Trinucleotide,1,1), "[", tcga96$`Mutation type`, "]", substr(tcga96$Trinucleotide,3,3))
tcga96 <- tcga96[ , grepl("Ovary", colnames(tcga96))]
colnames(tcga96) <- gsub(".*::", "", colnames(tcga96))

## Prep data in required format (see Vignette example):
tcga96 <- t(tcga96) #row=samples, column=trinucleotides
tcga96 <- list(nmf_matrix=tcga96) #required input format
View(tcga96$nmf_matrix)

## Identify optimal number of de novo signatures:
ov_trial <- estimateSignatures(tcga96) #start w/ default of 6
plotCophenetic(res=ov_trial) 

## Extract signature using optimal number:
ov_sig <- extractSignatures(
  mat = tcga96,
  n = 4, #identified w/ Elbow Method
  pConstant = 0.01
)

## Identify extracted signatures matching SBS3:
ov_db_check <- compareSignatures(nmfRes=ov_sig, sig_db="SBS")
# -Comparing against COSMIC signatures
# ...
# --Found Signature_4 most similar to SBS3
# Aetiology: Defects in DNA-DSB repair by HR [cosine-similarity: 0.927]

SBS3_NAME <- "Signature_4" #arbitrary name

res <- t(ov_sig$contributions[SBS3_NAME, , drop=FALSE])
colnames(res) <- "SBS3_contrib"
# write.csv(res, paste0(DIR_OUT,"20240714_tcga_sbs3.csv"), row.names=TRUE, quote=FALSE)
