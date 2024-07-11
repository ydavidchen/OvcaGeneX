# Prepare Clinical & Expression Matrix for TCGA-OV

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(data.table)
library(readxl)

COMMON_GENES <- read.table(paste0(paste0(DIR_OUT,"240501_shared_genes_nozero.txt")))$V1

# ---------------------------- Load & clean clinical ----------------------------
load_clin_hgssoc <- function(prim_key="bcr_sample_barcode", sec_key="bcr_patient_barcode") {
  #'@description Loads high grade-stage serous ovarian cancer from TCGA
  
  ## Patient: Liu 2018 survival:
  liu2018 <- read_excel(paste0(DIR,"TCGA_MISC/Liu2018CellTCGACDR.xlsx"), sheet=1, na=c("[Not Available]","[Not Applicable]"))
  liu2018 <- subset(liu2018, type=="OV")
  colnames(liu2018)[colnames(liu2018)=="age_at_initial_pathologic_diagnosis"] <- "AgeAtDx"
  colnames(liu2018)[colnames(liu2018)=="histological_grade"] <- "Grade"
  liu2018$raceWhite <- liu2018$race == "WHITE"
  liu2018 <- liu2018[ , c(sec_key,"AgeAtDx","Grade","raceWhite","tumor_status","OS.time","OS","PFI","PFI.time")]
  
  ## GDAC FireBrowse w/ additional covariates:
  gdac <- read.csv(paste0(DIR, "SynapseTCGAlive/OV/gdac_firebrowse_ov_cleaned.csv"))
  gdac$bcr_patient_barcode <- toupper(gdac$bcr_barcode)
  gdac$Neoadjuvant <- gdac$history_of_neoadjuvant_treatment=="yes"
  gdac$Radiation <- gdac$radiation_therapy=="yes"
  gdac$clinical_stage <- toupper(gsub("stage ", "", gdac$clinical_stage))
  gdac$simple_stage <- gsub("[A-C]", "", gdac$clinical_stage, ignore.case=FALSE)
  gdac <- gdac[ , c(sec_key,"ethnicity","simple_stage","Neoadjuvant","Radiation")]
  
  ## Synapse biosamples:
  samples <- read.csv(paste0(DIR, "SynapseTCGAlive/OV/nationwidechildrens.org_ov_bio.sample.tsv"), sep="\t")
  samples <- subset(samples, sample_type=="Primary Tumor") #exclude blood normals etc.
  samples <- samples[ , c(prim_key,"sample_type")]
  samples[sec_key] <- substr(samples$bcr_sample_barcode, 1, 12) #no bracket use
  
  ## Merge:
  res <- merge(liu2018, gdac, by=sec_key, all=TRUE) #patient-level
  res <- merge(res, samples, by=sec_key, all=TRUE) #add in 1st sample-level mapper
  
  ## Restrict:
  res <- subset(res, simple_stage %in% c("III","IV") & Grade %in% c("G2","G3"))
  return(res)
}

clin_tcga <- load_clin_hgssoc()

# ---------------------------- Load & Clean Expression ----------------------------
cts <- read.csv(paste0(DIR,"SynapseTCGAlive/OV/unc.edu_OV_IlluminaHiSeq_RNASeqV2.geneExp.tsv"), sep="\t", check.names=FALSE)
cts$gene_id <- gsub("\\|.*", "", cts$gene_id)

cts <- subset(cts, gene_id %in% COMMON_GENES)
sum(duplicated(cts$gene_id)) #1

cts <- aggregate(. ~ gene_id, data=cts, FUN=mean, na.rm=TRUE)
rownames(cts) <- cts$gene_id
cts$gene_id <- NULL
cts <- data.matrix(cts)
hist(cts)

anyDuplicated(substr(colnames(cts), 1, 16))
colnames(cts) <- substr(colnames(cts), 1, 16)

## Mutually subset:
cts <- cts[ , colnames(cts) %in% clin_tcga$bcr_sample_barcode]
clin_tcga <- subset(clin_tcga, bcr_sample_barcode %in% colnames(cts))
identical(colnames(cts), clin_tcga$bcr_sample_barcode) #optional check

# ---------------------------- Update Clinical Data w/ Gene Strata ----------------------------
cl_tcga <- stratify_by_gene(cts, GENE, "bcr_sample_barcode")
clin_tcga <- merge(clin_tcga, cl_tcga, by="bcr_sample_barcode")
rownames(clin_tcga) <- NULL

LEV <- c("High","Low")
clin_tcga$Group0 <- factor(clin_tcga$Group0, levels=LEV)
clin_tcga$GroupM <- factor(clin_tcga$GroupM, levels=LEV)

cts <- cts[rownames(cts)!= GENE, ]

# ---------------------------- Export ----------------------------
## 1. Export objects for downstream:
expr_tcga <- apply(cts, 1, StandardScale) #normalize for each GENE! will transpose...
expr_tcga <- t(expr_tcga)
sum(is.na(expr_tcga))
hist(expr_tcga)

expr_tcga <- winsorize(expr_tcga, -6, 6)

save(
  list = c("expr_tcga", "clin_tcga"),
  file = paste0(DIR_OUT, "240504a_tcgaov.RData"),
  compress = TRUE
)

## 2. Export raw counts:
write.table(cts, paste0(DIR_OUT,"240505a_hugo_counts_tcga.tsv"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

## 3. Export log2 counts:
cts <- log2(cts + 1)
hist(cts)

write.table(cts, paste0(DIR_OUT,"240505a_hugo_log2_expr_tcga.tsv"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
