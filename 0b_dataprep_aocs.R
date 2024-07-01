# Process AOCS+MOCOG (together referred to as "AOCS" or "AUS") RNA-seq Data
# Notes: 

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(data.table)
library(readxl)

COMMON_GENES <- read.table(paste0(paste0(DIR_OUT,"240501_shared_genes_nozero.txt")))$V1

# ---------------------------- AOCS & MOCOG Sample Maps ----------------------------
samps_mocog <- read.csv(paste0(DIR,PATH_AUS,"GSE211669_twoGPLs_curated.csv"), strip.white=TRUE) #all hgsc
samps_mocog$Reanalysis <- gsub("Re-analysis of ", "", samps_mocog$Reanalysis)
table(samps_mocog$SampleType)

samps_aocs <- read.csv(paste0(DIR,PATH_AUS,"GSE209964_series_matrix_curated.csv"), strip.white=TRUE)
table(samps_aocs$SampleType)

sum(samps_mocog$Reanalysis %in% samps_aocs$Accession) #75, as reported in GEO description
samps_aocs$SampleType[! samps_aocs$Accession %in% samps_mocog$Reanalysis]

samps_aocs <- subset(samps_aocs, SampleType=="PrimaryTumour")
samps_aocs <- subset(samps_aocs, CancerType=="High Grade Serous Ovarian Cancer")
samps_mocog <- subset(samps_mocog, SampleType=="Primary Tumor")
samps_mocog <- subset(samps_mocog, ! Reanalysis %in% samps_aocs$Accession)

COL_SAMPS <- c("Accession","Sample_title","case_id","SampleType")
samps_aocs <- samps_aocs[ , COL_SAMPS]
samps_mocog <- samps_mocog[ , COL_SAMPS]

# ---------------------------- Load & Clean Patient Metadata ----------------------------
load_aus_clinical <- function(path=paste0(DIR,PATH_AUS,"clinical_sele.xlsx")) {
  ## Clinical data:
  clin <- read_excel(path, sheet=1, trim_ws=TRUE, na=STRNAS)
  clin <- as.data.frame(clin)
  colnames(clin)[1] <- "patient"
  colnames(clin)[colnames(clin)=="FIGO stage"] <- "simple_stage"
  clin$simple_stage[clin$simple_stage=="IIIC"] <- "III"
  clin$Grade <- paste0("G", clin$Grade)
  clin$Grade[clin$Grade=="GNA"] <- NA
  
  ## Molecular annotation:
  molec <- read_excel(path, sheet=2, trim_ws=TRUE, na=STRNAS)
  molec <- as.data.frame(molec)

  molec$patient <- NA
  molec$patient[grepl("^AOCS",molec$Sample_ID)] <- substr(molec$Sample_ID[grepl("^AOCS",molec$Sample_ID)], 1, 8)
  molec$patient[grepl("^M",molec$Sample_ID)] <- gsub("\\-.*", "", molec$Sample_ID[grepl("^M",molec$Sample_ID)])
  
  molec <- subset(molec, Sample_type=="Primary tumor")
  
  ## Outer Join:
  df <- merge(clin, molec, by="patient")
  # df$patient <- gsub("-", "_", df$patient)
  # df$Sample_ID <- gsub("-", "_", df$Sample_ID)
  return(df)
  
}

clin_aus <- load_aus_clinical()

all(clin_aus$patient %in% c(samps_aocs$case_id,samps_mocog$case_id)) #important checkpoint

# ---------------------------- Load & Clean Expression ----------------------------
eAOCS <- read.table(paste0(DIR, PATH_AUS, "GSE209964_analysis_14_08_2019.HTSeq.all.raw.txt"), header=TRUE)
eMOCOG <- read.table(paste0(DIR, PATH_AUS, "GSE211669_20_04_2020.HTSeq.all.raw.txt"), header=TRUE)
stopifnot( identical(eAOCS$GeneID, eMOCOG$GeneID) ) #important checkpoint

rownames(eAOCS) <- eAOCS$GeneID
rownames(eMOCOG) <- eMOCOG$GeneID
eAOCS$GeneID <- eMOCOG$GeneID <- NULL

## Subset each expression data to selected/filtered samples:
eAOCS <- eAOCS[ , colnames(eAOCS) %in% samps_aocs$Sample_title]
stopifnot( identical(colnames(eAOCS), samps_aocs$Sample_title) )

eMOCOG <- eMOCOG[ , colnames(eMOCOG) %in% samps_mocog$Sample_title] 
stopifnot(identical(colnames(eMOCOG), samps_mocog$Sample_title))

samps_aus <- rbind(samps_aocs, samps_mocog)
clin_aus <- merge(samps_aus, clin_aus, by.x="case_id", by.y="patient", all.x=TRUE)
cts <- cbind(eAOCS, eMOCOG)
stopifnot( identical(colnames(cts), clin_aus$Sample_title) ) #important checkpoint
# colnames(cts) <- clin_aus$Accession
rm(samps_aocs, samps_mocog, samps_aus, eAOCS, eMOCOG)

## Gene annotation:
cts$GeneID <- gsub("\\|.*", "", rownames(cts))
ENSEMBL <- read.csv(paste0(DIR, PATH_AUS,"ensembl_dict_GSE211669_column1.csv"))
stopifnot(identical(cts$GeneID, ENSEMBL$ensembl)) #checkpoint
cts$GeneID <- ENSEMBL$HUGO

cts <- subset(cts, GeneID %in% COMMON_GENES)

## Handle duplicate HUGO symbols:
sum(duplicated((cts$GeneID))) #367
cts <- aggregate(. ~ GeneID, data=cts, FUN=mean, na.rm=TRUE) 
rownames(cts) <- cts$GeneID
cts$GeneID <- NULL
cts <- data.matrix(cts)

# ---------------------------- Update Clinical Data w/ Gene Strata ----------------------------
cl_aus <- stratify_by_gene(cts, GENE, "Accession")
clin_aus <- merge(clin_aus, cl_aus, by="Accession")

cts <- cts[rownames(cts)!= GENE, ]

# ----------------------------- Export -----------------------------
## 1. Export standard-normalized: 
expr_aus <- apply(cts, 1, StandardScale) #normalize for each GENE! will transpose...
expr_aus <- t(expr_aus)
sum(is.na(expr_aus))
hist(expr_aus)

# save(
#   list = c("expr_aus","clin_aus"),
#   file = paste0(DIR_OUT, "240504b_aocs_mocog.RData"),
#   compress = TRUE
# )

## 2. Export raw counts:
# write.table(cts, paste0(DIR_OUT,"240505b_hugo_expr_counts_aocs.tsv"), sep="\t", row.names=TRUE, quote=FALSE)

## 3. Export log2 counts:
cts <- log2(cts + 1)
hist(cts)
# write.table(cts, paste0(DIR_OUT,"240505b_hugo_expr_log2count_aocs.tsv"), sep="\t", row.names=TRUE, quote=FALSE)
