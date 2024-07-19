# Prepare Expression Matrices & Stratify Cohort
# Script Objectives:
# 1. Clean, assemble, restrict/filter, & load clinical + sample annotations.
# 2. Clean & filter expression & identify samples from RNAseq.
# 3. Export cleaned expression & clinical data for downstream use w/ convenience & consistency.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(readxl)

stratify_by_gene <- function(expr, geneName, thresh=0, idcol="sample") {
  #'@description Wrapper to stratify population by a gene at a threshold
  #'@param idCol Sample name to use in resultant dataframe
  df <- as.data.frame(expr[geneName, ])
  colnames(df) <- "mRNA"
  df$Group <- factor(ifelse(df$mRNA > thresh, "High", "Low"), levels=c("High","Low"))
  df[idcol] <- rownames(df)
  rownames(df) <- NULL
  return(df[ , c(idcol,"mRNA","Group")])
}

load_meta_tcga <- function(prim_key="bcr_sample_barcode", sec_key="bcr_patient_barcode") {
  #'@description Wrapper to assemble TCGA-OV sample metadata in isolated environment
  ## TCGA-CDR (Liu 2018 Cell) w/ some clinical & gold-standard survival profiles:
  liu2018 <- read_excel(paste0(DIR,"TCGA_MISC/Liu2018CellTCGACDR.xlsx"), sheet=1, na=c(STRNAS,"[Not Available]","[Not Applicable]","[Unknown]","[Discrepancy]"))
  liu2018 <- subset(liu2018, type=="OV")
  colnames(liu2018)[colnames(liu2018)=="age_at_initial_pathologic_diagnosis"] <- "AgeAtDx"
  colnames(liu2018)[colnames(liu2018)=="histological_grade"] <- "Grade"
  liu2018$raceWhite <- liu2018$race == "WHITE"
  liu2018 <- liu2018[ , c(sec_key,"AgeAtDx","Grade","raceWhite","tumor_status","OS.time","OS","PFI","PFI.time")]
  
  ## Manually cleaned GDAC FireBrowse w/ additional clinical data:
  gdac <- read.csv(paste0(DIR, "SynapseTCGAlive/OV/gdac_firebrowse_ov_cleaned.csv"), na.strings=STRNAS)
  gdac$bcr_patient_barcode <- toupper(gdac$bcr_barcode)
  gdac$Neoadjuvant <- gdac$history_of_neoadjuvant_treatment=="yes"
  gdac$Radiation <- gdac$radiation_therapy=="yes"
  gdac$clinical_stage <- toupper(gsub("stage ", "", gdac$clinical_stage))
  gdac$simple_stage <- gsub("[A-C]", "", gdac$clinical_stage, ignore.case=FALSE)
  gdac <- gdac[ , c(sec_key,"ethnicity","clinical_stage","simple_stage","Neoadjuvant","Radiation")]
  
  ## ASCAT ploidy & purity from Steele 2022 Nature:
  steele2022 <- read_excel(paste0(DIR,"TCGA_MISC/Steele2022_Nature_cnsig/TableS1_ascat_cn_for_all_platforms.xlsx"), sheet=3)
  steele2022 <- as.data.frame(steele2022)
  steele2022[prim_key] <- substr(steele2022$barcodeTumour, 1, 16)
  steele2022 <- steele2022[ , c(prim_key,"ploidy","purity")]
  
  ## Synapse bio-samples:
  samples <- read.csv(paste0(DIR, "SynapseTCGAlive/OV/nationwidechildrens.org_ov_bio.sample.tsv"), sep="\t", na.strings=STRNAS)
  samples <- subset(samples, sample_type=="Primary Tumor") #exclude blood normals etc.
  samples <- samples[ , c(prim_key,"sample_type")]
  samples[sec_key] <- substr(samples$bcr_sample_barcode, 1, 12) #no bracket use
  
  ## Outer joins:
  res <- merge(liu2018, gdac, by=sec_key, all=TRUE) #patient-level
  res <- merge(res, samples, by=sec_key, all=TRUE) #add in 1st sample-level mapper
  res <- merge(res, steele2022, by=prim_key, all.x=TRUE) #sample-level
  
  ## Restrict:
  res <- subset(res, simple_stage %in% c("III","IV") & Grade %in% c("G2","G3"))
  return(res)
}

load_meta_aocs <- function(path=paste0(DIR,PATH_AUS,"clinical_sele.xlsx")) {
  #'@description Wrapper to assemble published AOCS clinical & molecular profiles in isolated environment
  ## Clinical data (manually curated from Garsed et al. Suppl Tables):
  clin <- read_excel(path, sheet=1, trim_ws=TRUE, na=c(STRNAS,"Nil","Not known","Size unknown"))
  clin <- as.data.frame(clin)
  colnames(clin)[1] <- "patient"
  
  clin$simple_stage <- gsub("[A-C]", "", clin$`FIGO stage`)
  clin$`Residual disease` <- clin$`Residual disease` == "> 1 cm"
  
  clin <- subset(clin, ! is.na(Grade))
  clin$Grade <- paste0("G", clin$Grade)
  
  clin <- subset(clin, Neoadjuvant=="No")
  
  ## Molecular annotation:
  molec <- read_excel(path, sheet=2, trim_ws=TRUE, na=STRNAS)
  molec <- as.data.frame(molec)
  
  molec$patient <- NA
  molec$patient[grepl("^AOCS",molec$Sample_ID)] <- substr(molec$Sample_ID[grepl("^AOCS",molec$Sample_ID)], 1, 8)
  molec$patient[grepl("^M",molec$Sample_ID)] <- gsub("\\-.*", "", molec$Sample_ID[grepl("^M",molec$Sample_ID)])
  
  molec <- subset(molec, Sample_type=="Primary tumor")
  
  ## Sample map:
  sampmap <- load_curated_aocs_sampmap() #sample sheet

  ## Inner-join:
  df <- merge(clin, molec, by="patient")
  df <- merge(sampmap, df, by.x="case_id", by.y="patient")
  return(df)
}

load_curated_aocs_sampmap <- function(dir=paste0(DIR,PATH_AUS)) {
  #'@description Wrapper to assemble GEO sample sheets 
  #'@note Required for expression sample-name matching
  #'@note Sub-cohort labeled for name matching & ComBat
  
  ## ICGC sub-cohort (Patch 2015 Nature):
  samps_icgc <- read.csv(paste0(dir,"GSE209964_series_matrix_curated.csv"), strip.white=TRUE)
  samps_icgc <- subset(samps_icgc, SampleType=="PrimaryTumour")
  samps_icgc <- subset(samps_icgc, CancerType=="High Grade Serous Ovarian Cancer")
  
  ## MOCOG sub-cohort (Garsed 2022 Nature Genet):
  samps_mocog <- read.csv(paste0(dir,"GSE211669_twoGPLs_curated.csv"), strip.white=TRUE) #all hgsocs
  samps_mocog <- subset(samps_mocog, SampleType=="Primary Tumor")
  samps_mocog$Reanalysis <- gsub("Re-analysis of ", "", samps_mocog$Reanalysis)
  samps_mocog <- subset(samps_mocog, ! Reanalysis %in% samps_icgc$Accession)
  
  COL_SAMPS <- c("Accession","Sample_title","case_id","SampleType")
  samps_icgc <- samps_icgc[ , COL_SAMPS]
  samps_mocog <- samps_mocog[ , COL_SAMPS]
  return(rbind(
    cbind(samps_icgc, Subcohort="ICGC"), 
    cbind(samps_mocog, Subcohort="MOCOG")
  ))
}

# ---------------------------- TCGA-OV ---------------------------- 
## Load & pre-process expression matrix:
expr_tcga <- read.csv(paste0(DIR,"SynapseTCGAlive/OV/unc.edu_OV_IlluminaHiSeq_RNASeqV2.geneExp.tsv"), sep="\t", check.names=FALSE)
expr_tcga$gene_id <- gsub("\\|.*", "", expr_tcga$gene_id)
expr_tcga <- subset(expr_tcga, ! (duplicated(gene_id) | duplicated(gene_id, fromLast=TRUE)))
rownames(expr_tcga) <- expr_tcga$gene_id
expr_tcga$gene_id <- NULL
expr_tcga <- data.matrix(expr_tcga)

nrow(expr_tcga) #20500
expr_tcga <- drop_flat_vals(expr_tcga)
nrow(expr_tcga) #20161

## Handle samples & Mutually subset
clin_tcga <- load_meta_tcga()
anyDuplicated(substr(colnames(expr_tcga), 1, 16)) #checkpoint: 0
colnames(expr_tcga) <- substr(colnames(expr_tcga), 1, 16)

expr_tcga <- expr_tcga[ , colnames(expr_tcga) %in% clin_tcga$bcr_sample_barcode]
clin_tcga <- subset(clin_tcga, bcr_sample_barcode %in% colnames(expr_tcga))

## Scale data:
expr_tcga <- t(apply(expr_tcga, 1, StandardScale))
expr_tcga <- winsorize(expr_tcga, -10, 10)
hist(expr_tcga)

## Define groups & update annotation:
groups_tcga <- stratify_by_gene(expr_tcga, GENE, 0, "bcr_sample_barcode")
clin_tcga <- merge(clin_tcga, groups_tcga, by="bcr_sample_barcode")

# ---------------------------- AOCS ----------------------------
## Sample sheet labeled with sub-cohorts:
clin_aocs <- load_meta_aocs()

## Combine expression matrices of 2 sub-cohorts:
eICGC <- read.table(paste0(DIR, PATH_AUS, "GSE209964_analysis_14_08_2019.HTSeq.all.raw.txt"), header=TRUE)
rownames(eICGC) <- eICGC$GeneID
eICGC$GeneID <- NULL
eICGC <- data.matrix(eICGC)

eICGC <- eICGC[ , colnames(eICGC) %in% clin_aocs$Sample_title[clin_aocs$Subcohort=="ICGC"]]
dim(eICGC) #57773    69

eICGC <- drop_uniform_var(eICGC) #order matters
dim(eICGC) #55111    69

eMOCOG <- read.table(paste0(DIR, PATH_AUS, "GSE211669_20_04_2020.HTSeq.all.raw.txt"), header=TRUE)
rownames(eMOCOG) <- eMOCOG$GeneID
eMOCOG$GeneID <- NULL
eMOCOG <- data.matrix(eMOCOG)

eMOCOG <- eMOCOG[ , colnames(eMOCOG) %in% clin_aocs$Sample_title[clin_aocs$Subcohort=="MOCOG"]] 
dim(eMOCOG) #57773    53

eMOCOG <- drop_uniform_var(eMOCOG) #order matters
dim(eMOCOG) #51173    53

## Combine cohorts & handle gene names:
ENSEMBL <- read.csv(paste0(DIR,PATH_AUS,"ensembl_dict_GSE211669_column1.csv"))
ENSEMBL <- subset(ENSEMBL, chr %in% CHROMS)
ENSEMBL <- subset(ENSEMBL, category == "protein_coding")
ENSEMBL <- subset(ENSEMBL, ! (duplicated(HUGO) | duplicated(HUGO, fromLast=TRUE)))

expr_aocs <- merge(eICGC, eMOCOG, by="row.names")
colnames(expr_aocs)[1] <- "ensembl"
expr_aocs$ensembl <- gsub("\\|.*", "", expr_aocs$ensembl)
expr_aocs <- merge(expr_aocs, ENSEMBL[ , c("ensembl","HUGO")], by="ensembl")
rownames(expr_aocs) <- expr_aocs$HUGO #also ensures no duplicate
expr_aocs$HUGO <- expr_aocs$ensembl <- NULL
expr_aocs <- data.matrix(expr_aocs)
nrow(expr_aocs) #18864
stopifnot(identical(colnames(expr_aocs), clin_aocs$Sample_title)) #important checkpoint
rm(eICGC, eMOCOG, ENSEMBL)

## ComBat, Scale, & Winsorize:
mod_combat <- model.matrix( ~ 1, data=clin_aocs)
expr_aocs <- sva::ComBat(dat=expr_aocs, batch=clin_aocs$Subcohort, mod=mod_combat, par.prior=TRUE)

expr_aocs <- t(apply(expr_aocs, 1, StandardScale))
expr_aocs <- winsorize(expr_aocs, -10, 10)
hist(expr_aocs)

## Define gene strata & update clinical:
groups_aocs <- stratify_by_gene(expr_aocs, GENE, 0, "Sample_title")
clin_aocs <- merge(clin_aocs, groups_aocs, by="Sample_title")

# ---------------------- Find Common, Subset & Export ---------------------- 
## Find common:
gShared <- intersect(rownames(expr_tcga), rownames(expr_aocs))
length(gShared) #16512

expr_tcga <- expr_tcga[gShared, ]
expr_aocs <- expr_aocs[gShared, ]

## Important checkpoints:
stopifnot(identical(colnames(expr_tcga), clin_tcga$bcr_sample_barcode)) 
stopifnot(identical(colnames(expr_aocs), clin_aocs$Sample_title))

## Expression - RData for efficient storage & loading:
save(list=c("expr_tcga","expr_aocs"), file=paste0(DIR_OUT,"240504_processed_expression.RData"), compress=TRUE)

## Clinical - RData to preserve R datatypes & efficient loading:
save(list=c("clin_tcga", "clin_aocs"), file=paste0(DIR_OUT,"240504_curated_cohort_meta.RData"), compress=TRUE)
