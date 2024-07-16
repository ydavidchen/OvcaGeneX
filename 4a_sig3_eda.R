# EDA of SBS3
# Notes:
# - TCGA: Computed from Alexandrov 2020 trinucleotide tables & COSMIC look-up
# - AOCS: Available directly from Garsed 2022 Suppl Tab 9

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggplot2)

# ---------------------------- TCGA ----------------------------
load(paste0(DIR_OUT,"240504a_tcgaov.RData")); rm(expr_tcga)

sbs3_tcga <- read.csv(paste0(DIR_OUT,"20240714_tcga_sbs3.csv"), row.names=1)
sbs3_tcga$bcr_sample_barcode <- substr(rownames(sbs3_tcga), 1, 16)
sbs3_tcga$Signature3 <- MinMaxScale(sbs3_tcga$SBS3_contrib)

sbs3_tcga <- merge(sbs3_tcga, clin_tcga, by="bcr_sample_barcode")

ggplot(sbs3_tcga, aes(Group, Signature3)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

wilcox.test(Signature3 ~ Group, data=sbs3_tcga)

fit_sig3_tcga <- glm((Group=="High") ~ Signature3 + AgeAtDx + simple_stage, data=sbs3_tcga, family=binomial)
summary(fit_sig3_tcga)

# ---------------------------- AOCS ----------------------------
load(paste0(DIR_OUT,"240504b_aocs_mocog.RData")); rm(expr_aus)

sbs3_aocs <- readxl::read_excel(paste0(DIR,"ovarian_datasets/MOCOG/SupplementaryTables1-17.xlsx"), sheet=10, skip=2)
sbs3_aocs <- as.data.frame(sbs3_aocs)
rownames(sbs3_aocs) <- sbs3_aocs$SampleID
rownames(sbs3_aocs) <- gsub("AOCS-", "AOCS", rownames(sbs3_aocs)) #temp
rownames(sbs3_aocs) <- gsub("\\-.*", "", rownames(sbs3_aocs))
rownames(sbs3_aocs) <- gsub("AOCS", "AOCS-", rownames(sbs3_aocs)) #restore

sbs3_aocs <- sbs3_aocs[ , "sbs.SBS3", drop=FALSE]
sbs3_aocs$Signature3 <- MinMaxScale(sbs3_aocs$sbs.SBS3)

all(rownames(sbs3_aocs) %in% clin_aus$case_id) #checkpoint

sbs3_aocs <- merge(clin_aus, sbs3_aocs, by.x="case_id", by.y="row.names")
hist(sbs3_aocs$Signature3)

ggplot(sbs3_aocs, aes(Group, Signature3)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

wilcox.test(Signature3 ~ Group, data=sbs3_aocs)

fit_sig3_aocs <- glm((Group=="High") ~ Signature3 + AgeAtDx + simple_stage, data=sbs3_aocs, family=binomial)
summary(fit_sig3_aocs)
