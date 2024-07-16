# Selected Gene-Gene Visualizations & Tests
# Notes:

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggplot2)

GENES <- c("PARP1")

#------------------------------------ TCGA ------------------------------------ 
load(paste0(DIR_OUT,"240504a_tcgaov.RData"))

expr_tcga <- as.data.frame(t(expr_tcga[GENES, , drop=FALSE]))

clin_tcga <- merge(clin_tcga, expr_tcga, by.x="bcr_sample_barcode", by.y="row.names")

## Continuous:
ggplot(clin_tcga, aes(mRNA_ss, PARP1)) +
  geom_point() +
  geom_smooth(method="lm") +
  THEME_SCATTER

cor.test(clin_tcga$mRNA_ss, clin_tcga$PARP1)

## Binary:
ggplot(clin_tcga, aes(Group, PARP1)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

t.test(PARP1 ~ Group, data=clin_tcga)

fit_tcga <- glm(Group ~ PARP1 + AgeAtDx + simple_stage, data=clin_tcga, family=binomial)
summary(fit_tcga)

# ------------------------------------ AOCS ------------------------------------ 
load(paste0(DIR_OUT,"240504b_aocs_mocog.RData"))

expr_aus <- as.data.frame(t(expr_aus[GENES, , drop=FALSE]))

clin_aus <- merge(clin_aus, expr_aus, by.x="Sample_title", by.y="row.names")

## Continuous:
ggplot(clin_aus, aes(mRNA_ss, PARP1)) +
  geom_point() +
  geom_smooth(method="lm") +
  THEME_SCATTER

cor.test(clin_aus$mRNA_ss, clin_aus$PARP1)

## Binary:
ggplot(clin_aus, aes(Group, PARP1)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

t.test(PARP1 ~Group, data=clin_aus)

fit_aus <- glm(Group ~ PARP1 + AgeAtDx + simple_stage, data=clin_aus, family=binomial)
summary(fit_aus)
