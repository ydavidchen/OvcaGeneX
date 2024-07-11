# Selected Visualizations
# Notes:

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggplot2)

GENES <- c("PARP1", "MKI67")

#------------------------------------ TCGA ------------------------------------ 
load(paste0(DIR_OUT,"240504a_tcgaov.RData"))

expr_tcga <- as.data.frame(t(expr_tcga[GENES, ]))

clin_tcga <- merge(clin_tcga, expr_tcga, by.x="bcr_sample_barcode", by.y="row.names")

ggplot(clin_tcga, aes(Group0, MKI67)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

t.test(MKI67 ~ Group0, data=clin_tcga)
fit_tcga <- glm(GroupM ~ MKI67 + AgeAtDx + simple_stage + Grade, data=clin_tcga, family=binomial)
summary(fit_tcga)


ggplot(clin_tcga, aes(Group0, PARP1)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

t.test(PARP1 ~ Group0, data=clin_tcga)

fit_tcga <- glm(GroupM ~ PARP1 + AgeAtDx + simple_stage + Grade, data=clin_tcga, family=binomial)
summary(fit_tcga)

# ------------------------------------ AOCS ------------------------------------ 
load(paste0(DIR_OUT, "240504b_aocs_mocog.RData"))

expr_aus <- as.data.frame(t(expr_aus[GENES, ]))

clin_aus <- merge(clin_aus, expr_aus, by.x="Sample_title", by.y="row.names")

ggplot(clin_aus, aes(Group0, MKI67)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

t.test(MKI67 ~ Group0, data=clin_aus)

fit_aus <- glm(Group0 ~ MKI67 + AgeAtDx + simple_stage + Grade, data=clin_aus, family=binomial)
summary(fit_aus)


ggplot(clin_aus, aes(Group0, PARP1)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  THEME_BOX

fit_aus <- glm(Group0 ~ PARP1 + AgeAtDx + simple_stage + Grade, data=clin_aus, family=binomial)
summary(fit_aus)
