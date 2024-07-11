# Visualize & Test ESTIMATE Output
# Notes:
# - Use wrappers to maximize consistency between cohorts & avoid code duplication

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggplot2)
library(reshape2)

eda_wrapper <- function(df, groupVar="Group0") {
  #'@description Wrapper for ggplot2 by mRNA group for each cohort
  df_plt <- df[ , c(groupVar, "StromalScore", "ImmuneScore", "ESTIMATEScore")]
  df_plt <- melt(df_plt, value.name="Score")
  colnames(df_plt)[1] <- "Group"
  df_plt$variable <- gsub("Score", "", df_plt$variable)
  
  ggplot(df_plt, aes(variable, Score, color=Group)) + 
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(position=position_jitterdodge(jitter.width=0.1)) +
    THEME_BOX
}

load(paste0(DIR_OUT,"240506_estimate.RData")) #ESTIMATE results saved
est_tcga <- apply(est_tcga, 2, StandardScale)
est_aus <- apply(est_aus, 2, StandardScale)

## TCGA:
load(paste0(DIR_OUT,"240504a_tcgaov.RData"))
clin_tcga <- merge(clin_tcga, est_tcga, by.x="bcr_sample_barcode", by.y="row.names")
eda_wrapper(clin_tcga)

t.test(ESTIMATEScore ~ Group0, data=clin_tcga)

fit_tcga <- glm(Group0 ~ ESTIMATEScore + AgeAtDx + simple_stage, data=clin_tcga, family=binomial)
summary(fit_tcga)

## AOCS:
load(paste0(DIR_OUT, "240504b_aocs_mocog.RData"))
clin_aus <- merge(clin_aus, est_aus, by.x="Sample_title", by.y="row.names")
eda_wrapper(clin_aus)

t.test(ESTIMATEScore ~ Group0, data=clin_aus)

fit_aus <- glm(Group0 ~ ESTIMATEScore + AgeAtDx + simple_stage, data=clin_aus, family=binomial)
summary(fit_aus)

## Combined plot for publication:
CSELE <- c("Group0", "ESTIMATEScore","ImmuneScore","StromalScore")

df_comb <- rbind(
  cbind(Cohort="American (TCGA)", clin_tcga[ , CSELE]),
  cbind(Cohort="Australian (AOCS+MOCOG)", clin_aus[ , CSELE])
)
df_comb <- reshape2::melt(df_comb)
df_comb$variable <- gsub("Score","",df_comb$variable)

ggplot(df_comb, aes(variable, value, color=Group0)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.1)) +
  geom_vline(xintercept=1.5, linetype=2) +
  facet_wrap(~Cohort, ncol=2) +
  THEME_BOX
