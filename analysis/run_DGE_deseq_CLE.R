# Name: run_DGE_sunflower_inflo_combatseq.R
# Author: EY
# Date: 08/29/2023
# Version:4.2.1
# Description: Will run the differential gene expression on the infloresence stages accounting for variation with combatseq


setwd('/home/ely67071/cle_RNAseq/analysis')


library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
library(stringr)
library(tidyr)
source("Functions.R")


# read in and process data
adjusted_counts<-read.csv('adjusted_counts_combatseq.csv', row.names=1)
metadata<-read.csv('metadata.csv', row.names=1)
# create the factors of interest
metadata$dev_stage<-factor(metadata$cle_id)


# create the model with adjusted counts
adjusted_counts_deseq<-DESeqDataSetFromMatrix((adjusted_counts),colData = metadata, design=~0+cle_id)

# run the DGE analysis
DESeq_dataset_results_combatseq<-DESeq(adjusted_counts_deseq,parallel=TRUE)

# save this because it will be read in to WGCNA (WGCNA takes all genes into account, but uses the normalization of DESeq)
saveRDS(DESeq_dataset_results_combatseq, file='deseq_results/deseq_dataset_results_combatseq.RData')


resultsNames(DESeq_dataset_results_combatseq)
?results

# set up the pairwise contrasts 
#AtCLV3p vs control (cle1) at 5 hr
result_AtCLV3p_5h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE2","CLE1"),alpha=0.05,parallel=TRUE)

#AtCLV3 S5L
result_AtCLV3p_S5L_5h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE3","CLE1"),alpha=0.05,parallel=TRUE)

#Aster CLV3
result_AsterCLV3_5h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE4","CLE1"),alpha=0.05,parallel=TRUE)

#AsterCLV3 L5S vs control (cle1) at 5 hr
result_AsterCLV3_L5S_5h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE5","CLE1"),alpha=0.05,parallel=TRUE)




# same as above but at 24 hr
#AtCLV3p vs control (cle1) at 24 hr
result_AtCLV3p_24h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE7","CLE6"),alpha=0.05,parallel=TRUE)

#AtCLV3 S5L
result_AtCLV3p_S5L_24h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE8","CLE6"),alpha=0.05,parallel=TRUE)

#Aster CLV3
result_AsterCLV3_24h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE9","CLE6"),alpha=0.05,parallel=TRUE)


#AsterCLV3 L5S vs control (cle1) at 24 hr
result_AsterCLV3_L5S_24h<-results(DESeq_dataset_results_combatseq,contrast=c("cle_id","CLE10","CLE6"),alpha=0.05,parallel=TRUE)


# write to CSV file
write.csv(as.data.frame(result_AtCLV3p_5h), file='deseq_results/x5h/result_AtCLV3p_5h.csv')
write.csv(as.data.frame(result_AtCLV3p_S5L_5h), file='deseq_results/x5h/result_AtCLV3p_S5L_5h.csv')
write.csv(as.data.frame(result_AsterCLV3_5h), file='deseq_results/x5h/result_AsterCLV3_5h.csv')
write.csv(as.data.frame(result_AsterCLV3_L5S_5h), file='deseq_results/x5h/result_AsterCLV3_L5S_5h.csv')
write.csv(as.data.frame(result_AtCLV3p_24h), file='deseq_results/x24h/result_AtCLV3p_24h.csv')
write.csv(as.data.frame(result_AtCLV3p_S5L_24h), file='deseq_results/x24h/result_AtCLV3p_S5L_24h.csv')
write.csv(as.data.frame(result_AsterCLV3_24h), file='deseq_results/x24h/result_AsterCLV3_24h.csv')
write.csv(as.data.frame(result_AsterCLV3_L5S_24h), file='deseq_results/x24h/result_AsterCLV3_L5S_24h.csv')


