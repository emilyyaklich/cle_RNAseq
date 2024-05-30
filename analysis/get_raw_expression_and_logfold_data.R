# Name: get_raw_expression_and_logfold_data.R
# Author: EY
# Date: 05/30/2024
# Version:4.2.1
# Description: will plot raw expression from combatseq and the logfold change from DESeq

setwd('/home/ely67071/cle_RNAseq/analysis')


library(tidyverse)
library(stringr)


# read in adjusted count data 
adjusted_counts<-read.csv('adjusted_counts_combatseq.csv', row.names=1)


## For any list of gene IDs, extract the raw count data and plot it ##

# Define the list of gene IDs
gene_ids <- c("g51546.t1", "g66.t1", "g53294.t1", "g44552.t1")

# Initialize an empty list to store the results for each gene
all_results <- list()

# Loop over each gene ID
for (gene_id in gene_ids) {
  # Extract subset of counts for the current gene
  subset_gene <- adjusted_counts[rownames(adjusted_counts) == gene_id, ]
  
  # Reshape from wide to long format, extracting the numeric part
  df_long <- subset_gene %>%
    pivot_longer(cols = everything(), names_to = "column_name") %>%
    mutate(col_nums = as.numeric(str_extract(column_name, "(?<=\\.)(\\d+)(?=_.*$)")))
  
  # Calculate the mean value for each group
  result <- df_long %>%
    group_by(col_nums) %>%
    summarise(mean_value = mean(value))
  
  # Add dataset_id information
  result$dataset_id <- c("5h_control", "5h_AtCLV3", "5h_AtVar", "5h_AsterCLV3", "5h_AsterCLV3Var", "24h_control", "24h_AtCLV3", "24h_AtVar", "24h_AsterCLV3", "24h_AsterCLV3Var")
  result$dataset_id <- factor(result$dataset_id, levels = c("5h_control", "5h_AtCLV3", "5h_AtVar", "5h_AsterCLV3", "5h_AsterCLV3Var", "24h_control", "24h_AtCLV3", "24h_AtVar", "24h_AsterCLV3", "24h_AsterCLV3Var"))
  
  # Save the result to the list
  all_results[[gene_id]] <- result
  
  # Write the result to a CSV file
  write.csv(result, file = paste0(gene_id, '.csv'), row.names = FALSE)
  
  # Plot the results and save as PNG
  png(paste0("plots/", gene_id, "_exp.png"), res = 215, width = 2000, height = 1500)
  ggplot(result, aes(x = as.factor(dataset_id), y = mean_value)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = gene_id, x = "dataset_id", y = "Expression Counts (Raw)")
  dev.off()
}


## For any given gene in a list get the log fold change data and plot it ##

# read in the data for 5h
DEData_pairwise_cs_5h<-ImportCSVs('deseq_results/x5h/',0.05)

# read in the data for 24h
DEData_pairwise_cs_24h<-ImportCSVs('deseq_results/x24h/',0.05)

# Define the list of gene IDs
gene_ids <- c("g51546.t1")

# Loop over each gene ID
for (gene_id in gene_ids) {
  # Bind rows from different datasets for the current gene
  logFC_data <- bind_rows(
    DEData_pairwise_cs_5h$result_AtCLV3p_5h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AtCLV3p_5h"),
    DEData_pairwise_cs_5h$result_AtCLV3p_S5L_5h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AtCLV3p_S5L_5h"),
    DEData_pairwise_cs_5h$result_AsterCLV3_5h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AsterCLV3_5h"),
    DEData_pairwise_cs_5h$result_AsterCLV3_L5S_5h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AsterCLV3_L5S_5h"),
    DEData_pairwise_cs_24h$result_AtCLV3p_24h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AtCLV3p_24h"),
    DEData_pairwise_cs_24h$result_AtCLV3p_S5L_24h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AtCLV3p_S5L_24h"),
    DEData_pairwise_cs_24h$result_AsterCLV3_24h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AsterCLV3_24h"),
    DEData_pairwise_cs_24h$result_AsterCLV3_L5S_24h %>%
      filter_all(any_vars(. %in% gene_id)) %>%
      select(log2FoldChange, padj) %>%
      mutate(ID = "AsterCLV3_L5S_24h")
  )
  
  # Set factor levels (for axis ordering)
  logFC_data$ID <- factor(logFC_data$ID, levels = c("AtCLV3p_5h", "AtCLV3p_24h", "AtCLV3p_S5L_5h", "AtCLV3p_S5L_24h", "AsterCLV3_5h", "AsterCLV3_24h", "AsterCLV3_L5S_5h", "AsterCLV3_L5S_24h"))
  
  # Plot the results and save as PNG
  png(paste0("plots/", gene_id, "_logfold.png"), res = 215, width = 2000, height = 1500)
  ggplot(logFC_data, aes(x = as.factor(ID), y = log2FoldChange, fill = padj < 0.05)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("grey", "blue"), labels = c("padj > 0.05", "padj < 0.05")) +
    labs(title = gene_id, x = "dataset_id", y = "log2FoldChange")
  dev.off()
}

