# Name: compare logfold
# Author: EY
# Date: June 07 2024
# Version: Python 3.10
# Description: sort logfold to look for different patterns

import os
import pandas as pd

# Read the CSV files
directory = os.path.join(os.getcwd(), "deseq_results/x5h/")
files = [
    "result_AsterCLV3_5h.csv",
    "result_AsterCLV3_L5S_5h.csv",
    "result_AtCLV3p_5h.csv",
    "result_AtCLV3p_S5L_5h.csv"
]

dfs = {}
for file in files:
    file_path = os.path.join(directory, file)
    dfs[file] = pd.read_csv(file_path)


# Define a function to compare log2FoldChange values
def compare_changes(df1, df2, gene):
    gene_data1 = df1[df1.iloc[:, 0] == gene]
    gene_data2 = df2[df2.iloc[:, 0] == gene]
    if len(gene_data1) == 0 or len(gene_data2) == 0:
        return None
    pvalue1 = gene_data1.iloc[0]["padj"]
    pvalue2 = gene_data2.iloc[0]["padj"]
    if pvalue1 < 0.05 and pvalue2 < 0.05:
        fold_change1 = gene_data1.iloc[0]["log2FoldChange"]
        fold_change2 = gene_data2.iloc[0]["log2FoldChange"]
        return fold_change1, fold_change2
    else:
        return None


counter = 0
# Compare log2FoldChange values for each gene
gene_list = []
for gene in dfs[files[2]].iloc[:, 0]:
    changes = compare_changes(dfs[files[2]], dfs[files[1]], gene)
    if changes is not None:
        print(gene)
        counter += 1
        fold_change_AtCLV3p_5h, fold_change_AsterCLV3_L5S_5h = changes
        fold_change_AsterCLV3_5h = dfs[files[0]][dfs[files[0]].iloc[:, 0] == gene].iloc[0]["log2FoldChange"]
        fold_change_AtCLV3p_S5L_5h = dfs[files[3]][dfs[files[3]].iloc[:, 0] == gene].iloc[0]["log2FoldChange"]
        if (
                fold_change_AtCLV3p_5h > 0 and fold_change_AsterCLV3_L5S_5h > 0 and fold_change_AsterCLV3_5h < 0 and fold_change_AtCLV3p_S5L_5h < 0) or \
                (
                        fold_change_AtCLV3p_5h < 0 and fold_change_AsterCLV3_L5S_5h < 0 and fold_change_AsterCLV3_5h > 0 and fold_change_AtCLV3p_S5L_5h > 0):
            gene_list.append(gene)

# Print the list of genes
print("Genes where the change in result_AtCLV3p_5h.csv AND result_AsterCLV3_L5S_5h.csv is opposite to the others:")
print(gene_list)

# compare if log fold is 2x greater or less since none of them switch direction
gene_list_2x = []
for gene in dfs[files[2]].iloc[:, 0]:
    changes = compare_changes(dfs[files[2]], dfs[files[1]], gene)
    if changes is not None:
        fold_change_AtCLV3p_5h, fold_change_AsterCLV3_L5S_5h = changes
        fold_change_AsterCLV3_5h = dfs[files[0]][dfs[files[0]].iloc[:, 0] == gene].iloc[0]["log2FoldChange"]
        fold_change_AtCLV3p_S5L_5h = dfs[files[3]][dfs[files[3]].iloc[:, 0] == gene].iloc[0]["log2FoldChange"]
        if (abs(fold_change_AtCLV3p_5h) > abs(fold_change_AsterCLV3_5h) * 2 and abs(fold_change_AtCLV3p_5h) > abs(
                fold_change_AtCLV3p_S5L_5h) * 2 and abs(fold_change_AsterCLV3_L5S_5h) > abs(
                fold_change_AsterCLV3_5h) * 2 and abs(fold_change_AsterCLV3_L5S_5h) > abs(
                fold_change_AtCLV3p_S5L_5h) * 2) or \
                (abs(fold_change_AtCLV3p_5h) < abs(fold_change_AsterCLV3_5h) / 2 and abs(fold_change_AtCLV3p_5h) < abs(
                    fold_change_AtCLV3p_S5L_5h) / 2 and abs(fold_change_AsterCLV3_L5S_5h) < abs(
                    fold_change_AsterCLV3_5h) / 2 and abs(fold_change_AsterCLV3_L5S_5h) < abs(
                    fold_change_AtCLV3p_S5L_5h) / 2):
            print("fold_change_AtCLV3p_5h:", fold_change_AtCLV3p_5h)
            print("fold_change_AsterCLV3_L5S_5h", fold_change_AsterCLV3_L5S_5h)
            print("vs")
            print("fold_change_AsterCLV3_5h", fold_change_AsterCLV3_5h)
            print("fold_change_AtCLV3p_S5L_5h", fold_change_AtCLV3p_S5L_5h)
            gene_list_2x.append(gene)
print(gene_list_2x)
# Print the list of genes
print(
    "Genes where fold_change_AsterCLV3_5h and fold_change_AsterCLV3_L5S_5h are 50% larger or smaller than the other two:")
print(gene_list)
