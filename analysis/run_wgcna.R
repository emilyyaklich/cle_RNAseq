# run wGCNA


setwd('/home/ely67071/cle_RNAseq/')

# read in deseq results (output from run_DGE_deseq_sunflower_inflo.R)
# note this is reading in all of the data and is NOT filtered for genes that are
# differentially expressed 
deseq <- readRDS('deseq_results/deseq_dataset_results_combatseq.RData')
deseq$samples



# pre-filter for reads where at least 10 samples have a count of 10 or higher
keep<-rowSums(counts(deseq)>=10)>=10
length(which(keep==1))
deseq_filt<-deseq[keep,]

# transform data into a matrix
vsd <- vst(deseq,blind=TRUE)
vsd_matrix<-assay(vsd)

input_mat=t(vsd_matrix)
