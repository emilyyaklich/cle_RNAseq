# Name: analyze DGE up v down
# Author: EY (based off of code written by E. Dittmar)
# Date: 01/24/2023 Updated: 05/30/2024
# Version:4.1.2
# Description: Will analyze the output from the DESeq DGE cult in terms of
# overexpression and underexpression of genes 
# need Functions.R written by ED


setwd('/home/ely67071/cle_RNAseq/analysis')
#install.packages("devtools")
# install package used to plot upset
#devtools::install_github("krassowski/complex-upset", force=TRUE)
#install.packages("lubridate")
library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
library(ComplexUpset)
library(purrr)
library(naniar)
library(stringr)
source("Functions.R")



# (1) LOOK AT DATA BETWEEN @ 5H
# read in the data
DEData_5h<-ImportCSVs('deseq_results/x5h/',0.05)
# filter out significant results
mydataSig_5h<-lapply(DEData_5h,SigDEdf,PvaluesCol=7,CritP=0.05)

# different list of df for up or down expression
dataSigup_5h<-lapply(mydataSig_5h, MoreCritNum,column=3, critNum=0)
dataSigdown_5h<-lapply(mydataSig_5h, LessCritNum,column=3, critNum=0)


# (2A) Take steps to plot up and down on the same plot with "difference bar"

# add a column labeling "up" or "down" to each list of dataframes 
dataSigup_5h<-lapply(dataSigup_5h, transform, Direction="Up")
dataSigdown_5h<-lapply(dataSigdown_5h, transform, Direction="Down")

# combine the two lists of data frames 
mydataSig_5h<-mapply(rbind, dataSigup_5h, dataSigdown_5h,SIMPLIFY=FALSE)

mydataSig_5h$result_AtCLV3p_S5L_5h
# rename the dataframes in the list
names(mydataSig_5h) <-c("AsterCLV3_5h", "AsterCLV3_L5S_5h","AtCLV3p_5h", "AtCLV3p_S5L_5h")

# subset the dataframes to only contain the gene ID and the direction of expression
gene_direction_only<-lapply(mydataSig_5h, function(x) {x[c("Gene", "Direction")]})
gene_direction_only<-gene_direction_only %>% purrr::reduce(full_join,by=c('Gene'))
# combine all three columns of direction into one string with where direction 
# is outlined. For example. "Up Up Down" If the gene is not in a treatment NA is used
gene_direction_only$dif<-paste(gene_direction_only$Direction.x,gene_direction_only$Direction.y,gene_direction_only$Direction.x.x,gene_direction_only$Direction.y.y, gene_direction_only$Direction)

# loop through all of the genes present and add the genes to the "up", "down",
# or "difference" lists depending on the directionality of expression
difference=c()
up=c()
down=c()
for(i in 1:nrow(gene_direction_only)){
  if(grepl('Down',gene_direction_only[i,"dif"])&grepl("Up",gene_direction_only[i,"dif"])){
    difference<-append(difference,gene_direction_only[i,"Gene"])
  } else if (grepl('Up',gene_direction_only[i,"dif"])) {
    up<-append(up,gene_direction_only[i,"Gene"])
  } else {
    down<-append(down,gene_direction_only[i,"Gene"])
}
}
#gene_direction_only<-names(gene_direction_only)
# for each list created in the above loop, turn it into a df and then rbind them into one df
difference<-data.frame(difference)
difference<-cbind(difference,x="Difference")
colnames(difference)[1]<-"Gene"
up<-data.frame(up)
up<-cbind(up,x="Up")
colnames(up)[1]<-"Gene"
down<-data.frame(down)
down<-cbind(down,x="Down")
colnames(down)[1]<-"Gene"
# df with Gene ID as column 1 and direction of expression as column 2 
up_or_down<-rbind(down,up,difference)

# in our initial dataframe containing all treatments, add the directionality of the gene
mydataSig_5h$AsterCLV3_5h$Direction<-up_or_down$x[match(mydataSig_5h$AsterCLV3_5h$Gene,up_or_down$Gene)]
mydataSig_5h$AsterCLV3_L5S_5h$Direction<-up_or_down$x[match(mydataSig_5h$AsterCLV3_L5S_5h$Gene,up_or_down$Gene)]
mydataSig_5h$AtCLV3p_5h$Direction<-up_or_down$x[match(mydataSig_5h$AtCLV3p_5h$Gene,up_or_down$Gene)]
mydataSig_5h$AtCLV3p_S5L_5h$Direction<-up_or_down$x[match(mydataSig_5h$AtCLV3p_S5L_5h$Gene,up_or_down$Gene)]

# create a column in each df that contains a treatment variable
mydataSig_5h$AsterCLV3_5h<-cbind(mydataSig_5h$AsterCLV3_5h, Treatment="AsterCLV3_5h")
mydataSig_5h$AsterCLV3_L5S_5h<-cbind(mydataSig_5h$AsterCLV3_L5S_5h, Treatment="AsterCLV3_L5S_5h")
mydataSig_5h$AtCLV3p_5h<-cbind(mydataSig_5h$AtCLV3p_5h, Treatment="AtCLV3p_5h")
mydataSig_5h$AtCLV3p_S5L_5h<-cbind(mydataSig_5h$AtCLV3p_S5L_5h, Treatment="AtCLV3p_S5L_5h")


#mydataSig_treatment_full<-mydataSig_treatment %>% reduce(full_join,by=c('Treatment'))
# join the df by gene ID
#mydataSig_treatment_full<-mydataSig_treatment %>% purrr::reduce(full_join,by=c('Gene'))

# combine all the dataframes on the rows
mydataSig_5h_full<-bind_rows(mydataSig_5h$AsterCLV3_5h, mydataSig_5h$AsterCLV3_L5S_5h,mydataSig_5h$AtCLV3p_5h,mydataSig_5h$AtCLV3p_S5L_5h)
mydataSig_5h_full[,"Gene"]
comparisons<-c("AsterCLV3_5h", "AsterCLV3_L5S_5h", "AtCLV3p_5h","AtCLV3p_S5L_5h")
treatment<-c("Treatment")

# one-hot encode the treatment variables -- why did I do this?
mydataSig_5h_full<-mutate(mydataSig_5h_full, AsterCLV3_5h=ifelse(Treatment=="AsterCLV3_5h",1,0))
mydataSig_5h_full<-mutate(mydataSig_5h_full, AsterCLV3_L5S_5h=ifelse(Treatment=="AsterCLV3_L5S_5h",1,0))
mydataSig_5h_full<-mutate(mydataSig_5h_full, AtCLV3p_5h=ifelse(Treatment=="AtCLV3p_5h",1,0))
mydataSig_5h_full<-mutate(mydataSig_5h_full, AtCLV3p_S5L_5h=ifelse(Treatment=="AtCLV3p_S5L_5h",1,0))


mydataSig_5h_subset<-mydataSig_5h_full[c("Gene","Direction","AsterCLV3_5h","AsterCLV3_L5S_5h","AtCLV3p_5h","AtCLV3p_S5L_5h")]

#sig_overlapgraph<-lapply(mydataSig_treatment[1:3], function(x) {x$Gene})

#sig_overlapgraph_colors<-lapply(mydataSig_treatment[1:3], function(x) {x$Direction})

# get the treatment names from the columns
treatments<-colnames(mydataSig_5h_subset)[3:6]
# turn the binary to T or F
mydataSig_5h_subset[treatments]=mydataSig_5h_subset[treatments] ==1
# replace all F with NA
mydataSig_5h_subset<-mydataSig_5h_subset %>% replace_with_na_all(condition=~.x==FALSE)

# combine the columns where the gene ID is the same (keeping treatment identity intact)
mydataSig_5h_plot<-mydataSig_5h_subset %>% group_by(Gene)  %>% summarise_all(list(~ .[!is.na(.)][1]))
# convert the NA back to false (there has to be a better way to do this than
# converting the FALSE to NA, combining and then NA back to false...but I cannot figure it out)
mydataSig_5h_plot[is.na(mydataSig_5h_plot)] <- FALSE
# create a direction variable (stacked bar chart colored by column)
Direction<-colnames(mydataSig_5h_plot)[2]
# plot the upset 
png("plots/5h_up_AND_down.png", width=900, height=700, res=300)
ComplexUpset::upset(mydataSig_5h_plot,
                    treatments,name="Treatment",
                    base_annotations=list
                    ("Number of Intersecting Genes"=intersection_size
                      (counts=TRUE,mapping=aes(fill=Direction))
                      +scale_fill_manual(
                        values=c('Up'='#009E73','Down'='#D55E00','Difference'='#CC79A7'))+scale_y_continuous(expand=expansion(mult=c(0,0.1)))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()
                                                                                                                                                    ,axis.line=element_line(colour='black'))),
                    set_sizes=upset_set_size(geom=geom_bar(width = 0.4))+theme(axis.line.x = element_line(colour = 'black'),axis.ticks.x =element_line() ),
                    themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=15)))))
dev.off()


# colored upset...no legend
png("plots/5h_up_AND_down_white_stripe_no_legend.png", width=2700, height=2100, res=300)
ComplexUpset::upset(mydataSig_5h_plot,
                    treatments,name="Treatment",stripes="white",
                    base_annotations=list
                    ("Number of Intersecting Genes"=intersection_size
                      (counts=TRUE,mapping=aes(fill=Direction))
                      +scale_fill_manual(
                        values=c('Up'='goldenrod','Down'='#D55E00','Difference'='#40B0A6'), guide="none")+scale_y_continuous(expand=expansion(mult=c(0,0.1)))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()
                                                                                                                                                      ,axis.line=element_line(colour='black'))),
                    set_sizes=upset_set_size(geom=geom_bar(width = 0.4))+theme(axis.line.x = element_line(colour = 'black'),axis.ticks.x =element_line() ),
                    themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=15)))))
dev.off()


# plot the upset, but without the stacked bar chart 
png("plots/5h_complexupset.png", width=2700, height=2100, res=300)
ComplexUpset::upset(mydataSig_5h_plot,
                    treatments,name="Treatment",stripes="white",
                    base_annotations=list
                    ("Number of Intersecting Genes"=intersection_size
                      (counts=TRUE)
                      +scale_y_continuous(expand=expansion(mult=c(0,0.1)))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()
                                                                                                                                                    ,axis.line=element_line(colour='black'))),
                    set_sizes=upset_set_size(geom=geom_bar(width = 0.4))+theme(axis.line.x = element_line(colour = 'black'),axis.ticks.x =element_line() ),
                    themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=15)))))


dev.off()

genes_true <- mydataSig_5h_plot$Gene[mydataSig_5h_plot$Direction == "Difference" & mydataSig_5h_plot$AsterCLV3_5h == TRUE]
genes_true

genes_true <- mydataSig_5h_plot$Gene[mydataSig_5h_plot$AtCLV3p_5h == TRUE & mydataSig_5h_plot$AsterCLV3_5h == TRUE & mydataSig_5h_plot$AsterCLV3_L5S_5h==FALSE]
genes_true
genes_true <- mydataSig_5h_plot$Gene[mydataSig_5h_plot$AtCLV3p_5h == FALSE & mydataSig_5h_plot$AsterCLV3_5h == FALSE & mydataSig_5h_plot$AsterCLV3_L5S_5h==FALSE & mydataSig_5h_plot$AtCLV3p_S5L_5h==TRUE & mydataSig_5h_plot$Direction == "Up"]
genes_true





# (2) LOOK AT DATA BETWEEN @ 24H
# read in the data
DEData_24h <- ImportCSVs('deseq_results/x24h/', 0.05)
# filter out significant results
mydataSig_24h <- lapply(DEData_24h, SigDEdf, PvaluesCol=7, CritP=0.05)

# remove result_AtCLV3p_S5L_24h since it has no DE genes
mydataSig_24h <- mydataSig_24h[!(names(mydataSig_24h) %in% "result_AtCLV3p_S5L_24h")]

# different list of df for up or down expression
dataSigup_24h <- lapply(mydataSig_24h, MoreCritNum, column=3, critNum=0)
dataSigdown_24h <- lapply(mydataSig_24h, LessCritNum, column=3, critNum=0)
dataSigup_24h <- dataSigup_24h[!(names(dataSigup_24h) %in% "result_AsterCLV3_24h")]

# add a column labeling "up" or "down" to each list of dataframes 
dataSigup_24h <- lapply(dataSigup_24h, transform, Direction="Up")
dataSigdown_24h <- lapply(dataSigdown_24h, transform, Direction="Down")

# Get the common names in both lists
common_names <- intersect(names(dataSigup_24h), names(dataSigdown_24h))
# Initialize an empty list to store the result
direction_data_sig_24h <- list()

# Combine data frames with the same name
for (name in common_names) {
  direction_data_sig_24h[[name]] <- rbind(dataSigup_24h[[name]], dataSigdown_24h[[name]])
}
# Add the additional dataframe from dataSigdown_24h
direction_data_sig_24h$result_AsterCLV3_24h <- dataSigdown_24h$result_AsterCLV3_24h

# rename the dataframes in the list
names(direction_data_sig_24h) <- c("AsterCLV3_L5S_24h", "AtCLV3p_24h", "AsterCLV3_24h")

# subset the dataframes to only contain the gene ID and the direction of expression
gene_direction_only <- lapply(direction_data_sig_24h, function(x) { x[c("Gene", "Direction")] })
gene_direction_only <- gene_direction_only %>% purrr::reduce(full_join, by=c('Gene'))
# combine all three columns of direction into one string with where direction 
# is outlined. For example. "Up Up Down" If the gene is not in a treatment NA is used
gene_direction_only$dif <- paste(gene_direction_only$Direction.x, gene_direction_only$Direction.y, gene_direction_only$Direction)

# loop through all of the genes present and add the genes to the "up", "down",
# or "difference" lists depending on the directionality of expression

difference = c()
up = c()
down = c()
for (i in 1:nrow(gene_direction_only)) {
  if (grepl('Down', gene_direction_only[i, "dif"]) & grepl("Up", gene_direction_only[i, "dif"])) {
    difference <- append(difference, gene_direction_only[i, "Gene"])
  } else if (grepl('Up', gene_direction_only[i, "dif"])) {
    up <- append(up, gene_direction_only[i, "Gene"])
  } else {
    down <- append(down, gene_direction_only[i, "Gene"])
  }
}
#gene_direction_only <- names(gene_direction_only)
# for each list created in the above loop, turn it into a df and then rbind them into one df
#difference <- data.frame(difference)
#difference <- cbind(difference, x="Difference")
#colnames(difference)[1] <- "Gene"
# difference is EMPTY...none show a change in direction that is why it is commented out
up <- data.frame(up)
up <- cbind(up, x="Up")
colnames(up)[1] <- "Gene"
down <- data.frame(down)
down <- cbind(down, x="Down")
colnames(down)[1] <- "Gene"
# df with Gene ID as column 1 and direction of expression as column 2 
up_or_down <- rbind(down, up)

# in our initial dataframe containing all treatments, add the directionality of the gene
direction_data_sig_24h$AsterCLV3_24h$Direction <- up_or_down$x[match(direction_data_sig_24h$AsterCLV3_24h$Gene, up_or_down$Gene)]
direction_data_sig_24h$AsterCLV3_L5S_24h$Direction <- up_or_down$x[match(direction_data_sig_24h$AsterCLV3_L5S_24h$Gene, up_or_down$Gene)]
direction_data_sig_24h$AtCLV3p_24h$Direction <- up_or_down$x[match(direction_data_sig_24h$AtCLV3p_24h$Gene, up_or_down$Gene)]
#direction_data_sig_24h$AtCLV3p_S5L_24h$Direction <- up_or_down$x[match(direction_data_sig_24h$AtCLV3p_S5L_24h$Gene, up_or_down$Gene)]

# create a column in each df that contains a treatment variable
direction_data_sig_24h$AsterCLV3_24h <- cbind(direction_data_sig_24h$AsterCLV3_24h, Treatment="AsterCLV3_24h")
direction_data_sig_24h$AsterCLV3_L5S_24h <- cbind(direction_data_sig_24h$AsterCLV3_L5S_24h, Treatment="AsterCLV3_L5S_24h")
direction_data_sig_24h$AtCLV3p_24h <- cbind(direction_data_sig_24h$AtCLV3p_24h, Treatment="AtCLV3p_24h")
#direction_data_sig_24h$AtCLV3p_S5L_24h <- cbind(direction_data_sig_24h$AtCLV3p_S5L_24h, Treatment="AtCLV3p_S5L_24h")

#direction_data_sig_treatment_full <- direction_data_sig_treatment %>% reduce(full_join, by=c('Treatment'))
# join the df by gene ID
#direction_data_sig_treatment_full <- direction_data_sig_treatment %>% purrr::reduce(full_join, by=c('Gene'))

# combine all the dataframes on the rows
direction_data_sig_24h_full <- bind_rows(direction_data_sig_24h$AsterCLV3_24h, direction_data_sig_24h$AsterCLV3_L5S_24h, direction_data_sig_24h$AtCLV3p_24h)
direction_data_sig_24h_full[,"Gene"]
comparisons <- c("AsterCLV3_24h", "AsterCLV3_L5S_24h", "AtCLV3p_24h")
treatment <- c("Treatment")

# one-hot encode the treatment variables -- why did I do this?
direction_data_sig_24h_full <- mutate(direction_data_sig_24h_full, AsterCLV3_24h = ifelse(Treatment=="AsterCLV3_24h", 1, 0))
direction_data_sig_24h_full <- mutate(direction_data_sig_24h_full, AsterCLV3_L5S_24h = ifelse(Treatment=="AsterCLV3_L5S_24h", 1, 0))
direction_data_sig_24h_full <- mutate(direction_data_sig_24h_full, AtCLV3p_24h = ifelse(Treatment=="AtCLV3p_24h", 1, 0))

direction_data_sig_24h_subset <- direction_data_sig_24h_full[c("Gene", "Direction", "AsterCLV3_24h", "AsterCLV3_L5S_24h", "AtCLV3p_24h")]

#sig_overlapgraph <- lapply(direction_data_sig_treatment[1:3], function(x) { x$Gene })

#sig_overlapgraph_colors <- lapply(direction_data_sig_treatment[1:3], function(x) { x$Direction })

# get the treatment names from the columns
treatments <- colnames(direction_data_sig_24h_subset)[3:5]
# turn the binary to T or F
direction_data_sig_24h_subset[treatments] = direction_data_sig_24h_subset[treatments] == 1
# replace all F with NA
direction_data_sig_24h_subset <- direction_data_sig_24h_subset %>% replace_with_na_all(condition = ~.x == FALSE)

# combine the columns where the gene ID is the same (keeping treatment identity intact)
direction_data_sig_24h_plot <- direction_data_sig_24h_subset %>% group_by(Gene) %>% summarise_all(list(~ .[!is.na(.)][1]))
# convert the NA back to false (there has to be a better way to do this than
# converting the FALSE to NA, combining and then NA back to false...but I cannot figure it out)
direction_data_sig_24h_plot[is.na(direction_data_sig_24h_plot)] <- FALSE
# create a direction variable (stacked bar chart colored by column)
Direction <- colnames(direction_data_sig_24h_plot)[2]
# plot the upset 
png("plots/24h_up_AND_down.png", width=900, height=700, res=300)
ComplexUpset::upset(direction_data_sig_24h_plot,
                    treatments, name="Treatment",
                    base_annotations = list(
                      "Number of Intersecting Genes" = intersection_size(
                        counts = TRUE, mapping = aes(fill = Direction))
                      + scale_fill_manual(
                        values = c('Up' = '#009E73', 'Down' = '#D55E00', 'Difference' = '#CC79A7'))
                      + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
                      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = 'black'))
                    ),
                    set_sizes = upset_set_size(geom = geom_bar(width = 0.4)) + theme(axis.line.x = element_line(colour = 'black'), axis.ticks.x = element_line()),
                    themes = upset_modify_themes(list('intersections_matrix' = theme(text = element_text(size = 15)))))
dev.off()

# plot the upset, but without the stacked bar chart 
png("plots/24_complexupset.png", width=2700, height=2100, res=300)
ComplexUpset::upset(direction_data_sig_24h_plot,
                    treatments, name="Treatment", stripes="white",
                    base_annotations = list(
                      "Number of Intersecting Genes" = intersection_size(
                        counts = TRUE)
                      + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
                      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = 'black'))
                    ),
                    set_sizes = upset_set_size(geom = geom_bar(width = 0.4)) + theme(axis.line.x = element_line(colour = 'black'), axis.ticks.x = element_line()),
                    themes = upset_modify_themes(list('intersections_matrix' = theme(text = element_text(size = 15)))))
dev.off()


# can subset and look at differnet sets of genes using this
genes_true <- mydataSig_24h_plot$Gene[mydataSig_24h_plot$AtCLV3p_24h == TRUE & mydataSig_24h_plot$AsterCLV3_24h == TRUE & mydataSig_24h_plot$AsterCLV3_L5S_24h==FALSE]
