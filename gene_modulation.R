#Functions to query z-scores (Level 5 data), organize it in tidy formats, indexed by cell line

source('load_data.R')
source('cp_subset_controls.R')

library(cmapR)
library(data.table)
library(pheatmap)
library(repr)
library(ggplot2)
library(tidyr)



#function to grab signatures with the same pertubagem in two different times
get_timelapse <- function(time1, time2, n_replicates){

  #adding important columns in the sigmetrics dataset
  index <- match(sigmetrics$sig_id, siginfo$sig_id)
  sigmetrics$cell_id <- siginfo$cell_id[index]
  sigmetrics$pert_time <- siginfo$pert_time[index]
  sigmetrics$pert_itime <- siginfo$pert_itime[index]
  
  pert_time1 <- sigmetrics[distil_nsample>=n_replicates & pert_time == time1]$pert_id
  pert_time2 <- sigmetrics[distil_nsample>=n_replicates & pert_time == time2]$pert_id
  
  pert_time12 <- intersect(pert_time1, pert_time2)
  sig_time12 <- sigmetrics[pert_id %in% pert_time12 & (pert_time == time1 | pert_time == time2), .(sig_id, pert_id, pert_iname,pert_type, pert_itime)]
  
  return(sig_time12)
  
}

#getting a subset of signatures, based on the number of replicates and time of interest
get_onetime_sig <- function(n_replicates, time){
  #time must be in number of hours. The CMap dataset have data only for the following times:  1 h 2 h 3 h 4 h 6 h 24 h 48 h 72 h 96 h 120 h 144 h 168 h 
  
  #adding important columns in the sigmetrics dataset
  index <- match(sigmetrics$sig_id, siginfo$sig_id)
  sigmetrics$cell_id <- siginfo$cell_id[index]
  sigmetrics$pert_time <- siginfo$pert_time[index]
  sigmetrics$pert_itime <- siginfo$pert_itime[index]
  
  sigs_of_interest <- sigmetrics[distil_nsample>n_replicates & pert_time == time]
  sig_ids <- sigs_of_interest$sig_id
  
  return(sig_ids)
}
  
get_zscores <- function(gene_id, signatures){
  #gene_id must be the number in the "pr_gene_id" in the "geneinfo" table
  #Signatures must be a list only with the signatures (indexed by $sig_id in the metadata)
   
  #Querying the dataset, based on the gene of interest and the selected signatures
  ds <- cmapR::parse.gctx(ds_path,
                          rid=as.character(gene_id),
                          cid=signatures)
  
  tidy_ds <- as.data.frame(ds@mat) %>% mutate(pr_gene_id = row.names(ds@mat)) %>% tidyr::gather(key = sig_id, z_score, - pr_gene_id)
  
  return(tidy_ds)
}

#function returns the cell sample type (normal or tumor) for a list of signatures
get_cell_sample_types <- function(signatures){
  
  #getting the list of cell for each signature
  list_cells <- siginfo[sig_id %in% signatures$sig_id, .(sig_id, cell_id)]
  
  #getting the type of cell in each signature
  list_cells$sample_type <-  cellinfo[match(list_cells$cell_id, cellinfo$cell_id),"sample_type"] 
  
  return(list_cells)
}
  
#function returns a histogram comparing the distribution of z-scores for a given gene in two different time-frames
get_hist_zscore <- function(gene_id, signatures, by_cellsample_type = FALSE){. 
  #The signature table needs to contain the list of signatures coded as "sig_id" 
  #and the time of the experiment in a column called "pert_itime"
  
  #getting the gene expression data with the number of replicates and time of measure
  ds <- cmapR::parse.gctx(ds_path,
                          rid=as.character(gene_id),
                          cid=signatures$sig_id)
  
  my_gene <- geneinfo[pr_gene_id == gene_id]$pr_gene_symbol
  
  print(paste(i, "Getting info for the gene ", my_gene))
  #the signatures are the columns in the gene expression data. Let's make them become rows, and then we add a new column informing the time
  gene_expr <- as.data.frame(ds@mat)
  gene_expr <- t(gene_expr)
  
  #adding the column with the time
  gene_expr<- cbind(gene_expr, signatures$pert_itime)
  gene_expr <- as.data.frame(gene_expr)
  colnames(gene_expr) <- c("zscore", "time")
  gene_expr$zscore <- as.numeric(levels(gene_expr$zscore))[gene_expr$zscore]
  
  if (by_cellsample_type == TRUE){
    #getting the sample type for the signatures
    type_cell <- get_cell_sample_types(signatures)
    gene_expr <- cbind(gene_expr, sample_type = type_cell$sample_type)
    
    #some NAs values are coded as -666, so, in order to improve the grouping in plotting, let's make it only one group
    gene_expr$sample_type[gene_expr$sample_type == -666] <- NA
    #Now, let's plot a histogram with the z-score distribution of the gene of interest, subsetting it to different times and colors based on the sample_type of the cell
    hist_plot <- ggplot(gene_expr, aes(y=zscore, fill = sample_type))+
      geom_boxplot()+
      facet_wrap(~time)+
      ggtitle(my_gene)
    
  }else{
    #Now, let's plot a histogram with the z-score distribution of the gene of interest, subsetting it to different times
    hist_plot <- ggplot(gene_expr, aes(y=zscore))+
      geom_boxplot()+
      facet_wrap(~time)+
      ggtitle(my_gene)
  }
  
  print(hist_plot)
  
}
  

get_dist_zscore_onetime <- function(gene_id, signatures){
  #function returns a histogram comparing the distribution of z-scores for a given gene, regardless of the time of the experiment. 
  #The signature table needs to contain the list of signatures coded as "sig_id" 
  #and the time of the experiment in a column called "pert_itime"
  
  #getting the gene expression data with the number of replicates and time of measure
  ds <- get_zscores(gene_id, signatures)
  
  my_gene <- geneinfo[pr_gene_id == gene_id]$pr_gene_symbol
  
  #Now, let's plot a boxplot with the z-score distribution of the gene of interest
  options(repr.plot.width = 5.5, repr.plot.height = 4)
  hist(ds@mat, col = "dodgerblue", main = as.character(my_gene), xlab="z-score")
  abline(v=0.54, lty=2, col=2, lwd=1.3 )
  abline(v=-0.54, lty=2, col=2, lwd=1.3 )
}


#-------Functions to select signatures based on a threshold and cell line-----------------------#

#1. Create a list of signatures with up, down and non modulated z-scores.

get_mod_sigs <- function(gene_id, signatures, mod_thres = 2){
  #mod_thres - z-score used as a threshold to determine if a gene is modulated by the signature - the higher the more modulated. 
  #The negative value is used to threshold down-modulation.
  #non_mod_thres - z-score used to determine if a gene is NOT modulated by the signature - "non-regulation" increases when closer to zero
  
  tidy_ds <- get_zscores(gene_id, signatures)
 
  #getting the signatures that significantly modulated the gene, assuming that it is represented by an absolute z-score greater than mod_thres
  
  tidy_mod <- tidy_ds %>% mutate(mod = case_when(z_score> mod_thres ~ "up",
                                                z_score< (-mod_thres) ~ "dn",
                                                z_score > (-mod_thres) | z_score < mod_thres ~ "no"))
  
  return(tidy_mod)
  
}


#Returns the list of annotation of up and down modulation, indexed by the cell line
get_modsigs_per_cellid <- function(controls_stats, list_genes, signatures = siginfo){

  sig_mod <- data.frame()
  
  #sig_list_cell <- vector('list', nrow(controls_stats))

  for (i in 1:nrow(controls_stats)){

    cell <- as.character(controls_stats[i,2])
    threshold <- subset(controls_stats, variable == cell, mean)
    sigs <- subset(signatures, cell_id == cell, sig_id)


    sig_mod <- rbind(sig_mod, get_mod_sigs(gene_id = list_genes, signatures = sigs$sig_id, mod_thres = as.numeric(threshold)))
    i=i+1
  }
  
  #names(sig_list_cell) <- as.character(controls_stats$variable)
  # This list is organized as follows:
  # 1st level: cell line
  
  
  return(sig_mod)
  
}


#2. This function calls the above for each gene that you want to check for modulation and returns a list
##mod_list[[gene_id]][[type]][[sig_id]] - NOT NECESSARY ANYMORE!
# 
# get_list_modsigs <- function(list_genes, signatures, mod_thres = 3, non_mod_thres= 0.004){
#   
#   sig_list <- vector('list', nrow(signatures))
#   
#   for(i in 1:nrow(list_genes)){
#     sig_list[[i]] <- get_mod_sigs(list_genes[i,1], signatures$sig_id, mod_thres, non_mod_thres)
#     i = i+1 
#   }
#   
#   #removing the NULL elements
#   sig_list <-  sig_list[-which(sapply(sig_list, is.null))]
#   
#   return(sig_list)
# }


#UNNECESSARY
# make_df_modulators <- function(modulation_list, controls_stats, genes, type = "up_mod"){
# #options para "type" are: "up_mod", "dn_mod", "non_mod"
#   upperbound <- length(unlist(modulation_list))
#   
#   mod_df <- data.frame(sig_id = character(upperbound), gene_id = character(upperbound), cell_id = character(upperbound), time = factor(upperbound), stringsAsFactors=FALSE)
#   
#   new_row <- 1
# 
# #Structure of the list:  
# #mod_list[[cell_id]][[gene_id]][[type]][[sig_id]]
#   
# for (c in 1:nrow(controls_stats)){
#   
#   for (i in 1:length(modulation_list[[c]])){
#     print(paste("getting info for gene ", i))
#     #accessing only the modulated signatures
#     for(j in 1:length(names(modulation_list[[c]][[i]][[type]]))){
#       mod_df[new_row, 1] <- names(modulation_list[[c]][[i]][[type]][j])
#       mod_df[new_row, 2] <- as.character(genes[i,2])
#       mod_df[new_row, 3] <- as.character(controls_stats[c,2])
#       
#       new_row <-  new_row + 1
#       j <-  j+1
#     }
#     
#     i <-  i+1
#   }
#   
#   c <- c+1
# }
#   
#   mod_df$time <- siginfo[match(mod_df$sig_id, siginfo$sig_id),"pert_itime"]
#   mod_df <- mod_df[!apply(mod_df == "", 1, all),]
# 
#   return(mod_df)
# }

#cell_types <- get_cell_sample_types(mod_df)
#mod_df$sample_type <-  cell_types[match(mod_df$sig_id, cell_types$sig_id),"sample_type"]
#mod_df$time <- siginfo[match(mod_df$sig_id, siginfo$sig_id),"pert_itime"]


