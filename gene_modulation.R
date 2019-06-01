#Functions to query z-scores (Level 5 data), organize it in tidy formats, indexed by cell line

source('load_data.R')
source('cp_subset_controls.R')

library(cmapR)
library(data.table)
library(pheatmap)
library(repr)
library(ggplot2)
library(tidyr)



#' Function to query for expression data in the Connectivity Map dataset with output in tidy format.
#' 
#' @param gene_id The id of the gene of interest to query. It must be available from the `geneinfo` dataframe (loaded in session by the `load_data.R` script)  
#' @param signatures The id of the signatures of interest. It must be available from the `siginfo` dataframe (loaded in session by the `load_data.R` script)
#' @param n_replicates The minimal number of replicates a signature needs to have to be included in the analysis
#' @return Dataframe in tidy format with each row containing the signature of the experiment, a gene id and the corresponding z-score
#' @example 
#' get_z_score(1493, siginfo[10:30, 1])

get_zscores <- function(gene_id, signatures){
   
  #Querying the dataset, based on the gene of interest and the selected signatures
  ds <- cmapR::parse.gctx(ds_path,
                          rid=as.character(gene_id),
                          cid=signatures)
  
  tidy_ds <- as.data.frame(ds@mat) %>% mutate(pr_gene_id = row.names(ds@mat)) %>% tidyr::gather(key = sig_id, z_score, - pr_gene_id)
  
  return(tidy_ds)
}


#-------Functions to select signatures based on a threshold and cell line-----------------------#


#' Create a dataframe of signatures with up, down and non modulated z-scores.
#' 
#' @param gene_id The id of the gene of interest to query. It must be available from the `geneinfo` dataframe (loaded in session by the `load_data.R` script)  
#' @param signatures The id of the signatures of interest. It must be available from the `siginfo` dataframe (loaded in session by the `load_data.R` script)
#' @param mod_thres The threshold of z-score to consider up or down modulation. Absolute values higher that the threshold are considered modulated, while lower values are not-modulated
#' @return Dataframe in tidy format with each row containing the signature of the experiment, a gene id, the corresponding z-score, and if the gene is up, down or not modulated
#' @example 
#' get_mod_sigs(1493, siginfo[10:30, 1], 1.5)
#' 
get_mod_sigs <- function(gene_id, signatures, mod_thres = 2){
  #mod_thres - z-score used as a threshold to determine if a gene is modulated by the signature - the higher the more modulated. 
  #The negative value is used to threshold down-modulation.
  
  tidy_ds <- get_zscores(gene_id, signatures)
 
  #getting the signatures that significantly modulated the gene, assuming that it is represented by an absolute z-score greater than mod_thres
  
  tidy_mod <- tidy_ds %>% mutate(mod = case_when(z_score> mod_thres ~ "up",
                                                z_score< (-mod_thres) ~ "dn",
                                                z_score > (-mod_thres) | z_score < mod_thres ~ "no"))
  
  return(tidy_mod)
  
}


#' Create a dataframe of signatures with up, down and non modulated z-scores, indexed by cell line.
#' 
#' Each cell line may have a different z-score as threshold criteria for modulation. This function is useful to compute and organize the data considering these differences
#' 
#' @param controls_stats Dataframe with the threshold of modulation for each cell line. The column with the threshold information should be named "mean".  
#' @param list_genes List of genes to query against the Connectivity Map dataset. 
#' @param signatures The id of the signatures of interest. It must be available from the `siginfo` dataframe (loaded in session by the `load_data.R` script)
#' @return Dataframe in tidy format with each row containing the signature of the experiment, the cell line, a gene id, the corresponding z-score, and if the gene is up, down or not modulated
#' 

get_modsigs_per_cellid <- function(controls_stats, list_genes, signatures = siginfo){

  sig_mod <- data.frame()
  
  for (i in 1:nrow(controls_stats)){

    cell <- as.character(controls_stats[i,2])
    threshold <- subset(controls_stats, variable == cell, mean)
    sigs <- subset(signatures, cell_id == cell, sig_id)


    sig_mod <- rbind(sig_mod, get_mod_sigs(gene_id = list_genes, signatures = sigs$sig_id, mod_thres = as.numeric(threshold)))
    i=i+1
  }
  
  return(sig_mod)
  
}



#---deprecated functions - check if I want to keep them ---#

#' Selects signatures with the same pertubagens in two different time frames. The names of the objects of manipulation are inherited from 'load_data.R' 
#' This function is useful in case there is interest in comparing the effect of time in gene expression levels and provides an filtering stage based on particular stages of interest.
#' The CMap dataset have data only for the following times:  1 h 2 h 3 h 4 h 6 h 24 h 48 h 72 h 96 h 120 h 144 h 168 h
#' 
#' @param time1 Number representing one of the time frames to be considered in the analysis 
#' @param time2 Number representing the second of the time frames to be considered in the analysis
#' @param n_replicates The minimal number of replicates a signature needs to have to be included in the analysis
#' @return Dataframe with the signatures with pertubations in the two time frames provided 
#' @example 
#' get_timelapse(6,24,3)

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
