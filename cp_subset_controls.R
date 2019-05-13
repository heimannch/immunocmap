#List of controls to be used for subsetting signatures of controls
#Returns: list of cell lines with controls for compounds available (list_ctl).
#Subset of sig_ids that match the cell lines with available controls (sigs)



source('load_data.R')
library(stringr)
library(skimr)
library(plyr)
library(dplyr)

# 1. Selection of controls for analysis of COMPOUNDS!
# select consensus controls: "ctl_vehicle.cns" and "ctl_trt.cns"

sig_ctl <- siginfo %>% filter(stringr::str_detect(pert_type,"ctl_") & str_detect(pert_type,".cns"))
sig_ctl <- sig_ctl %>% filter(pert_type != "ctl_vector.cns" )

#select tissues that have samples from both normal and tumor cells

#adding the information about cells
sig_ctl <- merge(sig_ctl, cellinfo, by = "cell_id")

#list of tissues with controls of tumor and normal cells
aux_table <- sig_ctl %>% group_by(primary_site, sample_type) %>%
  filter(sample_type == "tumor" | sample_type == "normal") %>% group_by(primary_site)

#getting the list of sig_ids for controls
list_ctl <- subset(sig_ctl, primary_site %in% aux_table$primary_site, select = c(sig_id, pert_iname, cell_id, primary_site, sample_type, pert_time))

#-------------------------------------------------------------------------------------------------------------------------------#

#2. Selection of signatures that match the controls
# - select pertubations with compounds - "trt_cp"
# - select signatures that match set of controls by: cell id and pert_itime
# If two different pertubagens in the control are available, keep both
#Only include signatures with more than 3 replicates

#adding important columns in the sigmetrics dataset
index <- match(sigmetrics$sig_id, siginfo$sig_id)
siginfo$distil_nsample <- sigmetrics$distil_nsample[index]

sigs <- subset(siginfo, pert_type == "trt_cp" & cell_id %in% list_ctl$cell_id & pert_time %in% list_ctl$pert_time & distil_nsample>=3)

#------------------------------------------------------------------------------------------------------------------------------#
#3. Threshold of "modulating z-score" for each cell line and controls
#given a cell line, extract the associated consensus controls, and determine the maximum z-score for a given quantile

get_thresholds <- function(cell = "MCF7", list_ctls = list_ctl, pert = FALSE){
  #in case the threshold should be computed using only one kind of control, inform name (options: "DMSO", "UnTrt", "H2O")
  #The default ignores controls with H2O (consistently different results, and available only for a few cell lines)  
  
  if (pert == FALSE){
    ctl_cell <- subset(list_ctls, cell_id == cell & pert_iname != "H2O")
  } else{
    ctl_cell <- subset(list_ctls, cell_id == cell & pert_iname == pert)  
  }
  
  
  #getting the gene expression data for all controls associated with this cell_id
  ds <- cmapR::parse.gctx(ds_path,
                          cid=ctl_cell$sig_id)
  
  expr_data <- as.data.frame(ds@mat)
  quantile <- data.frame(sig_id= character(ncol(expr_data)), pert_iname = character(ncol(expr_data)), first = numeric(ncol(expr_data)), last = numeric(ncol(expr_data)),
                           stringsAsFactors = FALSE)
  
  for (i in 1:ncol(expr_data)){
    
    sig <- as.character(ctl_cell[i,1])
    
    quantile[i,1] <- sig  
    quantile[i,2] <- as.character(ctl_cell[i,1])
    quantile[i,3:4] <-  quantile(expr_data[,i] , c(.01, .99))
    
    
    i = i+1
  }
  return(quantile)
  
}

get_z_threshold <- function(cell_id, list_ctls = list_ctl, pert = FALSE){
  ctls_data <- get_thresholds(cell = cell_id, list_ctls, pert)
  
  #computing the mean of the 1% z-scores threshold across all controls
  zs <- as.data.frame(abs(c(ctls_data$first, ctls_data$last)))
  return(skim_to_wide(zs)) 
}

#Creating a list of the summary statistics of the 1% and 99% quantiles of z-scores of controls for a given cell line
get_z_list <- function(list_cells = "MCF7", tissue = FALSE, list_ctls = list_ctl, pert = FALSE){
#by default, it is returned a list with all cell lines with controls.
#The user needs to input a cell line name, but it can be changed by supplying a primary site name

#Geting the list of cell_id requested  
  if(tissue != FALSE){
    cells <- subset(list_ctls, primary_site %in% tissue, cell_id)
    cells <- unique(cells$cell_id)
  # }else if(any(list_cells == "all")){
  #   cells <- subset(cellinfo, cell_id %in% list_ctl$cell_id)
  #   cells <- unique(cells$cell_id)
  }else{
    cells <- unique(list_cells) 
  }
  
  
  z_control <- plyr::adply(matrix(cells), 1, get_z_threshold)
  z_control$variable <- cells
  z_control$X1 <- NULL
  
  return(z_control)
}


#------------------------------------------------------------------------------------------------------------------------------#
#4. Threshold of Signature Strength for each cell line

get_ss_threshold <- function(cell, list_ctls = list_ctl, pert = FALSE){
  
  if (pert == FALSE){
    ctl_cell <- subset(list_ctls, cell_id == cell & pert_iname != "H2O")
  } else{
    ctl_cell <- subset(list_ctls, cell_id == cell & pert_iname == pert)  
  }
  
  modulating <- data.frame(sig_id= character(nrow(ctl_cell)), 
                         ngenes_modulated_dn_lm = numeric(nrow(ctl_cell)),  ngenes_modulated_up_lm = numeric(nrow(ctl_cell)), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(ctl_cell)){
    
    sig <- as.character(ctl_cell[i,1])
    
    modulating[i,1] <- sig  
    lm_mod <-  subset(sigmetrics, sig_id == ctl_cell[i,1], c(ngenes_modulated_dn_lm, ngenes_modulated_up_lm))
    modulating[i,2:3] <- lm_mod
    
    i = i+1
  }
  return(modulating)
  
}


# test <- get_ss_threshold(cell = "A375")
# print(test)






