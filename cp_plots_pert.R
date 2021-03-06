#Code to make plots of differential expression across genes and compounds

source('gene_modulation.R')
library(pheatmap)
library(ggrepel)


#------------------------------------------------General Functions-----------------------------------------------------------#

#' Merges the siginfo and cellinfo dataframes based on the sig_id.
#' 
#' This function is useful for downstream analysis that consider data related to perturbations and cell info.
#' @param list_of_sigs dataframe with the signature of interest. May contain other information, but signatures need to be in a column named sig_id.
#' @return Dataframe with all columns present in the siginfo and cellinfo dataframes, including additional columns in the original dataset
#' 

make_full_table <- function(list_of_sigs){
  
  df_full <- merge(list_of_sigs, siginfo, by = "sig_id")
  df_full$time <- NULL
  df_full$sample_type <- NULL
  
  df_full <- merge(df_full, cellinfo, by = "cell_id")
  
  return(df_full)
}


#Returns the number of times and % that a gene is up and down modulated in a combination of cell line and compound or pertubagen class (PCL)
make_freq_mod <- function(tidy_mod, name_cp_pcl, by.pcl = FALSE, pcl_annot = FALSE){
  
  pert_info <- siginfo %>% 
    select(sig_id, pert_iname, cell_id)
  
  tidy_mod_full <- merge(tidy_mod, pert_info, by = "sig_id")
  tidy_mod_full <- merge(tidy_mod_full, tcga_genes, by = "pr_gene_id")
  
  if(by.pcl != FALSE){
    name_cp_pcl <- (subset(pcl_annot, pcl_id %in% name_cp_pcl, pert_iname))
    name_cp_pcl <- as.vector(name_cp_pcl$pert_iname)
  }
  
    freq_dn <- merge((tidy_mod_full %>% filter(pert_iname %in% name_cp_pcl) %>% 
                        group_by(cell_id, pr_gene_symbol, pert_iname) %>% 
                        summarise(total = n_distinct(sig_id))), 
                     tidy_mod_full %>% filter(mod == "dn" & pert_iname%in% name_cp_pcl)%>% 
                       group_by(cell_id, pr_gene_symbol, pert_iname) %>% 
                       summarise(n_dn = n_distinct(sig_id)), 
                     by= c("cell_id", "pr_gene_symbol", "pert_iname"), all.x = TRUE) %>% 
      mutate(freq_dn = n_dn/(total))
    
    
    freqs <- merge(freq_dn, 
                   (tidy_mod_full %>% filter(mod == "up" & pert_iname %in% name_cp_pcl) %>% 
                      group_by(cell_id, pr_gene_symbol, pert_iname) %>% 
                      summarise(n_up = n_distinct(sig_id))), 
                   by =c("cell_id", "pr_gene_symbol", "pert_iname"), all.x = TRUE) %>% 
      mutate(freq_up = n_up/(total)) %>% arrange(cell_id) 
    
  
  freqs[is.na(freqs)] <- 0
  
  #Computing the difference between frequency of up and down expression
  freqs <- freqs %>% mutate(diff = freq_up - freq_dn)
  return(freqs)
}

#the summary statistics of z-scores for a given gene across all compounds and returns the list of compounds 
#that generate the highest mean of up or down expression
get_top_pert_pergene <- function(tidy_mod, gene_symbol, pert_list = "all", stats_selector = "mean", by_cell = TRUE, n_pert = as.numeric(nrow(tidy_mod))){
  
  #pert_list: "all" selects all perturbations that are associated with the signatures in tidy_mod. User can define a vector with names of compounds 
  gene_id <- subset(tcga_genes, pr_gene_symbol == gene_symbol, pr_gene_id)
  
  #Including pert_iname and cell_id data on the tidy_mod table
  pert_cell_info <- siginfo %>% select(sig_id, pert_iname, cell_id)
  tidy_mod_full <- merge(pert_cell_info, tidy_mod, by= "sig_id")
  
  #selecting all compounds associated with the signatures provided
  if(pert_list == "all"){
    pert_list <- as.vector(tidy_mod_full$pert_iname)
  }
  
  #Organizing the data manipulation chosen with the "by_cell" parameter
  
  if(by_cell == TRUE){
    group_columns1 <- c("pr_gene_id", "cell_id", "pert_iname", "z_score")
    group_columns2 <- c("pr_gene_id", "cell_id", "pert_iname")
  }else{
    group_columns1 <- c("pr_gene_id", "pert_iname", "z_score")
    group_columns2 <- c("pr_gene_id", "pert_iname")
  }
  
  #Calculating the summary statistics for each combination of gene_id and pert_iname
  stats_pert <- tidy_mod_full %>% filter(pr_gene_id %in% gene_id & pert_iname %in% pert_list) %>% 
    select_(.dots = group_columns1) 
  
  stats_pert <- stats_pert %>% 
    group_by_(.dots = group_columns2) %>% skimr::skim_to_wide(z_score)
  
  #the columns are in character format, so let's change them to numeric
  stats_cols <- c("n", "mean", "sd", "p0", "p25", "p50", "p75", "p100")
  stats_pert[stats_cols] <- lapply(stats_pert[stats_cols], as.numeric)
  
  
  final_cols <- c(group_columns2, stats_cols)
  #Selecting the combination of cell lines and perturbations that have the highest absolute value for the selected statistic
  top_pert <- stats_pert %>% filter(n >2) %>% group_by_(match(group_columns2,names(tidy_mod_full))) %>% 
    arrange(desc(abs(get(stats_selector)))) %>% ungroup() %>% 
    select_(.dots = final_cols) %>% head(as.numeric(n_pert))
  
  return(top_pert)
}


#function to add annotation to an existing table (saved as pclsCustom.csv in the data folder)
add_annot <- function(annot_table, pert_id, pert_iname, pcl_id,  pcl_name, pcl_type, pcl_source, pcl_criteria){
  #example of addition of compound annotation:
  #add_annot(annot_table, "BRD-K51313569", "palbociclib", "CP_CDK46_INHIBITOR","CDK46 INHIBITOR", "CP","Drugbank","MoA")
  
  annot_table <- rbind(annot_table, c(pert_id, pert_iname, pcl_id, pcl_name, pcl_type, pcl_source, pcl_criteria))
  #update the original csv file, so this annotation is preserved for future use
  write.csv(annot_table, file = "data/pclsCustom.csv", quote = FALSE, sep = ",")
  return(annot_table)
}

#-------PLOTTING FUNCTIONS-------#


#This function makes a scatterplot of the number of signatures associated with a drug or PCL that are upmodulating a given gene, including labels for extremes.
plot_count_mod <- function(tidy_mod, name_cp_pcl, by.pcl = TRUE, pcl_annot = FALSE, x_scale, y_scale, x_label, y_label){
  #if  by.pcl = FALSE, provide a compound name to subset the data
  #provide x_scale and y_scale in the format: c(0,170)
  
  sigs_cp <- make_freq_mod(tidy_mod, name_cp_pcl, by.pcl, pcl_annot)
  sigs_cp <- merge(sigs_cp, tcga_genes, by = "pr_gene_symbol")
  
  sigs_cp <- sigs_cp %>%
          select(pr_gene_symbol, n_up, n_dn, Immune.Checkpoint) %>%
          group_by(pr_gene_symbol, Immune.Checkpoint) %>%
          summarise(su = sum(n_up), sdn = sum(n_dn))
   
  comp_plot <- ggplot(sigs_cp, aes(x = su, y = sdn, color = Immune.Checkpoint, label = pr_gene_symbol))+
    geom_point()+
    #geom_point(aes(gene_id = pr_gene_symbol))+
    geom_label_repel(data          = subset(sigs_cp, su > as.numeric(x_label) | sdn > as.numeric(y_label)),
                     nudge_y       = 10,
                     segment.size  = 0.2,
                     segment.color = "grey50")+ 
    scale_x_continuous(limits = x_scale)+
    scale_y_continuous(limits = y_scale)+
    xlab("# of signatures up-modulating")+
    ylab("# of signatures down-modulating")+
    ggtitle(paste("Number of signatures that up and down modulate a gene -", name_cp_pcl))+
    theme_light()
  
  print(comp_plot)
}


#plots a heatmap of up and down modulation
#by.cell = FALSE - does not segregate data across different cell lines, only different compounds

plot_heatmap_mod <- function(tidy_mod, name_cp_pcl, by.pcl = FALSE, by.cell = FALSE, pcl_annot = FALSE, filename = "heatmap.png", width = 12) {


  #Preparing the table for plotting
  
  freqs <- make_freq_mod(tidy_mod, name_cp_pcl, by.pcl, pcl_annot)
  
  
 if (by.cell == FALSE){
   #collapsing the information of frequencies for each gene
   
   freqs_gene <- freqs %>%
     select(pr_gene_symbol, pert_iname,total, n_up, n_dn) %>%
     group_by(pr_gene_symbol, pert_iname) %>%
     summarise(st = sum(total), su = sum(n_up), sdn = sum(n_dn)) %>% 
     mutate(diff = (su - sdn)/st)
     
   #This table is in the long format, but a heatmap requires it in the wide format
   freqs_wide <- freqs_gene %>% 
     select(pr_gene_symbol, pert_iname, diff) %>% 
     tidyr::spread(key = pr_gene_symbol, value = diff) %>% as.data.frame()
   
   freqs_wide <-  freqs_wide[,c(2:ncol(freqs_wide), 1)]
   
   rownames(freqs_wide) <- freqs_wide$pert_iname
 
   #organizing the annotation of perturbagen classes
   if(by.pcl != FALSE){
     
     pcl_df <- subset(pcl_annot, pert_iname %in% freqs_wide$pert_iname & pcl_id %in% name_cp_pcl, c(pert_iname, pcl_id)) %>% distinct()
     
     annot_df_rows <- freqs_wide %>% select(pert_iname) 
     
     annot_df_rows <- merge(annot_df_rows, pcl_df, by = "pert_iname")
     rnames <- (annot_df_rows$pert_iname)
     rownames(annot_df_rows) <- rnames
     
     annot_df_rows$pert_iname <- NULL
     
   } else { #only annotation of the cell lines and specific compounds
     annot_df_rows <- freqs_wide %>% 
       select(pert_iname)
     
     rnames <- (annot_df_rows$pert_iname)
     rownames(annot_df_rows) <- rnames
   }
   
 }else{ #you want to have a row for each combination of cell line and perturbation
   
  #This table is in the long format, but a heatmap requires it in the wide format
  freqs_wide <- freqs %>% 
    select(cell_id, pr_gene_symbol, pert_iname, diff) %>% 
    tidyr::spread(key = pr_gene_symbol, value = diff) 
  
  freqs_wide <-  freqs_wide[,c(3:ncol(freqs_wide), 1,2)]
  rnames <- paste(freqs_wide$cell_id, freqs_wide$pert_iname)
  rownames(freqs_wide) <- rnames
  

  #organizing the annotation of perturbagen classes
  if(by.pcl != FALSE){
    
    pcl_df <- subset(pcl_annot, pert_iname %in% freqs_wide$pert_iname & pcl_id %in% name_cp_pcl, c(pert_iname, pcl_id)) %>% distinct()
    
    annot_df_rows <- freqs_wide %>% 
      select(cell_id, pert_iname) 
    
    annot_df_rows <- merge(annot_df_rows, pcl_df, by = "pert_iname")
    rnames <- paste(annot_df_rows$cell_id, annot_df_rows$pert_iname)
    rownames(annot_df_rows) <- rnames
    
    annot_df_rows$pert_iname <- NULL
        
  } else { #only annotation of the cell lines and specific compounds
    annot_df_rows <- freqs_wide %>% 
      select(cell_id, pert_iname)
    
    rnames <- paste(annot_df_rows$cell_id, annot_df_rows$pert_iname)
    rownames(annot_df_rows) <- rnames
  }
 }

#organizing the annotation of gene classification
  
  #annotation of category of genes
  annot_df <- subset(tcga_genes, pr_gene_symbol %in% freqs$pr_gene_symbol, c(pr_gene_symbol, Super.Category, Immune.Checkpoint))
  rownames(annot_df) <-  annot_df$pr_gene_symbol
  annot_df$pr_gene_symbol <-  NULL

  #Now everything is ready for plotting the heatmap  
  
  p <- pheatmap(freqs_wide[,1:65], annotation_row = annot_df_rows, annotation_col = annot_df, angle_col = 45, filename = filename, width = width)
  return(pheatmap(freqs_wide[,1:65], annotation_row = annot_df_rows, annotation_col = annot_df, fontsize = 5, angle_col = 45, width = width))
  
}

#------Plots of distributions of z-scores------#
#This function DOES NOT take into consideration if a gene is up or down modulated. It will plot the distribution of z-scores
#of all signatures in a given combination of compound, cell line and gene.
#The function outputs a file with the boxplots and also print the plot.

plot_z_by_pcl <- function(gene_symbol, tidy_z, my_pcl, pcl_annot, cell_display = FALSE, cell_labels = FALSE, pcl_labels = FALSE){
  
  #selecting associated compounds
  pert_list <- subset(pcl_annot, pcl_id %in% my_pcl, c(pert_iname, pcl_id))
  
  gene_id <- subset(tcga_genes, pr_gene_symbol == gene_symbol, pr_gene_id)
  
  #organizing the list of cell lines that are present in the tidy_z df
  celllines <- subset(siginfo, sig_id %in% tidy_z$sig_id, cell_id) %>% distinct()
  
  ##control data for the selected genes and cell lines
  sigs_control <- list_ctl %>% filter(cell_id %in% celllines$cell_id) %>% select(sig_id)
  ctl_z <- get_zscores(gene_id, sigs_control$sig_id)
  
  
  #organizing the table for plotting
  ctl <- ctl_z %>% make_full_table %>% select(sig_id, cell_id, pr_gene_id, pert_iname, z_score)
  ctl$pcl_id <- "control"
  
  # pert_list <- stats_pert %>% group_by(cell_id, pr_gene_id) %>% select(cell_id, pcl_id, pr_gene_id, pert_iname, mean)
  
  #including the compound name in the modulation table
  pert_info <- siginfo %>% 
    select(sig_id, pert_iname, cell_id)
  
  tidy_z_full <- merge(tidy_z, pert_info, by = "sig_id")
  tidy_z_full <- merge(tidy_z_full, pcl_annot, by = "pert_iname")
  
  pert <- tidy_z_full %>%  filter(pr_gene_id %in% gene_id & pert_iname %in% pert_list$pert_iname & pcl_id %in% my_pcl) %>% 
    select(sig_id, cell_id, pcl_id, pr_gene_id, pert_iname, z_score) %>% distinct(sig_id, .keep_all = TRUE)
  
  plot_table <- rbind(as.data.frame(ctl), as.data.frame(pert))
  
  if(cell_labels !=FALSE){
    labels <- cell_labels
  }else{
    labels <- (plot_table$cell_id)
  }
  
  if(pcl_labels !=FALSE){
    labels_pcl <- pcl_labels
  }else{
    labels_pcl <- (plot_table$pcl_id)
  }
  
  if(cell_display != FALSE){
    #Plots are going to be divided by cell line and ordered in the provided order
    plot_table$cell_idf <- factor(plot_table$cell_id, cell_display)
    
    p <- ggplot(plot_table)+
      geom_boxplot(aes(x = factor(pert_iname), y=as.numeric(z_score)))+
      geom_jitter(aes(x = factor(pert_iname), y = as.numeric(z_score)))+
      facet_grid(cell_idf ~ pcl_id, scales = "free_x", space = "free",labeller = labeller(cell_idf = labels, pcl_id = labels_pcl))+
      geom_hline(yintercept = 0, color = "red", linetype = "dotted")+
      xlab("Compound name")+
      ylab("z-score")+
      theme_bw()+
      ggtitle(as.character(gene_symbol))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ggsave(paste(gene_symbol, ".png"), width = 15, height = 15)
    
  }else{
    p <- ggplot(plot_table)+
      geom_boxplot(aes(x = factor(pert_iname), y=as.numeric(z_score)))+
      geom_jitter(aes(x = factor(pert_iname), y = as.numeric(z_score)))+
      facet_grid(~ pcl_id, scales = "free_x", labeller = labeller(pcl_id = labels_pcl))+
      geom_hline(yintercept = 0, color = "red", linetype = "dotted")+
      xlab("Compound name")+
      ylab("z-score")+
      theme_bw()+
      ggtitle(as.character(gene_symbol))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ggsave(paste(gene_symbol, ".png"), width = 15, height = 15) 
  }
  print(p)
}


