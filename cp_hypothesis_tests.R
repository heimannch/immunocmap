#Module with functions for campare gene expression distributions between controls and treatment

source('cp_plots_pert.R')

#takes as input the gene_id of the gene and the cp_name, the name of the treatment to compare with the control

gene_ttest <-  function(mod_df, ctl_df, gene_id, cp_name){
  
  ctl_gene <- ctl_df %>% filter(pr_gene_id == gene_id)
  pert_gene <- mod_df %>% filter(pr_gene_id == gene_id & pert_iname == cp_name) #tirar dependencia de mod_list_full
  
  ttest <- (t.test(pert_gene$z_score, ctl_gene$z_score, alternative = "two.sided", paired = FALSE))
  
  
  return(ttest)
}


#In following sections, we will implement functions to do a batch of t-tests at once. 
#In order to better manipulate the results, it is simpler to convert the output to a dataframe

df_ttest <- function(list_ttest, all_cp){
  
  all_test_df <- (rbindlist(list_ttest))
  
  #the dataframe has repeated rows, so we delete them
  toDelete <- seq(1, nrow(all_test_df), 2)
  test_df <- all_test_df[toDelete,]
  
  test_df <- cbind(all_cp, test_df)
  
  return(test_df)
}


#The following function allows the user to run the t-test for a given gene across all drugs that are available in a give list of sig_id and modulation
test_per_gene <- function(gene, mod_df, ctl_df){
  
  all_cp <-as.data.frame(unique(mod_df$pert_iname)) 
  all_test <- apply(all_cp, 1, gene_ttest, mod_df = mod_df, ctl_df = ctl_df, gene_id = gene)
  
  all_test_df <- df_ttest(all_test, all_cp)
  all_test_df$pr_gene_id <- as.numeric(gene)

  return(all_test_df)
}

test_all <- function(gene_list, mod_df, ctl_df){
  all_tests <- data.frame()
  for (i in 1:length(gene_list)) {
    new_test <- test_per_gene(gene_list[i], mod_df, ctl_df)
    new_test$pr_gene_id <- gene_list[i]
    
    all_tests <- rbind(all_tests, new_test)
    
    print(i)
    i=i+1
  }
  colnames(all_tests)[1] <- "pert_iname"
  return(all_tests)
}

#-------FUNCTIONS TO PLOT HEATMAPS WITH THE RESULT OF TESTS------------------


test_wide <- function(tests_df){
  
  ph_df <- tests_df %>% select(pert_iname, pr_gene_id, statistic, p.value, conf.int)
  
  #taking the log of p-value, with different signs for positive or negative t stats
  
  ph_df <- ph_df %>% dplyr::mutate(logp = case_when( statistic > 0 ~ -log10(p.value),
                                              statistic < 0 ~ log10(p.value))
  )
  
  #organizing the data for plotting
  ph_wide <- ph_df %>% 
    select(pert_iname, pr_gene_id, logp) %>% 
    tidyr::spread(pert_iname, logp)
  
  ph_wide <- merge(tcga_genes[,1:2], ph_wide, by = "pr_gene_id")
  rownames(ph_wide) <- ph_wide$pr_gene_symbol
  ph_wide$pr_gene_id <- NULL
  ph_wide$pr_gene_symbol <- NULL
  
  return(ph_wide)
  
}

plot_heatmap_test <- function(test_tidy, wide_table, clustering_distance_rows = "euclidean",
                              clustering_distance_cols = "euclidean", filename = "heatmap.png", height = 25, width = 15){
  
  #annotation of category of genes
  annot_df <- subset(tcga_genes, pr_gene_symbol %in% rownames(wide_table), c(pr_gene_symbol, Super.Category, Immune.Checkpoint))
  rownames(annot_df) <-  annot_df$pr_gene_symbol
  annot_df$pr_gene_symbol <-  NULL
  
  #annotation of classes of drugs
  pcl_df <- merge(pcl_custom[,2:3], test_tidy[1], by= "pert_iname") %>% distinct()
  pcl_df <- subset(pcl_df, !duplicated(pcl_df [,1]) )
  rownames(pcl_df) <- pcl_df$pert_iname
  pcl_df$pert_iname <- NULL
  
  # Making sure that 0 is in the center of the color scheme
  breaks_list <-  seq(-10, 10, by = 0.5)
  my_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaks_list))
  
  #generating the heatmap
  wide_table %>% t() %>% 
    pheatmap::pheatmap(color = my_colors, breaks = breaks_list, clustering_distance_rows = clustering_distance_rows, clustering_distance_cols = clustering_distance_cols, 
                       annotation_row = pcl_df, annotation_col = annot_df, annotation_legend = FALSE)
  
   wide_table %>% t() %>% 
    pheatmap::pheatmap(color = my_colors, breaks = breaks_list, clustering_distance_rows = clustering_distance_rows, clustering_distance_cols = clustering_distance_cols, 
                       annotation_row = pcl_df, annotation_col = annot_df, filename = filename, height = height, width = width)
  
  
  
}
