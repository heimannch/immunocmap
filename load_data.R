#Module to load all the data in. Needs to be connected with the bassoon/cancerregulome server at ISB.

library(data.table)


#ds_path <- "data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
pcl_path <- "data/pcls.csv"
pcl_custom_path <- "data/pclsCustom.csv"
tcga_path <- "data/tcga_at_cmap.csv"
# 
# siginfo_path <- "data/GSE92742_Broad_LINCS_sig_info.txt"
# sig_metrics_path <- "data/GSE92742_Broad_LINCS_sig_metrics.txt"
# cellinfo_path <-  "data/GSE92742_Broad_LINCS_cell_info.txt"
# geneinfo_path <- "data/GSE92742_Broad_LINCS_gene_info.txt"
# geneinfo_lm_path <- "data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt"
# pertinfo_path <- "data/GSE92742_Broad_LINCS_pert_info.txt"
# pertmetrics_path <- "data/GSE92742_Broad_LINCS_pert_metrics.txt"
#instinfo_path <- "data/GSE92742_Broad_LINCS_inst_info.txt"
 
#At ISB
siginfo_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_sig_info.txt.gz"
sig_metrics_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_sig_metrics.txt.gz"
cellinfo_path <-  "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_cell_info.txt.gz"
geneinfo_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_gene_info.txt.gz"
geneinfo_lm_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz"
pertinfo_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_pert_info.txt.gz"
pertmetrics_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_pert_metrics.txt.gz"
instinfo_path <- "/Volumes/CancerRegulome14/users/gqin/CMap/PhI_GSE92742/GSE92742_Broad_LINCS_inst_info.txt.gz"
ds_path <- "Volumes/CancerRegulome14/users/gqin/CMap/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"

#   

siginfo <- suppressWarnings(data.table::fread(siginfo_path))
sigmetrics <- suppressWarnings(data.table::fread(sig_metrics_path))
cellinfo <- suppressWarnings(data.table::fread(cellinfo_path))
geneinfo <- suppressWarnings(data.table::fread(geneinfo_path))
geneinfo_lm <- suppressWarnings(data.table::fread(geneinfo_lm_path))
pertinfo <- suppressWarnings(data.table::fread(pertinfo_path))
pertmetrics <- suppressWarnings(data.table::fread(pertmetrics_path))

pcl_tidy <- read.csv(pcl_path , header = TRUE, row.names = 1)
pcl_custom <- read.csv(pcl_custom_path , header = TRUE, row.names = 1, stringsAsFactors = FALSE)
tcga_genes <- read.csv(tcga_path, header = TRUE, sep = "\t")
#instinfo <- suppressWarnings(data.table::fread(instinfo_path))

# 
# pertinfo2 <- suppressWarnings(data.table::fread("data/GSE70138_Broad_LINCS_pert_info.txt.gz"))
# siginfo2 <- suppressWarnings(data.table::fread("data/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"))
# cellinfo2 <- suppressWarnings(data.table::fread("data/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz"))
