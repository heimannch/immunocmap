#Module to load all the data in.
#!/usr/bin/env Rscript

library(data.table)


#User should set the UNIX path to files by running in the command line: 'export CMAP_HOME=path/to/files'


folder_path <- Sys.getenv("CMAP_HOME")


siginfo_path <- paste(folder_path, "/GSE92742_Broad_LINCS_sig_info.txt", sep = "")
sig_metrics_path <- paste(folder_path,"/GSE92742_Broad_LINCS_sig_metrics.txt.gz", sep = "")
cellinfo_path <-  paste(folder_path,"/GSE92742_Broad_LINCS_cell_info.txt", sep = "")
geneinfo_path <- paste(folder_path,"/GSE92742_Broad_LINCS_gene_info.txt.gz", sep = "")
geneinfo_lm_path <- paste(folder_path, "/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz", sep = "")
pertinfo_path <- paste(folder_path, "/GSE92742_Broad_LINCS_pert_info.txt", sep = "")
pertmetrics_path <- paste(folder_path,"/GSE92742_Broad_LINCS_pert_metrics.txt.gz", sep = "")
instinfo_path <- paste(folder_path,"/GSE92742_Broad_LINCS_inst_info.txt.gz", sep = "")
ds_path <- paste(folder_path,"/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", sep = "")

#Files saved locally
pcl_path <- "data/pcls.csv"
pcl_custom_path <- "data/pclsCustom.csv"
tcga_path <- "data/tcga_at_cmap.csv"


#Loading the files

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

##This last table contains info about all instances of the experiments, which is only necessary when dealing with levels 3 or below of the data.
#instinfo <- suppressWarnings(data.table::fread(instinfo_path))
