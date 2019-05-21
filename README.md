# Connectivity Map in the context of Immune Response

Code for analysis of immunomodulator genes in the Connectivity Map dataset

For a tutorial, please refer to the file ```cp_analysis_TNBC_cells.Rmd```. 

The code is developed entirely in **R** using the **R Studio** development environment. In addition, the following packages are used:
- ```cmapR```
- ```dplyr```
- ```data.table```
- ```ggplot2```
- ```ggrepel```
- ```pheatmap```
- ```plyr```
- ```skimr```

The ```cmapR``` package is not available in CRAN as for this moment. To install the development version from GitHub:
```devtools::install_github("cmap/cmapR")```

# Data

All the data used for this project is available for download at [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742). 
The files needed in this project are:
```GSE92742_Broad_LINCS_sig_info.txt```
```GSE92742_Broad_LINCS_sig_metrics.txt.gz```
```GSE92742_Broad_LINCS_cell_info.txt```
```GSE92742_Broad_LINCS_gene_info.txt.gz```
```GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz```
```GSE92742_Broad_LINCS_pert_info.txt```
```GSE92742_Broad_LINCS_pert_metrics.txt.gz```
```GSE92742_Broad_LINCS_inst_info.txt.gz```
```GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx```


The ```load_data.R``` file takes care of loading the data. Users need to set the path to the downloaded data by setting the variable ```CMAP_HOME``` to the desired location.
For example, in unix bash shell, this can be set by running: ```export CMAP_HOME=path/to/files```
In a R Session, this can be set with the command ```Sys.setenv(CMAP_HOME="path/to/file")```