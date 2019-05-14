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

The ```load_data.R``` file takes care of loading the data for users at the Institute of Systems Biology (ISB). For these users, the necessary files are available at the ```smb://bassonn/CancerRegulome14``` server. If you are connected to this server, there is no need to download the data.

One can set a different path to the data by setting the variable ```CMAP_HOME``` to the desired location. The present default is the folder at the ISB server.
