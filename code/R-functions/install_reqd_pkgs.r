# Author: Darren Tyson
# Date: 2022-10-11
# Code for installing all necessary packages to run all code in
# this repository; code for Diamantopoulou et al. 2022
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

req_pkgs <- c("tidyverse","scater","scran","edgeR","clusterProfiler","GSVA","foreach",
              "showtext","circlize","RColorBrewer","DT","cowplot","ggbeeswarm","ggpubr",
              "ggrepel","clusterProfiler","enrichplot","ComplexHeatmap","limma","colorblindr",
              "grid","gridExtra","ggalt","kableExtra","PCAtools")

inst_pkgs <- installed.packages()[,"Package"]

needed_pkgs <- setdiff(req_pkgs,inst_pkgs)
std_pkgs <- needed_pkgs[needed_pkgs %in% tools::CRAN_package_db()[,"Package"]]
bioc_pkgs <- needed_pkgs[needed_pkgs %in% BiocManager::available()]
other_pkgs <- setdiff(needed_pkgs,c(std_pkgs,bioc_pkgs))

if(!length(other_pkgs)==0) message(cat("Packages not found in CRAN or BioC still needed:",paste(other_pkgs, collapse=", ")))

if(length(std_pkgs)==0 & length(bioc_pkgs)==0) {
  message('All required CRAN and BioC packages have previously been installed')
} else {
    if(!is.null(std_pkgs)) install.packages(std_pkgs)
    if(!is.null(bioc_pkgs)) BiocManager::install(bioc_pkgs)
}

if("colorblindr" %in% needed_pkgs) devtools::install_github("clauswilke/colorblindr")

rm(list=c("req_pkgs","inst_pkgs","needed_pkgs","std_pkgs","bioc_pkgs","other_pkgs"))
