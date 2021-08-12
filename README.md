# BinfTools
Useful R functions for Multi-omics analyses

BinfTools - Quick Start Guide
=============================

For more detailed information about all of the functions, see the [BinfTools Wiki](https://github.com/kevincjnixon/BinfTools/wiki)

Installation
============
```
#Bioconductor packages don't install automatically on BinfTools install
BiocManager::install("SAGx")
BiocManager::install("GSVA")
BiocManager::install("fgsea")
BiocManager::install("gage")

devtools::install_github("kevincjnixon/gpGeneSets") # gprofiler genesets for mouse, human, and drosophila
devtools::install_github("kevincjnixon/BinfTools", build_vignettes=T) #This will take some time to build the vignette
```
Load the libraries
```
library(gpGeneSets)
library(BinfTools)
vignette("BinfTools") #Quick getting started vignette
```
