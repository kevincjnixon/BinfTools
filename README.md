# BinfTools
Useful R functions for Multi-omics analyses

BinfTools - Quick Start Guide
=============================
View the [BinfTools Vignette](https://drive.google.com/file/d/1XM9iBdcQBpgcXPtu0zAVzUyogzuk7DnS/view?usp=sharing).

For more detailed information about all of the functions, see the [BinfTools Wiki](https://github.com/kevincjnixon/BinfTools/wiki)

Installation
============
```
#Bioconductor packages don't install automatically on BinfTools install
BiocManager::install("SAGx") #This package isn't available for the latest version of R. It has been removed as a dependency
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
