---
title: "BinfTools"
output: html_document
bibliography: bt.bib
---

BinfTools is a package designed to help beginners with bioinformatics analyses - primarily with RNA-seq analyses. The idea behind this package is that differential expression analysis has already been performed using DESeq2 [@DESeq2] and we are then ready to run some analyses and make some figures.

# Installing BinfTools

Currently, BinfTools is only available for installation from GitHub:

```{r}
devtools::install_github("kevincjnixon/BinfTools")
library(BinfTools)
```
See the DESCRIPTION for the list of R dependencies. If these are already installed, installation should be very quick.

# Usage

Ideally, we will be using objects created from DESeq2, but objects from other RNA-seq analysis packages can be modified to work as input for BinfTools.

```{r}
###Example DESeq2 pipeline###
library(DESeq2)
#Create a sample table with sample names, filenames of count tables from HTSeq-Count, and conditions
sampleTable<-data.frame(Sample=c("S1","S2","S3","S4","S5","S6"),
  fileNames=c("S1.txt","S2.txt","S3.txt","S4.txt","S5.txt","S6.txt"),
  Condition=as.factor(c("WT","WT","WT","KO","KO","KO")))
#Create the DESeq data object
dds<-DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~Condition)
#Run DESeq
dds<-DESeq(dds)
#Get the results
res<-results(dds, contrasts=c("Condition","KO","WT"))
```

## Available Functions

### GenerateGSEA

This function will generate a ranked list of genes from the *res* object to run a PreRanked gene set enrichment analysis GSEA [@gsea05]
See *?GenerateGSEA* to view all options

```{r}
GenerateGSEA(res, fileName="GSEA.rnk", bystat=T, byFC=F)
```
### MA_plot()

This function is similar to the *plotMA()* function from DESeq2 [@DESeq2], however it is more customizable.
See *?MA_plot* to view all options

```{r}
#Pull top significant genes to label and colour
genes<-rownames(res[order(res$padj),c(1:5)])
#Make an MA plot
MA_plot(res, title="MA Plot", p=0.05, FC=1, lab=genes, col=genes)
```
### volcanoPlot()

This function will generate a volcano plot from the *res* object.
See *?volcanoPlot* to view all options

```{r}
volcanoPlot(res, title="Volcano Plot", p=0.05, FC=2, lab=genes, col=genes)
```

### GO_GEM()

This function will run a gene ontology analysis using gprofiler2's *gost()* function [@gprofiler2] and output the results in a table (.txt), a .gem file (for compatibility with Cytoscape's EnrichmentMap app [@em10; @em19]), and two figures in a single .pdf displaying the top 10 enriched and top 10 significant terms.
See *?GO_GEM* to view all options

```{r}
#Get upregulated gene names
genes<-rownames(subset(res, padj<0.05 & log2FoldChange > 1))
#Run GO analysis for Biological Process and output to folder "GO/GO_analysis*"
GO_GEM(genes, species="hsapiens", bg=rownames(res), source="GO:BP", prefix="GO/GO_analysis")
```

### zheat()

This function will generate a heatmap of z-score normalized gene counts using pheatmap [@pheat].
See *?zheat* to view all options

```{r}
#Get the names of all differentially expressed genes
genes<-rownames(subset(res, padj<0.05 & abs(log2FoldChange)>1))
#Generate a heatmap
zheat(genes=genes, counts=counts(dds, normalized=T), conditions=dds$Condition, con="WT", title="DEGs")
```

### count_plot()

This function Will generate a violin plot of normalized gene expression for a specified group of genes.
See *?count_plot* to view all options.

```{r}
#Check all upregulated genes
genes<-rownames(subset(res, padj<0.05 & log2FoldChange >1))
#Compare gene expression between the two conditions using z-score normalized counts:
count_plot(counts=counts(dds, normalized=T), scaling="zscore", genes=genes, condition=dds$Condition, title="Upregulated Genes", compare=list(c("WT","KO")))
```

### gsva_plot()

This function will run a single sample gene set enrichment analysis (ssGSEA) using the GSVA package [@GSVA] and generate a violin plot of normalized enrichment scores using custom gene sets. This is especially useful after running a GSEA using the rnk file generated using *GenerateGSEA()* and finding a cluster of related pathways in the EnrichmentMap [see @em19].
See *?gsva_plot* to view all options.

```{r}
#Import the gene sets associated with pathways of an enrichment map cluster
library(qusage)
geneSet<-read.gmt("pathway_geneSets.gmt")

#Now run the gsva using z-score normalized gene expression and make a plot:
gsva_plot(counts=t(scale(t(counts(dds, normalized=T)))), geneset=geneSet, method="ssgsea", condition=dds$Condition, title="ssGSEA", compare=list(c("WT","KO")))

```
## Other Dependencies Not Mentioned:
- ggplot2 [@ggplot2]
- ggpubr [@ggpub]
- rstatix [@rstatix]
- tidyr [@tidyr]
- dplyr [@dplyr]

# References
