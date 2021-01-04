---
title: "Introduction"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: bt.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=8,
  fig.height=8,
  fig.align = "center"
)
```

```{r setup}
#devtools::install_github("kevincjnixon/BinfTools")
library(BinfTools)
```
# Usage

Ideally, we will be using objects created from DESeq2, but objects from other RNA-seq analysis packages can be modified to work as input for BinfTools.
Here, we are using RNA-seq data from Nixon et al. [-@nixon] comparing gene expression in the *Drosophila* musrhoom body between Wild-Type (expressing mCherry-shRNA) and Bap60 knockdown.

```{r data}
###Load data generated from DESeq2:
head(norm_counts) #Normalized gene counts (from counts(dds, normalized=T))
head(res) #Results object (from as.data.frame(results(dds, contrast=c("Condition","WT","KO"))))
print(cond) #Character vector of conditions (from as.character(dds$Condition))
print(geneSet) #List of gene sets related to rhodopsin signaling imported from qusage::read.gmt("rhodopsin.gmt")
```

## Available Functions

### Converting results from Limma or EdgeR

If you prefer to use Limma or EdgeR for analysis, the functions *fromLimma()* or *fromEdgeR()* will convert the output of Limma's *topTable()* function or EdgeR's *topTags()* function to data frames compatible with BinfTools functions.

```{r converting}
### Converting Limma results:
#res<-limma::topTable(fit, number=Inf)
#res<-fromLimma(res)
###Converting EdgeR results:
#res<-edgeR::topTags(de, n=Inf)
#res<-fromEdgeR(res)
```

### GenerateGSEA

This function will generate a ranked list of genes from the *res* object to run a PreRanked gene set enrichment analysis GSEA [@gsea05]
See *?GenerateGSEA* to view all options

```{r GSEA, warning=FALSE}
GenerateGSEA(res, filename="GSEA.rnk", bystat=T, byFC=F)
```
### MA_Plot()

This function is similar to the *plotMA()* function from DESeq2 [@DESeq2], however it is more customizable.
See *?MA_Plot* to view all options

```{r MA_Plot, fig.cap="MA Plot", message=FALSE}
#Pull top significant genes to label and colour
genes<-rownames(res[order(res$padj),])[1:5]
#Make an MA plot
MA_Plot(res, title="MA Plot", p=0.05, FC=log(1.25,2), lab=genes, col=genes)
```
### volcanoPlot()

This function will generate a volcano plot from the *res* object.
See *?volcanoPlot* to view all options

```{r volcano, fig.cap="Volcano Plot diplaying upregulated (red) and downregulated (blue) genes", message=FALSE}
volcanoPlot(res, title="Volcano Plot", p=0.05, FC=log(1.25,2), lab=genes, col=genes)
```

### GO_GEM()

This function will run a gene ontology analysis using gprofiler2's *gost()* function [@gprofiler2] and output the results in a table (.txt), a .gem file (for compatibility with Cytoscape's EnrichmentMap app [@em10; @em19]), and two figures in a single .pdf displaying the top 10 enriched and top 10 significant terms.
See *?GO_GEM* to view all options

```{r GO, warning=FALSE, fig.cap=c("Top 10 Enriched GO terms - minimum 10 genes/term","Top 10 Significant GO terms - maximum 500 genes/term")}
#Get downregulated gene names
genes<-rownames(subset(res, padj<0.05 & log2FoldChange < log(1.25, 2)))
#create output directory
dir.create("GO")
#Run GO analysis for Biological Process and output to folder "GO/GO_analysis*"
GO_GEM(genes, species="dmelanogaster", bg=rownames(res), source="GO:BP", prefix="GO/GO_analysis")
```

### zheat()

This function will generate a heatmap of z-score normalized gene counts using pheatmap [@pheat].
See *?zheat* to view all options

```{r heatmap, fig.cap="Heatmap of DEGs"}
#Get the names of all differentially expressed genes
DEGs<-rownames(subset(res, padj<0.05 & abs(log2FoldChange)>log(1.25,2)))
#Pull top significant genes to label and colour
top_genes<-rownames(res[order(res$padj),])[1:5]
#Generate a heatmap
zheat(genes=DEGs, counts=norm_counts, conditions=cond, con="WT", title="DEGs", labgenes = top_genes)
```

### count_plot()

This function Will generate a violin plot of normalized gene expression for a specified group of genes.
See *?count_plot* to view all options.

```{r count_plot, fig.cap="Violin plot comparing gene expression of upregulated genes between conditions"}
#Check all upregulated genes
genes<-rownames(subset(res, padj<0.05 & log2FoldChange > log(1.25,2)))
#Compare gene expression between the two conditions using z-score normalized counts:
count_plot(counts=norm_counts, scaling="zscore", genes=genes, condition=cond, title="Upregulated Genes", compare=list(c("WT","KO")))
```

### gsva_plot()

This function will run a single sample gene set enrichment analysis (ssGSEA) using the GSVA package [@GSVA] and generate a violin plot of normalized enrichment scores using custom gene sets. This is especially useful after running a GSEA using the rnk file generated using *GenerateGSEA()* and finding a cluster of related pathways in the EnrichmentMap [see @em19].
See *?gsva_plot* to view all options.

```{r gsva, fig.cap="Violin plot comparing normalized enrichment scores of pathways involved in rhodopsin signaling between conditions"}
#Import the gene sets associated with pathways of an enrichment map cluster - In this case rhodopsin signaling
#library(qusage)
#geneSet<-read.gmt("rhodopsin.gmt")

#Now run the gsva using z-score normalized gene expression and make a plot:
gsva_plot(counts=t(scale(t(norm_counts))), geneset=geneSet, method="ssgsea", condition=cond, title="Rhodopsin-Mediated Signaling", compare=list(c("WT","KO")))

```

## Other Dependencies Not Mentioned:
- ggplot2 [@ggplot2]
- ggpubr [@ggpub]
- rstatix [@rstatix]
- tidyr [@tidyr]
- dplyr [@dplyr]

# References
